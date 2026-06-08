use std::io::{self, BufRead};
use std::marker::PhantomData;

use vita_core::tensor::Point3;
use vita_core::units::length::{Angstrom, Length};
use vita_core::{Element, Scalar};

use super::error::{Error, ErrorKind, ParseError};
use super::system::System;
use crate::Location;

/// Begins reading XYZ data from `reader`.
///
/// The reader is consumed: take a single structure with [`Reader::system`], or stream a
/// trajectory frame by frame with [`Reader::systems`]. To keep ownership of the reader,
/// pass a mutable reference (`&mut reader` is itself a [`BufRead`]).
#[inline]
pub fn read<R: BufRead>(reader: R) -> Reader<R> {
    Reader {
        lines: Lines::new(reader),
    }
}

/// A reader for XYZ data, returned by [`read`].
pub struct Reader<R> {
    lines: Lines<R>,
}

impl<R: BufRead> Reader<R> {
    /// Reads the first system from the input.
    ///
    /// Any frames beyond the first are ignored; use [`systems`](Reader::systems) to read
    /// a trajectory.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`](crate::Error::Parse) if the input is empty or the first
    /// frame is malformed, or [`Error::Io`](crate::Error::Io) if the reader fails.
    pub fn system<V: Scalar>(mut self) -> Result<System<V>, Error> {
        match parse_system(&mut self.lines)? {
            Some(system) => Ok(system),
            None => Err(text_error(
                self.lines.number.saturating_add(1),
                None,
                ErrorKind::UnexpectedEof,
            )),
        }
    }

    /// Reads the input as a trajectory: a stream of systems, one per frame.
    ///
    /// The returned iterator yields each frame lazily, and stops at the end of the input
    /// or at the first error.
    #[inline]
    pub fn systems<V: Scalar>(self) -> Systems<V, R> {
        Systems {
            lines: self.lines,
            done: false,
            marker: PhantomData,
        }
    }
}

/// An iterator over the systems in an XYZ trajectory, returned by [`Reader::systems`].
///
/// Each [`next`](Iterator::next) parses one frame. Iteration ends after the last frame,
/// or after the first error.
pub struct Systems<V, R> {
    lines: Lines<R>,
    done: bool,
    marker: PhantomData<fn() -> V>,
}

impl<V: Scalar, R: BufRead> Iterator for Systems<V, R> {
    type Item = Result<System<V>, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match parse_system(&mut self.lines) {
            Ok(Some(system)) => Some(Ok(system)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(error) => {
                self.done = true;
                Some(Err(error))
            }
        }
    }
}

impl<V: Scalar, R: BufRead> std::iter::FusedIterator for Systems<V, R> {}

/// Parses one system, or `Ok(None)` at a clean end of input (no further systems).
fn parse_system<V: Scalar, R: BufRead>(lines: &mut Lines<R>) -> Result<Option<System<V>>, Error> {
    let count = loop {
        match lines.next()? {
            None => return Ok(None),
            Some((_, line)) if line.trim().is_empty() => {}
            Some((number, line)) => {
                break parse_count(line).map_err(|kind| text_error(number, None, kind))?;
            }
        }
    };

    let comment: Box<str> = match lines.next()? {
        Some((_, line)) => Box::from(line),
        None => {
            return Err(text_error(
                lines.number.saturating_add(1),
                None,
                ErrorKind::UnexpectedEof,
            ));
        }
    };

    let mut elements = Vec::with_capacity(count as usize);
    let mut positions = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let (number, line) = match lines.next()? {
            Some(pair) => pair,
            None => {
                return Err(text_error(
                    lines.number.saturating_add(1),
                    None,
                    ErrorKind::UnexpectedEof,
                ));
            }
        };
        let (element, position) = parse_atom::<V>(line, number)?;
        elements.push(element);
        positions.push(position);
    }

    Ok(Some(System::from_parts(
        comment,
        elements.into_boxed_slice(),
        positions.into_boxed_slice(),
    )))
}

/// Parses the atom-count line: the whole line, ignoring surrounding whitespace.
///
/// Parses first as `u64` to distinguish a syntax error (not a non-negative integer at all)
/// from a range error (valid integer but above vita-core's `SiteId` limit of `u32::MAX`).
fn parse_count(line: &str) -> Result<u32, ErrorKind> {
    let text = line.trim();
    let n = text
        .parse::<u64>()
        .map_err(|_| ErrorKind::AtomCount { found: text.into() })?;
    u32::try_from(n).map_err(|_| ErrorKind::AtomCountRange { count: n })
}

/// Parses one atom line into its element and position.
fn parse_atom<V: Scalar>(
    line: &str,
    number: u32,
) -> Result<(Element, Point3<Length<V, Angstrom>>), Error> {
    let mut fields = line.split_whitespace();
    let (Some(symbol), Some(x), Some(y), Some(z), None) = (
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
    ) else {
        return Err(text_error(
            number,
            None,
            ErrorKind::FieldCount {
                found: line.split_whitespace().count(),
            },
        ));
    };

    let element = Element::from_symbol(symbol).ok_or_else(|| {
        text_error(
            number,
            Some(column_of(line, symbol)),
            ErrorKind::ElementSymbol {
                found: symbol.into(),
            },
        )
    })?;
    let px = coordinate::<V>(line, x, number)?;
    let py = coordinate::<V>(line, y, number)?;
    let pz = coordinate::<V>(line, z, number)?;

    Ok((
        element,
        Point3::new(Length::new(px), Length::new(py), Length::new(pz)),
    ))
}

/// Parses one coordinate field, reporting its column on failure.
fn coordinate<V: Scalar>(line: &str, field: &str, number: u32) -> Result<V, Error> {
    parse_coordinate::<V>(field).ok_or_else(|| {
        text_error(
            number,
            Some(column_of(line, field)),
            ErrorKind::Coordinate {
                found: field.into(),
            },
        )
    })
}

/// Parses a coordinate token, accepting only a finite real number.
fn parse_coordinate<V: Scalar>(field: &str) -> Option<V> {
    field
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .map(V::from_f64)
}

/// Returns the 1-based byte column at which `field` begins within `line`.
///
/// `field` must be a sub-slice of `line`, as produced by splitting it.
fn column_of(line: &str, field: &str) -> u32 {
    (field.as_ptr() as usize - line.as_ptr() as usize) as u32 + 1
}

/// Builds a parse error at a text location.
fn text_error(line: u32, column: Option<u32>, kind: ErrorKind) -> Error {
    Error::Parse(ParseError::new(Location::Text { line, column }, kind))
}

/// A line reader that strips line terminators and tracks the 1-based line number.
struct Lines<R> {
    reader: R,
    buffer: String,
    number: u32,
}

impl<R: BufRead> Lines<R> {
    #[inline]
    fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: String::new(),
            number: 0,
        }
    }

    /// Reads the next line, returning its 1-based number and its content without the
    /// trailing line terminator, or `None` at end of input.
    fn next(&mut self) -> io::Result<Option<(u32, &str)>> {
        self.buffer.clear();
        if self.reader.read_line(&mut self.buffer)? == 0 {
            return Ok(None);
        }
        self.number = self.number.saturating_add(1);
        let mut line = self.buffer.as_str();
        if let Some(rest) = line.strip_suffix('\n') {
            line = rest;
        }
        if let Some(rest) = line.strip_suffix('\r') {
            line = rest;
        }
        Ok(Some((self.number, line)))
    }
}

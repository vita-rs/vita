use std::io::{self, BufRead};
use std::marker::PhantomData;

use vita_core::tensor::Point3;
use vita_core::units::length::{Angstrom, Length};
use vita_core::{Element, Scalar};

use super::{Error, ErrorKind, ParseError, System};
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

    let element = Element::from_symbol(&canonical_symbol(symbol)).ok_or_else(|| {
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

/// Converts `s` to canonical element-symbol case: initial cap, rest lowercase.
fn canonical_symbol(s: &str) -> String {
    let mut chars = s.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => {
            let mut out = String::with_capacity(s.len());
            out.push(first.to_ascii_uppercase());
            for c in chars {
                out.push(c.to_ascii_lowercase());
            }
            out
        }
    }
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
    let offset = field.as_ptr() as usize - line.as_ptr() as usize;
    offset.saturating_add(1).min(u32::MAX as usize) as u32
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

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{HasElements, HasPositions, HasSites, SiteId};

    const FRAME: &str = "2\nwater\nO 1 2 3\nH 4 5 6\n";

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    fn point(x: f64, y: f64, z: f64) -> Point3<Length<f64, Angstrom>> {
        Point3::new(Length::new(x), Length::new(y), Length::new(z))
    }

    fn parse(input: &str) -> Result<System, Error> {
        read(input.as_bytes()).system()
    }

    fn parse_all(input: &str) -> Vec<Result<System, Error>> {
        read(input.as_bytes()).systems().collect()
    }

    fn parse_error(input: &str) -> ParseError {
        match parse(input).unwrap_err() {
            Error::Parse(e) => e,
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    fn text(line: u32, column: Option<u32>) -> Location {
        Location::Text { line, column }
    }

    #[test]
    fn empty_input_is_an_unexpected_eof_error() {
        assert_eq!(
            parse_error(""),
            ParseError::new(text(1, None), ErrorKind::UnexpectedEof),
        );
    }

    #[test]
    fn a_zero_atom_frame_parses_to_an_empty_system() {
        let system = parse("0\ncomment\n").unwrap();
        assert_eq!(system.site_count(), 0);
        assert_eq!(system.comment(), "comment");
    }

    #[test]
    fn parses_the_declared_atom_count() {
        assert_eq!(parse(FRAME).unwrap().site_count(), 2);
    }

    #[test]
    fn parses_each_atom_element() {
        let system = parse(FRAME).unwrap();
        assert_eq!(system.element(s(1)), elem("O"));
        assert_eq!(system.element(s(2)), elem("H"));
    }

    #[test]
    fn parses_each_atom_position() {
        let system = parse(FRAME).unwrap();
        assert_eq!(system.position::<Angstrom>(s(1)), point(1.0, 2.0, 3.0));
        assert_eq!(system.position::<Angstrom>(s(2)), point(4.0, 5.0, 6.0));
    }

    #[test]
    fn keeps_the_comment_verbatim() {
        let system = parse("1\n  spaced comment  \nH 0 0 0\n").unwrap();
        assert_eq!(system.comment(), "  spaced comment  ");
    }

    #[test]
    fn keeps_an_empty_comment() {
        assert_eq!(parse("1\n\nH 0 0 0\n").unwrap().comment(), "");
    }

    #[test]
    fn treats_element_symbols_case_insensitively() {
        assert_eq!(parse("1\nc\nhe 0 0 0\n").unwrap().element(s(1)), elem("He"));
        assert_eq!(parse("1\nc\nFE 0 0 0\n").unwrap().element(s(1)), elem("Fe"));
    }

    #[test]
    fn accepts_arbitrary_whitespace_between_fields() {
        let system = parse("1\nc\nH  1\t2   3\n").unwrap();
        assert_eq!(system.position::<Angstrom>(s(1)), point(1.0, 2.0, 3.0));
    }

    #[test]
    fn accepts_surrounding_whitespace_around_the_count() {
        assert_eq!(
            parse("  2  \nc\nH 0 0 0\nO 1 1 1\n").unwrap().site_count(),
            2
        );
    }

    #[test]
    fn accepts_scientific_notation_for_coordinates() {
        let system = parse("1\nc\nH 1e3 -2.5 1.5E-2\n").unwrap();
        assert_eq!(
            system.position::<Angstrom>(s(1)),
            point(1000.0, -2.5, 0.015)
        );
    }

    #[test]
    fn reads_coordinates_into_the_f32_scalar() {
        let system = read("1\nc\nH 1.5 2.5 3.5\n".as_bytes())
            .system::<f32>()
            .unwrap();
        assert_eq!(
            system.position::<Angstrom>(s(1)),
            Point3::new(Length::new(1.5_f32), Length::new(2.5), Length::new(3.5)),
        );
    }

    #[test]
    fn accepts_carriage_return_line_endings() {
        let system = parse("1\r\ncomment\r\nH 0 0 0\r\n").unwrap();
        assert_eq!(system.comment(), "comment");
    }

    #[test]
    fn skips_blank_lines_before_a_frame() {
        assert_eq!(parse("\n\n1\nc\nH 0 0 0\n").unwrap().comment(), "c");
    }

    #[test]
    fn a_non_numeric_count_is_an_error() {
        assert_eq!(
            parse_error("abc\n"),
            ParseError::new(
                text(1, None),
                ErrorKind::AtomCount {
                    found: "abc".into()
                }
            ),
        );
    }

    #[test]
    fn a_count_beyond_u32_is_out_of_range() {
        assert_eq!(
            parse_error("5000000000\n"),
            ParseError::new(
                text(1, None),
                ErrorKind::AtomCountRange {
                    count: 5_000_000_000
                },
            ),
        );
    }

    #[test]
    fn a_missing_comment_is_an_unexpected_eof() {
        assert_eq!(
            parse_error("1\n"),
            ParseError::new(text(2, None), ErrorKind::UnexpectedEof),
        );
    }

    #[test]
    fn too_few_atoms_is_an_unexpected_eof() {
        assert_eq!(
            parse_error("2\nc\nH 0 0 0\n"),
            ParseError::new(text(4, None), ErrorKind::UnexpectedEof),
        );
    }

    #[test]
    fn an_atom_line_without_four_fields_is_an_error() {
        assert_eq!(
            parse_error("1\nc\nH 0 0\n"),
            ParseError::new(text(3, None), ErrorKind::FieldCount { found: 3 }),
        );
        assert_eq!(
            parse_error("1\nc\nH 0 0 0 0\n"),
            ParseError::new(text(3, None), ErrorKind::FieldCount { found: 5 }),
        );
    }

    #[test]
    fn a_blank_line_among_the_atoms_is_an_error() {
        assert_eq!(
            parse_error("1\nc\n\n"),
            ParseError::new(text(3, None), ErrorKind::FieldCount { found: 0 }),
        );
    }

    #[test]
    fn an_unknown_element_symbol_is_an_error() {
        assert_eq!(
            parse_error("1\nc\nXx 0 0 0\n"),
            ParseError::new(
                text(3, Some(1)),
                ErrorKind::ElementSymbol { found: "Xx".into() },
            ),
        );
    }

    #[test]
    fn an_invalid_coordinate_is_an_error() {
        assert_eq!(
            parse_error("1\nc\nH x 0 0\n"),
            ParseError::new(
                text(3, Some(3)),
                ErrorKind::Coordinate { found: "x".into() },
            ),
        );
    }

    #[test]
    fn a_non_finite_coordinate_is_rejected() {
        assert_eq!(
            parse_error("1\nc\nH inf 0 0\n"),
            ParseError::new(
                text(3, Some(3)),
                ErrorKind::Coordinate {
                    found: "inf".into()
                },
            ),
        );
    }

    #[test]
    fn systems_yields_every_frame_in_order() {
        let frames = parse_all("1\nfirst\nH 0 0 0\n1\nsecond\nO 1 1 1\n");
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[0].as_ref().unwrap().comment(), "first");
        assert_eq!(frames[1].as_ref().unwrap().comment(), "second");
    }

    #[test]
    fn system_reads_the_first_frame_and_ignores_the_rest() {
        let system = parse("1\nfirst\nH 0 0 0\n1\nsecond\nO 1 1 1\n").unwrap();
        assert_eq!(system.comment(), "first");
    }

    #[test]
    fn systems_on_empty_input_yields_no_frames() {
        assert!(parse_all("").is_empty());
    }

    #[test]
    fn systems_stops_at_the_first_malformed_frame() {
        let frames = parse_all("1\nok\nH 0 0 0\nbad\ncomment\n");
        assert_eq!(frames.len(), 2);
        assert!(frames[0].is_ok());
        assert!(frames[1].is_err());
    }

    #[test]
    fn blank_lines_between_frames_are_skipped() {
        let frames = parse_all("1\na\nH 0 0 0\n\n1\nb\nO 1 1 1\n");
        assert_eq!(frames.len(), 2);
        assert!(frames.iter().all(Result::is_ok));
    }

    #[test]
    fn trailing_blank_lines_yield_no_extra_frame() {
        let frames = parse_all("1\na\nH 0 0 0\n\n\n");
        assert_eq!(frames.len(), 1);
    }

    #[test]
    fn systems_is_fused_after_completion() {
        let mut systems = read("1\na\nH 0 0 0\n".as_bytes()).systems::<f64>();
        assert!(systems.next().is_some());
        assert!(systems.next().is_none());
        assert!(systems.next().is_none());
    }

    #[test]
    fn parsing_is_deterministic() {
        assert_eq!(parse(FRAME).unwrap(), parse(FRAME).unwrap());
    }
}

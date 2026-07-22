use std::fmt;

use vita_core::{Element, Isotope};

use crate::algorithm::composition::{Composition, Constituent};

/// What made a formula token invalid.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ErrorKind {
    /// A letter run is neither an element symbol nor the `D`/`T` alias.
    UnknownSymbol,
    /// A count or part coefficient is zero.
    ZeroQuantity,
    /// A count, scaled or summed, exceeds `u32`, or a net charge leaves
    /// `i32`.
    Overflow,
    /// A bracketed mass number is missing, exceeds `u16`, or falls below its
    /// element's atomic number.
    InvalidMassNumber,
    /// A `[` has no matching `]`.
    UnclosedBracket,
    /// A character fits no formula token, or trails its part's charge.
    UnexpectedCharacter,
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            ErrorKind::UnknownSymbol => "unknown element symbol",
            ErrorKind::ZeroQuantity => "zero quantity",
            ErrorKind::Overflow => "quantity out of range",
            ErrorKind::InvalidMassNumber => "invalid mass number",
            ErrorKind::UnclosedBracket => "unclosed bracket",
            ErrorKind::UnexpectedCharacter => "unexpected character",
        };
        f.write_str(message)
    }
}

/// A formula parse failure.
pub type ParseError = crate::algorithm::notation::ParseError<ErrorKind>;

/// Reads a molecular formula into the [`Composition`] it denotes.
///
/// Reads the notation's whole established vocabulary, of which
/// [`write`](super::write) emits one canonical dialect: units in any order,
/// repeated symbols accumulating (`CH3COOH` counts C2H4O2), the deuterium
/// and tritium aliases `D` and `T`, dot-separated parts whose counts and
/// charges sum (`Na+.Cl-`), and a part-leading coefficient scaling its
/// part's counts and charge (`CuSO4.5H2O`, `2Na+`). Nothing outside the
/// vocabulary is repaired: case and whitespace bear meaning, so `h2o` and
/// blanks are rejected. Symbols match greedily — `Co` is cobalt, never
/// carbon and oxygen — and a real symbol outranks an alias: `Db` is
/// dubnium. A charge closes its part.
///
/// # Errors
///
/// Returns the [`ParseError`] locating and classifying the first offending
/// token.
///
/// # Complexity
///
/// O(L + U · log U) time and O(U) space, over the input's `L` bytes and `U`
/// parsed units; the log factor merges repeated constituents.
pub fn parse(text: &str) -> Result<Composition, ParseError> {
    let bytes = text.as_bytes();
    let mut units: Vec<(Constituent, u32, usize)> = Vec::new();
    let mut charge = 0i32;
    let mut coefficient = 1u32;
    let mut fresh = true;
    let mut charged = false;
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'.' => {
                coefficient = 1;
                fresh = true;
                charged = false;
                i += 1;
            }
            _ if charged => return Err(ParseError::new(i, ErrorKind::UnexpectedCharacter)),
            b'0'..=b'9' if fresh => {
                let (end, value) = number(bytes, i);
                coefficient = match value {
                    Some(0) => return Err(ParseError::new(i, ErrorKind::ZeroQuantity)),
                    Some(coefficient) => coefficient,
                    None => return Err(ParseError::new(i, ErrorKind::Overflow)),
                };
                fresh = false;
                i = end;
            }
            b'+' | b'-' => {
                let sign = if bytes[i] == b'+' { 1i64 } else { -1i64 };
                let (end, value) = number(bytes, i + 1);
                let magnitude = if end == i + 1 {
                    1
                } else {
                    value.ok_or(ParseError::new(i, ErrorKind::Overflow))?
                };
                let total = i64::from(magnitude)
                    .checked_mul(i64::from(coefficient))
                    .and_then(|delta| i64::from(charge).checked_add(sign * delta))
                    .ok_or(ParseError::new(i, ErrorKind::Overflow))?;
                charge =
                    i32::try_from(total).map_err(|_| ParseError::new(i, ErrorKind::Overflow))?;
                charged = true;
                fresh = false;
                i = end;
            }
            b'[' => {
                let (constituent, end) = nuclide(text, bytes, i)?;
                let (count, end) = count(bytes, end)?;
                units.push((constituent, scale(count, coefficient, i)?, i));
                fresh = false;
                i = end;
            }
            b'A'..=b'Z' => {
                let (constituent, end) = symbol(text, bytes, i)?;
                let (count, end) = count(bytes, end)?;
                units.push((constituent, scale(count, coefficient, i)?, i));
                fresh = false;
                i = end;
            }
            _ => return Err(ParseError::new(i, ErrorKind::UnexpectedCharacter)),
        }
    }

    units.sort_unstable_by_key(|&(constituent, _, at)| (constituent, at));
    let mut merged: Vec<(Constituent, u32)> = Vec::new();
    for (constituent, unit_count, at) in units {
        match merged.last_mut() {
            Some((last, total)) if *last == constituent => {
                *total = total
                    .checked_add(unit_count)
                    .ok_or(ParseError::new(at, ErrorKind::Overflow))?;
            }
            _ => merged.push((constituent, unit_count)),
        }
    }
    Ok(Composition::from_counts(merged, charge))
}

fn symbol(text: &str, bytes: &[u8], start: usize) -> Result<(Constituent, usize), ParseError> {
    let mut run = start + 1;
    while run < bytes.len() && bytes[run].is_ascii_lowercase() && run - start < 3 {
        run += 1;
    }
    for length in (1..=run - start).rev() {
        if let Some(element) = Element::from_symbol(&text[start..start + length]) {
            return Ok((Constituent::Element(element), start + length));
        }
        if length == 1 {
            let mass_number = match bytes[start] {
                b'D' => 2,
                b'T' => 3,
                _ => return Err(ParseError::new(start, ErrorKind::UnknownSymbol)),
            };
            return Ok((heavy_hydrogen(mass_number), start + 1));
        }
    }
    Err(ParseError::new(start, ErrorKind::UnknownSymbol))
}

fn nuclide(text: &str, bytes: &[u8], bracket: usize) -> Result<(Constituent, usize), ParseError> {
    let digits = bracket + 1;
    let (digits_end, value) = number(bytes, digits);
    if digits_end == digits {
        return Err(ParseError::new(digits, ErrorKind::InvalidMassNumber));
    }
    let mut end = digits_end;
    while end < bytes.len() && bytes[end] != b']' {
        end += 1;
    }
    if end == bytes.len() {
        return Err(ParseError::new(bracket, ErrorKind::UnclosedBracket));
    }
    let element = Element::from_symbol(&text[digits_end..end])
        .ok_or(ParseError::new(digits_end, ErrorKind::UnknownSymbol))?;
    let mass_number = value
        .and_then(|value| u16::try_from(value).ok())
        .ok_or(ParseError::new(digits, ErrorKind::InvalidMassNumber))?;
    let isotope = Isotope::new(element, mass_number)
        .ok_or(ParseError::new(digits, ErrorKind::InvalidMassNumber))?;
    Ok((Constituent::Nuclide(isotope), end + 1))
}

fn count(bytes: &[u8], start: usize) -> Result<(u32, usize), ParseError> {
    let (end, value) = number(bytes, start);
    if end == start {
        return Ok((1, start));
    }
    match value {
        Some(0) => Err(ParseError::new(start, ErrorKind::ZeroQuantity)),
        Some(count) => Ok((count, end)),
        None => Err(ParseError::new(start, ErrorKind::Overflow)),
    }
}

fn scale(count: u32, coefficient: u32, at: usize) -> Result<u32, ParseError> {
    count
        .checked_mul(coefficient)
        .ok_or(ParseError::new(at, ErrorKind::Overflow))
}

fn heavy_hydrogen(mass_number: u16) -> Constituent {
    let hydrogen = Element::new(1).expect("atomic number one is hydrogen");
    let isotope = Isotope::new(hydrogen, mass_number).expect("the alias mass numbers exceed one");
    Constituent::Nuclide(isotope)
}

fn number(bytes: &[u8], start: usize) -> (usize, Option<u32>) {
    let mut end = start;
    let mut value = Some(0u32);
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        let digit = u32::from(bytes[end] - b'0');
        value = value
            .and_then(|value| value.checked_mul(10))
            .and_then(|value| value.checked_add(digit));
        end += 1;
    }
    (end, value)
}

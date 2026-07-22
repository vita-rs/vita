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

#[cfg(test)]
mod tests {
    use super::*;

    use crate::algorithm::notation::formula::write;

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    fn natural(symbol: &str) -> Constituent {
        Constituent::Element(elem(symbol))
    }

    fn nuclide(symbol: &str, mass_number: u16) -> Constituent {
        Constituent::Nuclide(Isotope::new(elem(symbol), mass_number).unwrap())
    }

    fn water() -> Composition {
        Composition::from_counts([(natural("H"), 2), (natural("O"), 1)], 0)
    }

    fn sulfate() -> Composition {
        Composition::from_counts([(natural("S"), 1), (natural("O"), 4)], -2)
    }

    fn heavy_benzene() -> Composition {
        Composition::from_counts(
            [(natural("C"), 6), (natural("H"), 5), (nuclide("H", 2), 1)],
            0,
        )
    }

    #[test]
    fn an_empty_string_parses_to_the_empty_composition() {
        assert_eq!(parse(""), Ok(Composition::from_counts([], 0)));
    }

    #[test]
    fn a_lone_symbol_parses_to_one_atom() {
        assert_eq!(
            parse("O"),
            Ok(Composition::from_counts([(natural("O"), 1)], 0))
        );
    }

    #[test]
    fn a_count_follows_its_symbol() {
        assert_eq!(
            parse("O2"),
            Ok(Composition::from_counts([(natural("O"), 2)], 0))
        );
    }

    #[test]
    fn a_two_letter_symbol_parses_whole() {
        assert_eq!(
            parse("Cl"),
            Ok(Composition::from_counts([(natural("Cl"), 1)], 0))
        );
    }

    #[test]
    fn a_bracketed_nuclide_parses() {
        assert_eq!(
            parse("[13C]"),
            Ok(Composition::from_counts([(nuclide("C", 13), 1)], 0))
        );
    }

    #[test]
    fn a_count_multiplies_a_nuclide() {
        assert_eq!(
            parse("[2H]2"),
            Ok(Composition::from_counts([(nuclide("H", 2), 2)], 0))
        );
    }

    #[test]
    fn deuterium_aliases_to_hydrogen_two() {
        assert_eq!(
            parse("D"),
            Ok(Composition::from_counts([(nuclide("H", 2), 1)], 0))
        );
    }

    #[test]
    fn tritium_aliases_to_hydrogen_three() {
        assert_eq!(
            parse("T"),
            Ok(Composition::from_counts([(nuclide("H", 3), 1)], 0))
        );
    }

    #[test]
    fn an_alias_takes_a_count() {
        let heavy_water = Composition::from_counts([(nuclide("H", 2), 2), (natural("O"), 1)], 0);
        assert_eq!(parse("D2O"), Ok(heavy_water));
    }

    #[test]
    fn a_sign_alone_parses_to_a_unit_charge() {
        assert_eq!(parse("+"), Ok(Composition::from_counts([], 1)));
        assert_eq!(parse("-"), Ok(Composition::from_counts([], -1)));
    }

    #[test]
    fn a_charge_magnitude_follows_its_sign() {
        assert_eq!(
            parse("Fe+3"),
            Ok(Composition::from_counts([(natural("Fe"), 1)], 3))
        );
    }

    #[test]
    fn a_dot_separates_summed_parts() {
        let expected = Composition::from_counts(
            [
                (natural("Na"), 1),
                (natural("Cl"), 1),
                (natural("H"), 2),
                (natural("O"), 1),
            ],
            0,
        );
        assert_eq!(parse("NaCl.H2O"), Ok(expected));
    }

    #[test]
    fn charges_sum_across_parts() {
        let expected = Composition::from_counts([(natural("Na"), 1), (natural("Cl"), 1)], 0);
        assert_eq!(parse("Na+.Cl-"), Ok(expected));
    }

    #[test]
    fn repeated_symbols_accumulate() {
        let acetic_acid =
            Composition::from_counts([(natural("C"), 2), (natural("H"), 4), (natural("O"), 2)], 0);
        assert_eq!(parse("CH3COOH"), Ok(acetic_acid));
    }

    #[test]
    fn a_coefficient_multiplies_its_part() {
        let hydrate = Composition::from_counts(
            [
                (natural("Cu"), 1),
                (natural("S"), 1),
                (natural("O"), 9),
                (natural("H"), 10),
            ],
            0,
        );
        assert_eq!(parse("CuSO4.5H2O"), Ok(hydrate));
    }

    #[test]
    fn a_coefficient_scales_its_parts_charge() {
        assert_eq!(
            parse("2Na+"),
            Ok(Composition::from_counts([(natural("Na"), 2)], 2))
        );
    }

    #[test]
    fn an_error_kind_displays_its_condition() {
        assert_eq!(
            ErrorKind::UnknownSymbol.to_string(),
            "unknown element symbol"
        );
    }

    #[test]
    fn an_unknown_symbol_is_rejected() {
        assert_eq!(
            parse("Q"),
            Err(ParseError::new(0, ErrorKind::UnknownSymbol))
        );
    }

    #[test]
    fn an_unknown_bracketed_symbol_is_rejected() {
        assert_eq!(
            parse("[13Qq]"),
            Err(ParseError::new(3, ErrorKind::UnknownSymbol))
        );
    }

    #[test]
    fn a_zero_count_is_rejected() {
        assert_eq!(
            parse("C0"),
            Err(ParseError::new(1, ErrorKind::ZeroQuantity))
        );
    }

    #[test]
    fn a_zero_coefficient_is_rejected() {
        assert_eq!(
            parse("0H2O"),
            Err(ParseError::new(0, ErrorKind::ZeroQuantity))
        );
    }

    #[test]
    fn an_overflowing_count_is_rejected() {
        assert_eq!(
            parse("C4294967296"),
            Err(ParseError::new(1, ErrorKind::Overflow))
        );
    }

    #[test]
    fn an_overflowing_scaled_count_is_rejected() {
        assert_eq!(
            parse("5C4000000000"),
            Err(ParseError::new(1, ErrorKind::Overflow))
        );
    }

    #[test]
    fn an_overflowing_merged_count_is_rejected() {
        assert_eq!(
            parse("C4294967295C"),
            Err(ParseError::new(11, ErrorKind::Overflow))
        );
    }

    #[test]
    fn a_missing_mass_number_is_rejected() {
        assert_eq!(
            parse("[C]"),
            Err(ParseError::new(1, ErrorKind::InvalidMassNumber))
        );
    }

    #[test]
    fn a_mass_number_below_the_atomic_number_is_rejected() {
        assert_eq!(
            parse("[1C]"),
            Err(ParseError::new(1, ErrorKind::InvalidMassNumber))
        );
    }

    #[test]
    fn an_overflowing_mass_number_is_rejected() {
        assert_eq!(
            parse("[100000H]"),
            Err(ParseError::new(1, ErrorKind::InvalidMassNumber))
        );
    }

    #[test]
    fn an_unclosed_bracket_is_rejected() {
        assert_eq!(
            parse("[13C"),
            Err(ParseError::new(0, ErrorKind::UnclosedBracket))
        );
    }

    #[test]
    fn an_overflowing_charge_magnitude_is_rejected() {
        assert_eq!(
            parse("+4294967296"),
            Err(ParseError::new(0, ErrorKind::Overflow))
        );
    }

    #[test]
    fn an_overflowing_net_charge_is_rejected() {
        assert_eq!(
            parse("+2000000000.+2000000000"),
            Err(ParseError::new(12, ErrorKind::Overflow))
        );
    }

    #[test]
    fn a_lowercase_letter_is_rejected() {
        assert_eq!(
            parse("h2o"),
            Err(ParseError::new(0, ErrorKind::UnexpectedCharacter))
        );
    }

    #[test]
    fn a_unit_after_a_charge_is_rejected() {
        assert_eq!(
            parse("Na+K"),
            Err(ParseError::new(3, ErrorKind::UnexpectedCharacter))
        );
    }

    #[test]
    fn whitespace_is_rejected() {
        assert_eq!(
            parse("H2 O"),
            Err(ParseError::new(2, ErrorKind::UnexpectedCharacter))
        );
    }

    #[test]
    fn the_longest_symbol_wins() {
        assert_eq!(
            parse("Co"),
            Ok(Composition::from_counts([(natural("Co"), 1)], 0))
        );
    }

    #[test]
    fn a_real_symbol_outranks_an_alias() {
        assert_eq!(
            parse("Db"),
            Ok(Composition::from_counts([(natural("Db"), 1)], 0))
        );
    }

    #[test]
    fn an_empty_part_contributes_nothing() {
        assert_eq!(parse(".H2O"), Ok(water()));
    }

    #[test]
    fn a_bare_coefficient_part_contributes_nothing() {
        assert_eq!(parse("2.H2O"), Ok(water()));
    }

    #[test]
    fn the_most_negative_charge_parses() {
        assert_eq!(
            parse("-2147483648"),
            Ok(Composition::from_counts([], i32::MIN))
        );
    }

    #[test]
    fn a_salt_with_nuclides_counts_and_charges_parses() {
        let expected = Composition::from_counts(
            [
                (nuclide("H", 2), 2),
                (natural("O"), 1),
                (natural("Na"), 1),
                (natural("Cl"), 1),
            ],
            0,
        );
        assert_eq!(parse("[2H]2O.Na+.Cl-"), Ok(expected));
    }

    #[test]
    fn parsing_is_independent_of_unit_order() {
        assert_eq!(parse("OH2"), parse("H2O"));
    }

    #[test]
    fn parsing_recovers_what_write_produced() {
        let compositions = [
            Composition::from_counts([], 0),
            Composition::from_counts([], -3),
            water(),
            sulfate(),
            heavy_benzene(),
            Composition::from_counts([(natural("C"), 4_000_000_000)], 0),
        ];
        for composition in compositions {
            assert_eq!(parse(&write(&composition)), Ok(composition));
        }
    }
}

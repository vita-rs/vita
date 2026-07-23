use crate::algorithm::composition::{Composition, Constituent};

/// Writes a composition as its Hill-order molecular formula.
///
/// A carbon-bearing composition leads with C, then H, then the remaining
/// symbols alphabetically; a carbon-free one orders every symbol
/// alphabetically. Within an element the natural mixture precedes its
/// nuclides in ascending mass number, a nuclide renders bracketed as `[13C]`,
/// a count above one follows its symbol, and a nonzero net charge ends the
/// formula sign-first: `+`, `-2`. The empty composition writes the empty
/// string.
///
/// Every composition has exactly one rendering, and [`parse`](super::parse)
/// reads it back: `parse(&write(c)) == Ok(c)`.
///
/// # Complexity
///
/// O(N · log N + L) time and O(N + L) space, over the composition's `N`
/// distinct constituents and the output's `L` bytes; the log factor orders
/// Hill.
pub fn write(composition: &Composition) -> String {
    let mut pairs: Vec<(Constituent, u32)> = composition.iter().collect();
    let carbonaceous = pairs
        .iter()
        .any(|&(constituent, _)| constituent.element().symbol() == "C");
    pairs.sort_unstable_by_key(|&(constituent, _)| {
        (hill_rank(constituent, carbonaceous), constituent)
    });

    let mut formula = String::new();
    for (constituent, count) in pairs {
        match constituent {
            Constituent::Element(element) => formula.push_str(element.symbol()),
            Constituent::Nuclide(isotope) => {
                formula.push('[');
                formula.push_str(&isotope.mass_number().to_string());
                formula.push_str(isotope.element().symbol());
                formula.push(']');
            }
        }
        if count > 1 {
            formula.push_str(&count.to_string());
        }
    }

    let charge = composition.charge();
    if charge != 0 {
        formula.push(if charge > 0 { '+' } else { '-' });
        let magnitude = charge.unsigned_abs();
        if magnitude > 1 {
            formula.push_str(&magnitude.to_string());
        }
    }
    formula
}

fn hill_rank(constituent: Constituent, carbonaceous: bool) -> (u8, &'static str) {
    let symbol = constituent.element().symbol();
    match (carbonaceous, symbol) {
        (true, "C") => (0, symbol),
        (true, "H") => (1, symbol),
        _ => (2, symbol),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, Isotope};

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    fn natural(symbol: &str) -> Constituent {
        Constituent::Element(elem(symbol))
    }

    fn nuclide(symbol: &str, mass_number: u16) -> Constituent {
        Constituent::Nuclide(Isotope::new(elem(symbol), mass_number).unwrap())
    }

    #[test]
    fn an_empty_composition_writes_the_empty_string() {
        assert_eq!(write(&Composition::from_counts([], 0)), "");
    }

    #[test]
    fn a_bare_charge_writes_only_the_charge() {
        assert_eq!(write(&Composition::from_counts([], -1)), "-");
    }

    #[test]
    fn a_lone_element_writes_its_bare_symbol() {
        assert_eq!(
            write(&Composition::from_counts([(natural("O"), 1)], 0)),
            "O"
        );
    }

    #[test]
    fn a_count_above_one_follows_its_symbol() {
        assert_eq!(
            write(&Composition::from_counts([(natural("O"), 2)], 0)),
            "O2"
        );
    }

    #[test]
    fn a_carbonaceous_formula_leads_with_carbon_then_hydrogen() {
        let acetic_acid =
            Composition::from_counts([(natural("O"), 2), (natural("H"), 4), (natural("C"), 2)], 0);
        assert_eq!(write(&acetic_acid), "C2H4O2");
    }

    #[test]
    fn the_remaining_symbols_run_alphabetically() {
        let composition = Composition::from_counts(
            [
                (natural("N"), 1),
                (natural("Cl"), 1),
                (natural("C"), 1),
                (natural("O"), 1),
                (natural("H"), 1),
                (natural("Br"), 1),
            ],
            0,
        );
        assert_eq!(write(&composition), "CHBrClNO");
    }

    #[test]
    fn a_carbonless_formula_orders_every_symbol_alphabetically() {
        let sulfuric_acid =
            Composition::from_counts([(natural("S"), 1), (natural("O"), 4), (natural("H"), 2)], 0);
        assert_eq!(write(&sulfuric_acid), "H2O4S");
    }

    #[test]
    fn a_nuclide_writes_its_bracketed_mass_number() {
        assert_eq!(
            write(&Composition::from_counts([(nuclide("C", 13), 1)], 0)),
            "[13C]"
        );
    }

    #[test]
    fn the_natural_mixture_precedes_its_nuclides() {
        let heavy_benzene = Composition::from_counts(
            [(natural("C"), 6), (natural("H"), 5), (nuclide("H", 2), 1)],
            0,
        );
        assert_eq!(write(&heavy_benzene), "C6H5[2H]");
    }

    #[test]
    fn nuclides_ascend_by_mass_number() {
        let composition =
            Composition::from_counts([(nuclide("C", 13), 1), (nuclide("C", 12), 1)], 0);
        assert_eq!(write(&composition), "[12C][13C]");
    }

    #[test]
    fn a_unit_charge_writes_only_its_sign() {
        let ammonium = Composition::from_counts([(natural("N"), 1), (natural("H"), 4)], 1);
        assert_eq!(write(&ammonium), "H4N+");
    }

    #[test]
    fn a_larger_charge_writes_its_magnitude() {
        let sulfate = Composition::from_counts([(natural("S"), 1), (natural("O"), 4)], -2);
        assert_eq!(write(&sulfate), "O4S-2");
    }

    #[test]
    fn a_neutral_composition_writes_no_charge() {
        let water = Composition::from_counts([(natural("H"), 2), (natural("O"), 1)], 0);
        assert_eq!(write(&water), "H2O");
    }

    #[test]
    fn hydrogen_takes_no_precedence_without_carbon() {
        let borane = Composition::from_counts([(natural("H"), 3), (natural("B"), 1)], 0);
        assert_eq!(write(&borane), "BH3");
    }

    #[test]
    fn the_formula_is_independent_of_input_order() {
        let forward = Composition::from_counts([(natural("H"), 2), (natural("O"), 1)], 0);
        let backward = Composition::from_counts([(natural("O"), 1), (natural("H"), 2)], 0);
        assert_eq!(write(&forward), write(&backward));
    }
}

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

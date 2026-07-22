use std::cmp::Ordering;
use std::ops::{Add, AddAssign};

use vita_core::{Element, Isotope};

use crate::algorithm::utils::SortedMap;

/// An atomic species as a composition counts it: an element at its natural
/// isotopic mixture, or one specific nuclide.
///
/// The two precisions are distinct species — deuterated benzene counts five
/// natural hydrogens and one deuterium — and order by element, the natural
/// mixture before any of its nuclides, then by ascending mass number.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Constituent {
    /// An element at its natural isotopic mixture.
    Element(Element),
    /// One specific nuclide.
    Nuclide(Isotope),
}

impl Constituent {
    /// The element, at either isotopic precision.
    pub fn element(self) -> Element {
        match self {
            Constituent::Element(element) => element,
            Constituent::Nuclide(isotope) => isotope.element(),
        }
    }

    /// The mass number of a specific nuclide.
    ///
    /// Returns `None` for the natural isotopic mixture.
    pub fn mass_number(self) -> Option<u16> {
        match self {
            Constituent::Element(_) => None,
            Constituent::Nuclide(isotope) => Some(isotope.mass_number()),
        }
    }
}

impl Ord for Constituent {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.element(), self.mass_number()).cmp(&(other.element(), other.mass_number()))
    }
}

impl PartialOrd for Constituent {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl From<Element> for Constituent {
    fn from(element: Element) -> Self {
        Constituent::Element(element)
    }
}

impl From<Isotope> for Constituent {
    fn from(isotope: Isotope) -> Self {
        Constituent::Nuclide(isotope)
    }
}

/// What a molecule is made of: how many atoms of each constituent, and the
/// net charge, with connectivity forgotten.
///
/// The coarsest projection of a molecule — atoms counted by species, bonds
/// dropped — so constitutional isomers share one composition, and the
/// molecular formula is its rendering. Counts live over [`Constituent`]s, so
/// natural-mixture atoms and isotopically labelled ones tally separately. The
/// net charge is the sum of the formal charges, on which every Lewis form of
/// a molecule agrees.
///
/// The empty composition is [`Default`]; [`AddAssign`] pools fragments —
/// assembling pieces or balancing totals — and
/// [`empirical`](Self::empirical) reduces to the smallest whole-number ratio.
///
/// Obtain via [`elemental`](super::elemental), [`isotopic`](super::isotopic),
/// or, from raw counts, [`from_counts`](Self::from_counts).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Composition {
    counts: SortedMap<Constituent, u32>,
    charge: i32,
}

impl Composition {
    /// Builds a composition from `(constituent, count)` pairs and a net
    /// charge, summing repeated constituents and dropping zero counts.
    ///
    /// Reads back what [`iter`](Self::iter) and [`charge`](Self::charge)
    /// write out, and assembles compositions no molecule produced.
    ///
    /// # Complexity
    ///
    /// O(K · log K) time and O(K) space, over the `K` pairs.
    pub fn from_counts(counts: impl IntoIterator<Item = (Constituent, u32)>, charge: i32) -> Self {
        let mut pairs: Vec<(Constituent, u32)> =
            counts.into_iter().filter(|&(_, count)| count > 0).collect();
        pairs.sort_unstable_by_key(|&(constituent, _)| constituent);
        Composition {
            counts: SortedMap::from_pairs(
                pairs
                    .chunk_by(|a, b| a.0 == b.0)
                    .map(|run| (run[0].0, run.iter().map(|&(_, count)| count).sum())),
            ),
            charge,
        }
    }

    /// The number of distinct constituents.
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns `true` if the composition counts no atoms.
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// The number of atoms counted as `constituent`, or `0` if it is absent.
    ///
    /// Exact over the isotopic precision: the natural mixture and its
    /// nuclides are distinct species. For an element's total across both, use
    /// [`element_count`](Self::element_count).
    pub fn count(&self, constituent: Constituent) -> u32 {
        self.counts.get(&constituent).copied().unwrap_or(0)
    }

    /// The number of atoms of `element` at any isotopic precision, or `0` if
    /// it is absent.
    pub fn element_count(&self, element: Element) -> u32 {
        self.counts
            .iter()
            .filter(|(constituent, _)| constituent.element() == element)
            .map(|(_, &count)| count)
            .sum()
    }

    /// The total number of atoms, multiplicities included.
    pub fn atom_count(&self) -> u32 {
        self.counts.iter().map(|(_, &count)| count).sum()
    }

    /// The net charge.
    pub fn charge(&self) -> i32 {
        self.charge
    }

    /// Iterates `(constituent, count)` pairs in ascending constituent order.
    pub fn iter(&self) -> impl Iterator<Item = (Constituent, u32)> + '_ {
        self.counts
            .iter()
            .map(|(&constituent, &count)| (constituent, count))
    }

    /// The empirical composition: counts and charge divided by their greatest
    /// common divisor, the smallest whole numbers in the same proportions.
    ///
    /// Acetic acid's C2H4O2 reduces to CH2O, oxalate's C2O4 with charge −2 to
    /// CO2 with charge −1; a composition with coprime counts returns
    /// unchanged.
    ///
    /// # Complexity
    ///
    /// O(N · log N + N · log C) time and O(N) space, over the `N` distinct
    /// constituents, where `C` bounds the counts and the charge magnitude;
    /// Euclid's algorithm contributes the log C.
    pub fn empirical(&self) -> Composition {
        let divisor = self
            .counts
            .iter()
            .map(|(_, &count)| count)
            .fold(self.charge.unsigned_abs(), gcd);
        if divisor <= 1 {
            return self.clone();
        }
        Composition {
            counts: SortedMap::from_pairs(
                self.counts
                    .iter()
                    .map(|(&constituent, &count)| (constituent, count / divisor)),
            ),
            charge: (i64::from(self.charge) / i64::from(divisor)) as i32,
        }
    }
}

impl Default for Composition {
    fn default() -> Self {
        Self::from_counts([], 0)
    }
}

impl AddAssign<&Composition> for Composition {
    fn add_assign(&mut self, other: &Composition) {
        *self =
            Composition::from_counts(self.iter().chain(other.iter()), self.charge + other.charge);
    }
}

impl Add<&Composition> for &Composition {
    type Output = Composition;

    fn add(self, other: &Composition) -> Composition {
        Composition::from_counts(self.iter().chain(other.iter()), self.charge + other.charge)
    }
}

fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    fn natural(symbol: &str) -> Constituent {
        Constituent::Element(elem(symbol))
    }

    fn nuclide(symbol: &str, mass_number: u16) -> Constituent {
        Constituent::Nuclide(Isotope::new(elem(symbol), mass_number).unwrap())
    }

    fn empty() -> Composition {
        Composition::from_counts([], 0)
    }

    fn water() -> Composition {
        Composition::from_counts([(natural("H"), 2), (natural("O"), 1)], 0)
    }

    fn sulfate() -> Composition {
        Composition::from_counts([(natural("S"), 1), (natural("O"), 4)], -2)
    }

    fn acetic_acid() -> Composition {
        Composition::from_counts([(natural("C"), 2), (natural("H"), 4), (natural("O"), 2)], 0)
    }

    fn oxalate() -> Composition {
        Composition::from_counts([(natural("C"), 2), (natural("O"), 4)], -2)
    }

    fn heavy_benzene() -> Composition {
        Composition::from_counts(
            [(natural("C"), 6), (natural("H"), 5), (nuclide("H", 2), 1)],
            0,
        )
    }

    #[test]
    fn an_empty_composition_has_no_constituents() {
        assert_eq!(empty().len(), 0);
    }

    #[test]
    fn an_empty_composition_is_empty() {
        assert!(empty().is_empty());
    }

    #[test]
    fn an_empty_composition_counts_no_atoms() {
        assert_eq!(empty().atom_count(), 0);
    }

    #[test]
    fn an_empty_composition_yields_no_pairs() {
        assert!(empty().iter().next().is_none());
    }

    #[test]
    fn empirical_of_an_empty_composition_is_empty() {
        assert_eq!(empty().empirical(), empty());
    }

    #[test]
    fn the_default_composition_is_empty_and_chargeless() {
        assert_eq!(Composition::default(), empty());
    }

    #[test]
    fn constituent_element_spans_both_precisions() {
        assert_eq!(natural("C").element(), elem("C"));
        assert_eq!(nuclide("C", 13).element(), elem("C"));
    }

    #[test]
    fn constituent_mass_number_is_some_for_a_nuclide() {
        assert_eq!(nuclide("H", 2).mass_number().unwrap(), 2);
    }

    #[test]
    fn constituent_mass_number_is_none_for_the_natural_mixture() {
        assert!(natural("H").mass_number().is_none());
    }

    #[test]
    fn an_element_converts_into_a_natural_constituent() {
        assert_eq!(Constituent::from(elem("C")), natural("C"));
    }

    #[test]
    fn an_isotope_converts_into_a_nuclide_constituent() {
        let isotope = Isotope::new(elem("H"), 2).unwrap();
        assert_eq!(Constituent::from(isotope), nuclide("H", 2));
    }

    #[test]
    fn from_counts_sums_repeated_constituents() {
        let composition = Composition::from_counts([(natural("H"), 1), (natural("H"), 1)], 0);
        assert_eq!(composition.count(natural("H")), 2);
    }

    #[test]
    fn len_counts_distinct_constituents() {
        assert_eq!(heavy_benzene().len(), 3);
    }

    #[test]
    fn count_returns_the_multiplicity_of_a_present_constituent() {
        assert_eq!(water().count(natural("H")), 2);
    }

    #[test]
    fn element_count_totals_the_mixture_and_its_nuclides() {
        assert_eq!(heavy_benzene().element_count(elem("H")), 6);
    }

    #[test]
    fn atom_count_totals_all_multiplicities() {
        assert_eq!(water().atom_count(), 3);
    }

    #[test]
    fn charge_returns_the_net_charge() {
        assert_eq!(sulfate().charge(), -2);
    }

    #[test]
    fn iter_yields_each_constituent_with_its_count() {
        let pairs: Vec<(Constituent, u32)> = water().iter().collect();
        assert_eq!(pairs, vec![(natural("H"), 2), (natural("O"), 1)]);
    }

    #[test]
    fn empirical_reduces_counts_to_their_smallest_ratio() {
        let expected =
            Composition::from_counts([(natural("C"), 1), (natural("H"), 2), (natural("O"), 1)], 0);
        assert_eq!(acetic_acid().empirical(), expected);
    }

    #[test]
    fn adding_compositions_sums_counts_and_charges() {
        let expected = Composition::from_counts(
            [(natural("H"), 2), (natural("O"), 5), (natural("S"), 1)],
            -2,
        );
        assert_eq!(&water() + &sulfate(), expected);
    }

    #[test]
    fn add_assign_accumulates_in_place() {
        let mut composition = water();
        composition += &water();
        assert_eq!(composition.count(natural("H")), 4);
    }

    #[test]
    fn count_is_zero_for_an_absent_constituent() {
        assert_eq!(water().count(natural("C")), 0);
    }

    #[test]
    fn count_keeps_the_mixture_and_its_nuclides_apart() {
        assert_eq!(heavy_benzene().count(natural("H")), 5);
        assert_eq!(heavy_benzene().count(nuclide("H", 2)), 1);
    }

    #[test]
    fn element_count_is_zero_for_an_absent_element() {
        assert_eq!(water().element_count(elem("C")), 0);
    }

    #[test]
    fn from_counts_drops_zero_counts() {
        assert!(Composition::from_counts([(natural("H"), 0)], 0).is_empty());
    }

    #[test]
    fn a_composition_with_atoms_is_not_empty() {
        assert!(!water().is_empty());
    }

    #[test]
    fn empirical_leaves_coprime_counts_unchanged() {
        assert_eq!(water().empirical(), water());
    }

    #[test]
    fn empirical_divides_the_charge_with_the_counts() {
        let expected = Composition::from_counts([(natural("C"), 1), (natural("O"), 2)], -1);
        assert_eq!(oxalate().empirical(), expected);
    }

    #[test]
    fn empirical_reduces_a_bare_charge_to_a_single_unit() {
        assert_eq!(
            Composition::from_counts([], -4).empirical(),
            Composition::from_counts([], -1)
        );
    }

    #[test]
    fn empirical_reduces_the_most_negative_charge_to_a_single_unit() {
        assert_eq!(
            Composition::from_counts([], i32::MIN).empirical(),
            Composition::from_counts([], -1)
        );
    }

    #[test]
    fn a_bare_charge_composition_is_empty() {
        assert!(Composition::from_counts([], 1).is_empty());
    }

    #[test]
    fn the_natural_mixture_orders_before_its_nuclides() {
        assert!(natural("H") < nuclide("H", 1));
    }

    #[test]
    fn constituents_order_by_element_before_precision() {
        assert!(nuclide("H", 2) < natural("C"));
    }

    #[test]
    fn nuclides_order_by_ascending_mass_number() {
        assert!(nuclide("H", 2) < nuclide("H", 3));
    }

    #[test]
    fn the_empirical_of_a_doubled_fragment_recovers_the_fragment() {
        assert_eq!((&water() + &water()).empirical(), water());
    }

    #[test]
    fn the_composition_is_independent_of_input_order() {
        let forward = Composition::from_counts([(natural("H"), 2), (natural("O"), 1)], 0);
        let backward = Composition::from_counts([(natural("O"), 1), (natural("H"), 2)], 0);
        assert_eq!(forward, backward);
    }

    #[test]
    fn addition_is_commutative() {
        assert_eq!(&water() + &sulfate(), &sulfate() + &water());
    }

    #[test]
    fn addition_is_associative() {
        let left = &(&water() + &sulfate()) + &heavy_benzene();
        let right = &water() + &(&sulfate() + &heavy_benzene());
        assert_eq!(left, right);
    }

    #[test]
    fn the_default_composition_is_the_additive_identity() {
        assert_eq!(&water() + &Composition::default(), water());
    }

    #[test]
    fn iter_is_ordered_by_constituent() {
        let composition = Composition::from_counts(
            [
                (natural("O"), 1),
                (nuclide("H", 2), 1),
                (nuclide("C", 13), 1),
                (natural("H"), 1),
                (natural("C"), 1),
            ],
            0,
        );
        let keys: Vec<Constituent> = composition
            .iter()
            .map(|(constituent, _)| constituent)
            .collect();
        assert_eq!(
            keys,
            vec![
                natural("H"),
                nuclide("H", 2),
                natural("C"),
                nuclide("C", 13),
                natural("O"),
            ]
        );
    }
}

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

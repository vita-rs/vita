//! Structural fingerprints: hashed substructure multisets that measure
//! similarity.
//!
//! A [`Fingerprint`] counts hashed local substructures; the overlap between
//! two is their similarity. Three generators enumerate the substructures —
//! [`circular`] atom environments of growing radius, [`path`] linear walks,
//! [`pair`] atom pairs at their topological distance — and a caller's own
//! enumeration enters through [`Fingerprint::from_codes`]. Comparison runs
//! through [`FeatureVector`] — Tanimoto, Tversky, cosine — and
//! [`Fingerprint::fold`] projects onto a fixed-width [`FoldedFingerprint`]
//! where compactness outweighs counts.

mod circular;
mod feature;
mod folded;
mod pair;
mod path;
mod similarity;

pub use circular::circular;
pub use feature::Fingerprint;
pub use folded::FoldedFingerprint;
pub use pair::pair;
pub use path::path;
pub use similarity::FeatureVector;

use vita_core::SiteId;

use crate::algorithm::utils::{AdjacencyList, FxHashMap};
use crate::{BondId, HasBonds};

/// A fixed, seedless start state for the feature hash — the golden-ratio
/// constant, arbitrary but stable so codes reproduce across runs.
const SEED: u64 = 0x9e37_79b9_7f4a_7c15;

/// Avalanches a 64-bit value with the splitmix64 finalizer.
///
/// A seedless bijection with strong avalanche: a one-bit change of the input flips
/// about half the output bits. Feature codes fold their content through it, so they
/// reproduce across runs and machines — the property that lets a fingerprint persist
/// and compare against one built elsewhere.
#[inline]
fn mix(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Folds a sequence of words into one feature code.
///
/// A Merkle–Damgård chain over [`mix`]: each word avalanches into the running
/// state, so the code depends on the words and their order. Generators present
/// the words in a canonical order — sorted, or direction-folded — to make the
/// code independent of how the molecule happens to be indexed.
fn combine(words: impl IntoIterator<Item = u64>) -> u64 {
    words
        .into_iter()
        .fold(SEED, |state, word| mix(state ^ word))
}

/// A molecule's sites and bonds in order, with the adjacency over their indices —
/// the dense `0..n` form the generators consume.
struct Indexed {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
    adjacency: AdjacencyList,
}

/// Indexes a molecule into contiguous vertices.
///
/// Site `sites[i]` becomes vertex `i` and bond `bonds[e]` becomes edge `e`, the
/// indices [`Indexed::adjacency`] is built over.
fn index<M: HasBonds>(mol: &M) -> Indexed {
    let sites: Vec<SiteId> = mol.sites().collect();
    let position: FxHashMap<SiteId, usize> =
        sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let bonds: Vec<BondId> = mol.bonds().collect();
    let adjacency = AdjacencyList::build(
        sites.len(),
        bonds.iter().enumerate().map(|(edge, &bond)| {
            let (a, b) = mol.bond_endpoints(bond);
            (edge, position[&a], position[&b])
        }),
    );
    Indexed {
        sites,
        bonds,
        adjacency,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasBonds for Mol {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.endpoints[i]
        }
    }

    fn empty() -> Mol {
        Mol {
            sites: Vec::new(),
            bonds: Vec::new(),
            endpoints: Vec::new(),
        }
    }

    fn chain() -> Mol {
        Mol {
            sites: vec![s(5), s(2), s(9)],
            bonds: vec![b(7), b(4)],
            endpoints: vec![(s(5), s(2)), (s(2), s(9))],
        }
    }

    #[test]
    fn an_empty_molecule_indexes_to_nothing() {
        let indexed = index(&empty());
        assert!(indexed.sites.is_empty());
        assert!(indexed.bonds.is_empty());
        assert_eq!(indexed.adjacency.len(), 0);
    }

    #[test]
    fn mix_of_zero_is_zero() {
        assert_eq!(mix(0), 0);
    }

    #[test]
    fn combining_no_words_yields_the_seed() {
        assert_eq!(combine([]), SEED);
    }

    #[test]
    fn mix_is_stable_across_builds() {
        assert_eq!(mix(1), 0x5692_161d_100b_05e5);
    }

    #[test]
    fn combine_is_stable_across_builds() {
        assert_eq!(combine([1, 2]), 0x64b1_7449_5128_d35e);
    }

    #[test]
    fn combine_depends_on_the_words() {
        assert_ne!(combine([1]), combine([2]));
    }

    #[test]
    fn combine_depends_on_word_order() {
        assert_ne!(combine([1, 2]), combine([2, 1]));
    }

    #[test]
    fn combine_depends_on_word_count() {
        assert_ne!(combine([1]), combine([1, 1]));
    }

    #[test]
    fn index_preserves_site_order() {
        assert_eq!(index(&chain()).sites, vec![s(5), s(2), s(9)]);
    }

    #[test]
    fn index_preserves_bond_order() {
        assert_eq!(index(&chain()).bonds, vec![b(7), b(4)]);
    }

    #[test]
    fn index_wires_the_adjacency_over_the_indices() {
        let indexed = index(&chain());
        assert_eq!(indexed.adjacency.neighbors(0), [(0, 1)]);
        assert_eq!(indexed.adjacency.neighbors(1), [(0, 0), (1, 2)]);
        assert_eq!(indexed.adjacency.neighbors(2), [(1, 1)]);
    }
}

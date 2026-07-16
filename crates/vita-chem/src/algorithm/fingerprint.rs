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

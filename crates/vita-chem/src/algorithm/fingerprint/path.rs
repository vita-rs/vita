use vita_core::SiteId;

use super::{Fingerprint, Indexed, combine, index};
use crate::algorithm::utils::{AdjacencyList, BitSet};
use crate::{BondId, HasBonds};

/// The path fingerprint of a molecule — features from its linear substructures.
///
/// Each feature hashes a simple path of up to `length` bonds by the alternating
/// sequence of `site_seed` and `bond_seed` values along it, read in whichever
/// direction is smaller so both ends give the same code. A feature's multiplicity
/// is the number of paths realising it. Single atoms are the concern of
/// [`circular`](super::circular) at radius `0`, so paths run from one bond.
///
/// Independent of input order: codes derive from the colouring along the
/// direction-folded path, and each undirected path is emitted once, from its
/// lower-indexed end.
///
/// # Complexity
///
/// O(V · Δ^L · L + P · log P) time and O(V + L + P) space, over the molecule's
/// `V` sites, maximum degree `Δ`, path length `L`, and the `P` paths emitted.
pub fn path<M: HasBonds>(
    mol: &M,
    site_seed: impl Fn(SiteId) -> u64,
    bond_seed: impl Fn(BondId) -> u64,
    length: usize,
) -> Fingerprint {
    let Indexed {
        sites,
        bonds,
        adjacency,
    } = index(mol);
    let site_seeds: Vec<u64> = sites.iter().map(|&site| site_seed(site)).collect();
    let bond_seeds: Vec<u64> = bonds.iter().map(|&bond| bond_seed(bond)).collect();

    let mut walk = Walk {
        adjacency: &adjacency,
        site_seeds: &site_seeds,
        bond_seeds: &bond_seeds,
        length,
        visited: BitSet::zeros(sites.len()),
        trail: Vec::new(),
        codes: Vec::new(),
    };
    for (start, &seed) in site_seeds.iter().enumerate() {
        walk.visited.set(start);
        walk.trail.push(seed);
        walk.extend(start, start);
        walk.trail.pop();
        walk.visited.toggle(start);
    }
    Fingerprint::from_codes(walk.codes)
}

/// A depth-first walk enumerating the simple paths from a start vertex, carrying
/// the seed sequence of the current path and the codes emitted so far.
struct Walk<'a> {
    adjacency: &'a AdjacencyList,
    site_seeds: &'a [u64],
    bond_seeds: &'a [u64],
    length: usize,
    visited: BitSet,
    trail: Vec<u64>,
    codes: Vec<u64>,
}

impl Walk<'_> {
    /// Extends the current path from `current` (started at `start`), emitting each
    /// path of at least one bond once, from its lower-indexed end.
    fn extend(&mut self, start: usize, current: usize) {
        let bonds = self.trail.len() / 2;
        if bonds >= 1 && start < current {
            self.codes.push(code_of(&self.trail));
        }
        if bonds == self.length {
            return;
        }
        for &(edge, next) in self.adjacency.neighbors(current) {
            if !self.visited.test(next) {
                self.visited.set(next);
                self.trail.push(self.bond_seeds[edge]);
                self.trail.push(self.site_seeds[next]);
                self.extend(start, next);
                self.trail.pop();
                self.trail.pop();
                self.visited.toggle(next);
            }
        }
    }
}

/// The code of a path from its alternating seed sequence, taken in whichever
/// direction reads smaller so that either end gives the same feature.
fn code_of(trail: &[u64]) -> u64 {
    let reversed: Vec<u64> = trail.iter().rev().copied().collect();
    let canonical: &[u64] = if trail <= reversed.as_slice() {
        trail
    } else {
        reversed.as_slice()
    };
    combine(canonical.iter().copied())
}

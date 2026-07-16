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

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::fingerprint::FeatureVector;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        colors: Vec<u64>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
    }

    impl Mol {
        fn color(&self, site: SiteId) -> u64 {
            self.colors[self.sites.iter().position(|&x| x == site).unwrap()]
        }
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

    fn mol(atoms: &[(u32, u64)], bonds: &[(u32, u32, u32)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, _)| s(id)).collect(),
            colors: atoms.iter().map(|&(_, color)| color).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c)| (s(a), s(c))).collect(),
        }
    }

    fn fingerprint(m: &Mol, length: usize) -> Fingerprint {
        path(m, |site| m.color(site), |_| 0, length)
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            colors: m.colors.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().map(|&(a, c)| (c, a)).collect(),
        }
    }

    fn chain(colors: [u64; 3]) -> Mol {
        mol(
            &[(1, colors[0]), (2, colors[1]), (3, colors[2])],
            &[(1, 1, 2), (2, 2, 3)],
        )
    }

    fn triangle() -> Mol {
        mol(
            &[(1, 1), (2, 1), (3, 1)],
            &[(1, 1, 2), (2, 2, 3), (3, 1, 3)],
        )
    }

    fn star() -> Mol {
        mol(
            &[(1, 0), (2, 2), (3, 3), (4, 4), (5, 5)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
        )
    }

    #[test]
    fn an_empty_molecule_has_an_empty_fingerprint() {
        assert!(fingerprint(&mol(&[], &[]), 2).is_empty());
    }

    #[test]
    fn a_bondless_molecule_has_an_empty_fingerprint() {
        assert!(fingerprint(&mol(&[(1, 7), (2, 8)], &[]), 2).is_empty());
    }

    #[test]
    fn each_simple_path_up_to_the_length_becomes_a_feature() {
        let print = fingerprint(&chain([1, 2, 3]), 2);
        assert_eq!(print.len(), 3);
        assert_eq!(print.cardinality(), 3);
    }

    #[test]
    fn equivalent_paths_count_with_multiplicity() {
        let print = fingerprint(&chain([1, 2, 1]), 1);
        assert_eq!(print.len(), 1);
        assert_eq!(print.cardinality(), 2);
    }

    #[test]
    fn bond_colours_enter_the_path() {
        let m = mol(&[(1, 1), (2, 2)], &[(1, 1, 2)]);
        let plain = path(&m, |x| m.color(x), |_| 0, 1);
        let recoloured = path(&m, |x| m.color(x), |_| 1, 1);
        assert_ne!(plain, recoloured);
    }

    #[test]
    fn a_path_reads_the_same_from_both_ends() {
        assert_eq!(
            fingerprint(&mol(&[(1, 1), (2, 2)], &[(1, 1, 2)]), 1),
            fingerprint(&mol(&[(1, 2), (2, 1)], &[(1, 1, 2)]), 1),
        );
    }

    #[test]
    fn paths_longer_than_the_length_are_skipped() {
        assert_eq!(fingerprint(&chain([1, 2, 3]), 1).cardinality(), 2);
    }

    #[test]
    fn a_path_never_revisits_a_site() {
        assert_eq!(fingerprint(&triangle(), 3).cardinality(), 6);
    }

    #[test]
    fn a_zero_length_yields_an_empty_fingerprint() {
        assert!(fingerprint(&chain([1, 2, 3]), 0).is_empty());
    }

    #[test]
    fn a_length_beyond_the_longest_path_adds_nothing() {
        assert_eq!(fingerprint(&triangle(), 3), fingerprint(&triangle(), 2));
    }

    #[test]
    fn a_palindromic_path_is_counted_once() {
        assert_eq!(fingerprint(&chain([1, 2, 1]), 2).cardinality(), 3);
    }

    #[test]
    fn the_fingerprint_is_independent_of_input_order() {
        let m = star();
        assert_eq!(
            path(&m, |x| m.color(x), |bond| u64::from(bond.get()), 2),
            path(
                &reversed(&m),
                |x| m.color(x),
                |bond| u64::from(bond.get()),
                2
            ),
        );
    }
}

use vita_core::SiteId;

use super::{Fingerprint, combine};
use crate::HasBonds;
use crate::algorithm::topology::path::distances;

/// The atom-pair fingerprint of a molecule (Carhart).
///
/// Each feature hashes an unordered pair of atoms with the topological distance
/// between them — `(distance, lesser seed, greater seed)` — using only the atoms'
/// `site_seed` coloring; bonds enter solely through the distance. A feature's
/// multiplicity is the number of pairs realising it. Pairs in different connected
/// components, or farther apart than `max_distance`, are skipped.
///
/// Independent of input order: each unordered pair contributes once, and ordering
/// the two seeds makes the code blind to which atom is which.
///
/// # Complexity
///
/// O(V · (V + E) + V² · log V) time and O(V²) space, over the molecule's `V`
/// sites and `E` bonds — the all-pairs distances, then a scan of the `V²` pairs.
pub fn pair<M: HasBonds>(
    mol: &M,
    site_seed: impl Fn(SiteId) -> u64,
    max_distance: usize,
) -> Fingerprint {
    let matrix = distances(mol);
    let sites: Vec<SiteId> = mol.sites().collect();
    let seeds: Vec<u64> = sites.iter().map(|&site| site_seed(site)).collect();

    let mut codes: Vec<u64> = Vec::new();
    for i in 0..sites.len() {
        for j in (i + 1)..sites.len() {
            if let Some(distance) = matrix.get(sites[i], sites[j]) {
                if (1..=max_distance).contains(&distance) {
                    let lesser = seeds[i].min(seeds[j]);
                    let greater = seeds[i].max(seeds[j]);
                    codes.push(combine([distance as u64, lesser, greater]));
                }
            }
        }
    }
    Fingerprint::from_codes(codes)
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondId;
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

    fn fingerprint(m: &Mol, max_distance: usize) -> Fingerprint {
        pair(m, |site| m.color(site), max_distance)
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

    fn star() -> Mol {
        mol(
            &[(1, 0), (2, 2), (3, 3), (4, 4), (5, 5)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
        )
    }

    #[test]
    fn an_empty_molecule_has_an_empty_fingerprint() {
        assert!(fingerprint(&mol(&[], &[]), 3).is_empty());
    }

    #[test]
    fn a_single_site_forms_no_pair() {
        assert!(fingerprint(&mol(&[(1, 7)], &[]), 3).is_empty());
    }

    #[test]
    fn each_pair_within_reach_becomes_a_feature() {
        let print = fingerprint(&chain([1, 2, 3]), 2);
        assert_eq!(print.len(), 3);
        assert_eq!(print.cardinality(), 3);
    }

    #[test]
    fn equivalent_pairs_count_with_multiplicity() {
        let print = fingerprint(&chain([1, 2, 1]), 2);
        assert_eq!(print.len(), 2);
        assert_eq!(print.cardinality(), 3);
    }

    #[test]
    fn the_distance_parts_pairs_of_equal_colors() {
        let print = fingerprint(&chain([1, 1, 1]), 2);
        assert_eq!(print.len(), 2);
        assert_eq!(print.cardinality(), 3);
    }

    #[test]
    fn a_pair_is_blind_to_which_atom_bears_which_color() {
        assert_eq!(
            fingerprint(&mol(&[(1, 1), (2, 2)], &[(1, 1, 2)]), 1),
            fingerprint(&mol(&[(1, 2), (2, 1)], &[(1, 1, 2)]), 1),
        );
    }

    #[test]
    fn pairs_in_different_components_are_skipped() {
        let print = fingerprint(&mol(&[(1, 1), (2, 2), (3, 3)], &[(1, 1, 2)]), 5);
        assert_eq!(print.cardinality(), 1);
    }

    #[test]
    fn pairs_beyond_the_maximum_distance_are_skipped() {
        assert_eq!(fingerprint(&chain([1, 2, 3]), 1).cardinality(), 2);
    }

    #[test]
    fn a_zero_maximum_yields_an_empty_fingerprint() {
        assert!(fingerprint(&chain([1, 2, 3]), 0).is_empty());
    }

    #[test]
    fn the_maximum_distance_is_inclusive() {
        assert_eq!(fingerprint(&chain([1, 2, 3]), 2).cardinality(), 3);
    }

    #[test]
    fn the_fingerprint_is_independent_of_input_order() {
        let m = star();
        assert_eq!(fingerprint(&m, 2), fingerprint(&reversed(&m), 2));
    }
}

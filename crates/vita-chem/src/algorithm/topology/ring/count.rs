use vita_core::SiteId;

use crate::HasBonds;
use crate::algorithm::utils::{DisjointSet, FxHashMap};

/// The number of independent rings in a molecule — its cycle rank.
///
/// The cycle rank μ = E − V + C, for `E` bonds, `V` sites, and `C` connected
/// components, is the size of any minimum cycle basis
/// ([`Rings::len`](super::Rings::len)), obtained here without enumerating the
/// rings. Equivalently, it counts the bonds that close a cycle — those left
/// over once a spanning forest is drawn. An acyclic molecule, the empty one
/// included, has rank zero.
///
/// # Complexity
///
/// O(V + E) time and O(V) space, over the molecule's `V` sites and `E` bonds.
pub fn count<M: HasBonds>(mol: &M) -> usize {
    let index: FxHashMap<SiteId, usize> =
        mol.sites().enumerate().map(|(i, site)| (site, i)).collect();
    let mut forest = DisjointSet::new(index.len());
    let mut rank = 0;
    for bond in mol.bonds() {
        let (a, b) = mol.bond_endpoints(bond);
        if !forest.union(index[&a], index[&b]) {
            rank += 1;
        }
    }
    rank
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondId;
    use crate::topology::connectivity::components;

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
            sites: vec![],
            bonds: vec![],
            endpoints: vec![],
        }
    }

    fn single() -> Mol {
        Mol {
            sites: vec![s(1)],
            bonds: vec![],
            endpoints: vec![],
        }
    }

    fn chain() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
        }
    }

    fn forest() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(4), s(5)), (s(5), s(6))],
        }
    }

    fn triangle() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3))],
        }
    }

    fn square() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4)), (s(1), s(4))],
        }
    }

    fn tadpole() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3)), (s(1), s(4))],
        }
    }

    fn fused() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(1), s(4)),
                (s(3), s(5)),
                (s(5), s(6)),
                (s(4), s(6)),
            ],
        }
    }

    fn two_triangles() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(4), s(6)),
            ],
        }
    }

    fn cyclic_and_acyclic() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(4), s(5)),
                (s(5), s(6)),
            ],
        }
    }

    fn k4() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6)],
            endpoints: vec![
                (s(1), s(2)),
                (s(1), s(3)),
                (s(1), s(4)),
                (s(2), s(3)),
                (s(2), s(4)),
                (s(3), s(4)),
            ],
        }
    }

    fn cube() -> Mol {
        Mol {
            sites: (1..=8).map(s).collect(),
            bonds: (1..=12).map(b).collect(),
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(1), s(4)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(8)),
                (s(5), s(8)),
                (s(1), s(5)),
                (s(2), s(6)),
                (s(3), s(7)),
                (s(4), s(8)),
            ],
        }
    }

    #[test]
    fn empty_molecule_has_rank_zero() {
        assert_eq!(count(&empty()), 0);
    }

    #[test]
    fn single_site_has_rank_zero() {
        assert_eq!(count(&single()), 0);
    }

    #[test]
    fn a_single_cycle_has_rank_one() {
        assert_eq!(count(&triangle()), 1);
        assert_eq!(count(&square()), 1);
    }

    #[test]
    fn a_tree_has_rank_zero() {
        assert_eq!(count(&chain()), 0);
    }

    #[test]
    fn a_forest_has_rank_zero() {
        assert_eq!(count(&forest()), 0);
    }

    #[test]
    fn a_cycle_with_an_acyclic_tail_has_rank_one() {
        assert_eq!(count(&tadpole()), 1);
    }

    #[test]
    fn fused_rings_have_rank_two() {
        assert_eq!(count(&fused()), 2);
    }

    #[test]
    fn two_disjoint_cycles_have_rank_two() {
        assert_eq!(count(&two_triangles()), 2);
    }

    #[test]
    fn a_cyclic_and_an_acyclic_component_sum_their_ranks() {
        assert_eq!(count(&cyclic_and_acyclic()), 1);
    }

    #[test]
    fn counts_all_independent_cycles_of_a_polycyclic_graph() {
        assert_eq!(count(&k4()), 3);
        assert_eq!(count(&cube()), 5);
    }

    #[test]
    fn rank_equals_edges_minus_sites_plus_components() {
        let mol = cyclic_and_acyclic();
        let e = mol.bonds().count();
        let v = mol.sites().count();
        let c = components(&mol).len();
        assert_eq!(count(&mol), e + c - v);
    }

    #[test]
    fn count_is_independent_of_input_order() {
        let shuffled = Mol {
            sites: vec![s(6), s(5), s(4), s(3), s(2), s(1)],
            bonds: vec![b(7), b(6), b(5), b(4), b(3), b(2), b(1)],
            endpoints: vec![
                (s(4), s(6)),
                (s(5), s(6)),
                (s(3), s(5)),
                (s(1), s(4)),
                (s(3), s(4)),
                (s(2), s(3)),
                (s(1), s(2)),
            ],
        };
        assert_eq!(count(&fused()), count(&shuffled));
    }
}

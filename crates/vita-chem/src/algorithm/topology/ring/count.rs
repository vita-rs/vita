use vita_core::HasSites;

use crate::HasBonds;
use crate::topology::connectivity::components;

/// Number of independent rings in a molecule.
///
/// Equals the cycle rank of the molecular graph, μ = E − V + C, where C is the
/// number of connected components. This is the size of any minimum cycle basis
/// — the value [`Rings::len`](super::Rings::len) reports — but is obtained
/// without enumerating the rings. Acyclic and empty molecules have zero.
///
/// # Complexity
///
/// O(V + E) time.
pub fn count<M: HasBonds + HasSites>(mol: &M) -> usize {
    let v = mol.sites().count();
    let e = mol.bonds().count();
    let c = components(mol).len();
    e + c - v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BondId;
    use vita_core::SiteId;

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

    fn triangle() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3))],
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

    #[test]
    fn empty_molecule_has_no_rings() {
        assert_eq!(count(&empty()), 0);
    }

    #[test]
    fn single_site_has_no_rings() {
        assert_eq!(count(&single()), 0);
    }

    #[test]
    fn chain_has_no_rings() {
        assert_eq!(count(&chain()), 0);
    }

    #[test]
    fn triangle_has_one_ring() {
        assert_eq!(count(&triangle()), 1);
    }

    #[test]
    fn tadpole_has_one_ring() {
        assert_eq!(count(&tadpole()), 1);
    }

    #[test]
    fn fused_has_two_rings() {
        assert_eq!(count(&fused()), 2);
    }

    #[test]
    fn two_triangles_has_two_rings() {
        assert_eq!(count(&two_triangles()), 2);
    }

    #[test]
    fn count_equals_rings_len() {
        for mol in [triangle(), tadpole(), fused(), two_triangles()] {
            assert_eq!(count(&mol), super::super::rings(&mol).len());
        }
    }
}

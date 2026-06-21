use std::collections::{HashMap, HashSet};

use vita_core::{HasSites, SiteId};

use crate::{BondId, HasBonds};

/// The bridgehead atoms of a molecule.
///
/// A bridgehead is an atom where a bridge meets a ring: an endpoint of the
/// bonds two rings share, when they share two or more. Ortho-fused rings
/// (sharing one bond) and spiro rings (sharing one site) have none; a bridged
/// system such as bicyclo[2.2.2]octane has two.
///
/// # Complexity
///
/// O(V² · E) time, dominated by the minimum cycle basis it builds.
pub fn bridgeheads<M: HasBonds + HasSites>(mol: &M) -> impl Iterator<Item = SiteId> {
    let basis = super::rings(mol);
    let cycles: Vec<HashSet<BondId>> = basis
        .iter()
        .map(|r| r.bonds().iter().copied().collect())
        .collect();

    let mut heads: HashSet<SiteId> = HashSet::new();
    for (i, a) in cycles.iter().enumerate() {
        for b in &cycles[i + 1..] {
            let shared: Vec<BondId> = a.intersection(b).copied().collect();
            if shared.len() < 2 {
                continue;
            }
            let mut incident: HashMap<SiteId, usize> = HashMap::new();
            for bond in shared {
                let (u, v) = mol.bond_endpoints(bond);
                *incident.entry(u).or_default() += 1;
                *incident.entry(v).or_default() += 1;
            }
            heads.extend(
                incident
                    .into_iter()
                    .filter(|&(_, n)| n == 1)
                    .map(|(s, _)| s),
            );
        }
    }

    let mut heads: Vec<SiteId> = heads.into_iter().collect();
    heads.sort_unstable();
    heads.into_iter()
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

    fn mol(sites: &[u32], bonds: &[(u32, u32, u32)]) -> Mol {
        Mol {
            sites: sites.iter().map(|&n| s(n)).collect(),
            bonds: bonds.iter().map(|&(id, _, _)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, u, v)| (s(u), s(v))).collect(),
        }
    }

    fn chain() -> Mol {
        mol(&[1, 2, 3], &[(1, 1, 2), (2, 2, 3)])
    }

    fn square() -> Mol {
        mol(&[1, 2, 3, 4], &[(1, 1, 2), (2, 2, 3), (3, 3, 4), (4, 1, 4)])
    }

    fn fused() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 3, 4),
                (4, 1, 4),
                (5, 3, 5),
                (6, 5, 6),
                (7, 4, 6),
            ],
        )
    }

    fn spiro() -> Mol {
        mol(
            &[1, 2, 3, 4, 5],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 1, 3),
                (4, 3, 4),
                (5, 4, 5),
                (6, 3, 5),
            ],
        )
    }

    fn bicyclo222() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6, 7, 8],
            &[
                (1, 1, 3),
                (2, 3, 4),
                (3, 4, 2),
                (4, 1, 5),
                (5, 5, 6),
                (6, 6, 2),
                (7, 1, 7),
                (8, 7, 8),
                (9, 8, 2),
            ],
        )
    }

    fn bridged_square() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 3, 4),
                (4, 1, 4),
                (5, 2, 5),
                (6, 5, 6),
                (7, 4, 6),
            ],
        )
    }

    #[test]
    fn acyclic_has_no_bridgeheads() {
        assert_eq!(bridgeheads(&chain()).count(), 0);
    }

    #[test]
    fn single_ring_has_no_bridgeheads() {
        assert_eq!(bridgeheads(&square()).count(), 0);
    }

    #[test]
    fn ortho_fused_has_no_bridgeheads() {
        assert_eq!(bridgeheads(&fused()).count(), 0);
    }

    #[test]
    fn spiro_has_no_bridgeheads() {
        assert_eq!(bridgeheads(&spiro()).count(), 0);
    }

    #[test]
    fn bicyclo222_bridgeheads_are_the_bridge_ends() {
        assert_eq!(
            bridgeheads(&bicyclo222()).collect::<Vec<_>>(),
            vec![s(1), s(2)]
        );
    }

    #[test]
    fn bridged_square_bridgeheads_are_the_bridge_ends() {
        assert_eq!(
            bridgeheads(&bridged_square()).collect::<Vec<_>>(),
            vec![s(2), s(4)]
        );
    }
}

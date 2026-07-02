use std::cmp::Ordering;

use vita_core::SiteId;

use super::rings;
use crate::{BondId, HasBonds};

/// The bridgehead atoms of a molecule, in ascending order.
///
/// A bridgehead is an atom where a bridge joins a ring. When two rings of the
/// minimum cycle basis share two or more bonds, those bonds form a bridge whose
/// ends — the atoms incident to exactly one shared bond — are the bridgeheads.
/// Rings sharing a single bond (ortho-fused) or a single atom (spiro) have none,
/// so naphthalene yields nothing while bicyclo[2.2.2]octane yields two.
///
/// # Complexity
///
/// O(V² · E) time, dominated by the minimum cycle basis it builds, and
/// O(V + E) space.
pub fn bridgeheads<M: HasBonds>(mol: &M) -> impl Iterator<Item = SiteId> {
    let basis = rings(mol);
    let ring_bonds: Vec<Vec<BondId>> = basis
        .iter()
        .map(|ring| {
            let mut bonds = ring.bonds().to_vec();
            bonds.sort_unstable();
            bonds
        })
        .collect();

    let mut heads: Vec<SiteId> = Vec::new();
    for (i, a) in ring_bonds.iter().enumerate() {
        for b in &ring_bonds[i + 1..] {
            let shared = shared_bonds(a, b);
            if shared.len() < 2 {
                continue;
            }
            let mut ends: Vec<SiteId> = shared
                .iter()
                .flat_map(|&bond| {
                    let (u, v) = mol.bond_endpoints(bond);
                    [u, v]
                })
                .collect();
            ends.sort_unstable();
            heads.extend(
                ends.chunk_by(|x, y| x == y)
                    .filter(|run| run.len() == 1)
                    .map(|run| run[0]),
            );
        }
    }

    heads.sort_unstable();
    heads.dedup();
    heads.into_iter()
}

/// The bonds common to two ascending bond slices.
fn shared_bonds(a: &[BondId], b: &[BondId]) -> Vec<BondId> {
    let mut shared = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                shared.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    shared
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
            sites: vec![],
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

    fn spiro() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(3), s(5)),
            ],
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
                (s(4), s(6)),
                (s(5), s(6)),
                (s(3), s(5)),
            ],
        }
    }

    fn theta() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6)],
            endpoints: vec![
                (s(1), s(3)),
                (s(2), s(3)),
                (s(1), s(4)),
                (s(2), s(4)),
                (s(1), s(5)),
                (s(2), s(5)),
            ],
        }
    }

    fn k24() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8)],
            endpoints: vec![
                (s(1), s(3)),
                (s(2), s(3)),
                (s(1), s(4)),
                (s(2), s(4)),
                (s(1), s(5)),
                (s(2), s(5)),
                (s(1), s(6)),
                (s(2), s(6)),
            ],
        }
    }

    fn bicyclo222() -> Mol {
        Mol {
            sites: (1..=8).map(s).collect(),
            bonds: (1..=9).map(b).collect(),
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(8)),
                (s(1), s(4)),
                (s(4), s(5)),
                (s(5), s(8)),
                (s(1), s(6)),
                (s(6), s(7)),
                (s(7), s(8)),
            ],
        }
    }

    fn disconnected() -> Mol {
        Mol {
            sites: (1..=11).map(s).collect(),
            bonds: (1..=12).map(b).collect(),
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(8)),
                (s(1), s(4)),
                (s(4), s(5)),
                (s(5), s(8)),
                (s(1), s(6)),
                (s(6), s(7)),
                (s(7), s(8)),
                (s(9), s(10)),
                (s(10), s(11)),
                (s(9), s(11)),
            ],
        }
    }

    #[test]
    fn empty_molecule_has_no_bridgeheads() {
        assert_eq!(bridgeheads(&empty()).count(), 0);
    }

    #[test]
    fn single_ring_has_no_bridgeheads() {
        assert_eq!(bridgeheads(&triangle()).count(), 0);
    }

    #[test]
    fn bridgeheads_are_the_atoms_where_the_bridges_meet_the_rings() {
        assert_eq!(
            bridgeheads(&bicyclo222()).collect::<Vec<_>>(),
            vec![s(1), s(8)]
        );
    }

    #[test]
    fn an_acyclic_molecule_has_no_bridgeheads() {
        assert_eq!(bridgeheads(&chain()).count(), 0);
    }

    #[test]
    fn ortho_fused_rings_have_no_bridgeheads() {
        assert_eq!(bridgeheads(&fused()).count(), 0);
    }

    #[test]
    fn spiro_rings_have_no_bridgeheads() {
        assert_eq!(bridgeheads(&spiro()).count(), 0);
    }

    #[test]
    fn rings_sharing_exactly_two_bonds_yield_bridgeheads() {
        assert_eq!(bridgeheads(&theta()).collect::<Vec<_>>(), vec![s(1), s(2)]);
    }

    #[test]
    fn unbridged_components_add_no_bridgeheads() {
        assert_eq!(
            bridgeheads(&disconnected()).collect::<Vec<_>>(),
            vec![s(1), s(8)],
        );
    }

    #[test]
    fn a_bridgehead_shared_by_several_ring_pairs_is_listed_once() {
        assert_eq!(bridgeheads(&k24()).collect::<Vec<_>>(), vec![s(1), s(2)]);
    }

    #[test]
    fn output_is_independent_of_input_order() {
        let reordered = Mol {
            sites: (1..=8).rev().map(s).collect(),
            bonds: (1..=9).rev().map(b).collect(),
            endpoints: vec![
                (s(7), s(8)),
                (s(6), s(7)),
                (s(1), s(6)),
                (s(5), s(8)),
                (s(4), s(5)),
                (s(1), s(4)),
                (s(3), s(8)),
                (s(2), s(3)),
                (s(1), s(2)),
            ],
        };
        assert_eq!(
            bridgeheads(&bicyclo222()).collect::<Vec<_>>(),
            bridgeheads(&reordered).collect::<Vec<_>>(),
        );
    }
}

use vita_core::SiteId;

use crate::HasBonds;
use crate::algorithm::utils::{AdjacencyList, FxHashMap, SortedMap};
use crate::topology::ring::membership;

/// The structural role of an atom in a molecule's Bemis–Murcko framework.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Role {
    /// In a ring.
    Ring,
    /// On a path between ring systems, but in no ring.
    Linker,
    /// On an acyclic branch hanging off the framework.
    SideChain,
}

/// The Bemis–Murcko framework of a molecule.
///
/// Assigns every atom a structural [`Role`]: a ring atom, a linker joining ring
/// systems, or a side-chain atom on an acyclic branch. The framework itself —
/// the ring and linker atoms, given by [`sites`](Self::sites) — is what remains
/// once the side chains are stripped; an acyclic molecule has none.
///
/// Obtain via [`framework`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Framework {
    roles: SortedMap<SiteId, Role>,
}

impl Framework {
    /// Returns the structural role of `site`.
    ///
    /// Returns `None` if `site` is absent from the molecule.
    pub fn role(&self, site: SiteId) -> Option<Role> {
        self.roles.get(&site).copied()
    }

    /// The framework sites — the ring and linker atoms left once the side chains
    /// are stripped — in ascending order.
    pub fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.of_role(|role| role != Role::SideChain)
    }

    /// The linker sites — framework atoms on a path between ring systems but in
    /// no ring — in ascending order.
    pub fn linkers(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.of_role(|role| role == Role::Linker)
    }

    /// The side-chain sites — atoms on the acyclic branches stripped from the
    /// framework — in ascending order.
    pub fn side_chains(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.of_role(|role| role == Role::SideChain)
    }

    /// The sites whose role satisfies `keep`, in ascending order.
    fn of_role<'a>(
        &'a self,
        keep: impl Fn(Role) -> bool + 'a,
    ) -> impl Iterator<Item = SiteId> + 'a {
        self.roles
            .iter()
            .filter(move |&(_, &role)| keep(role))
            .map(|(&site, _)| site)
    }
}

/// Bemis–Murcko framework of a molecule.
///
/// Strips the acyclic side chains — the atoms outside the 2-core — by repeatedly
/// removing atoms of degree below two. Each surviving atom is a ring atom or a
/// linker between rings, told apart by [ring membership](membership).
///
/// # Complexity
///
/// O(V · log V + E · log E) time and O(V + E) space, over the molecule's `V`
/// sites and `E` bonds.
pub fn framework<M: HasBonds>(mol: &M) -> Framework {
    let rings = membership(mol);

    let sites: Vec<SiteId> = mol.sites().collect();
    let n = sites.len();
    let index: FxHashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let adjacency = AdjacencyList::build(
        n,
        mol.bonds().enumerate().map(|(e, bond)| {
            let (a, b) = mol.bond_endpoints(bond);
            (e, index[&a], index[&b])
        }),
    );

    let peeled = peel_side_chains(&adjacency);

    let mut is_ring = vec![false; n];
    for site in rings.sites() {
        is_ring[index[&site]] = true;
    }

    let roles = sites.iter().enumerate().map(|(i, &site)| {
        let role = if peeled[i] {
            Role::SideChain
        } else if is_ring[i] {
            Role::Ring
        } else {
            Role::Linker
        };
        (site, role)
    });

    Framework {
        roles: SortedMap::from_pairs(roles),
    }
}

/// Flags each site stripped as a side chain: the atoms outside the 2-core of the
/// molecular graph.
///
/// Repeatedly removes atoms of degree below two until only the 2-core — the ring
/// and linker atoms — remains. The flag at index `u` refers to the site whose
/// index is `u` in `adjacency`.
fn peel_side_chains(adjacency: &AdjacencyList) -> Vec<bool> {
    let n = adjacency.len();
    let mut degree: Vec<usize> = (0..n).map(|u| adjacency.neighbors(u).len()).collect();
    let mut peeled = vec![false; n];
    let mut worklist: Vec<usize> = (0..n).filter(|&u| degree[u] < 2).collect();

    while let Some(u) = worklist.pop() {
        if peeled[u] {
            continue;
        }
        peeled[u] = true;
        for &(_, v) in adjacency.neighbors(u) {
            if !peeled[v] {
                degree[v] -= 1;
                if degree[v] == 1 {
                    worklist.push(v);
                }
            }
        }
    }

    peeled
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondId;

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

    fn ring_with_tail() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3)), (s(1), s(4))],
        }
    }

    fn linked_rings() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(5), s(7)),
                (s(3), s(4)),
                (s(4), s(5)),
            ],
        }
    }

    fn decorated() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(5), s(7)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(4), s(8)),
            ],
        }
    }

    #[test]
    fn empty_molecule_has_no_framework_sites() {
        assert_eq!(framework(&empty()).sites().count(), 0);
    }

    #[test]
    fn single_atom_is_a_side_chain() {
        assert_eq!(framework(&single()).role(s(1)), Some(Role::SideChain));
    }

    #[test]
    fn ring_atom_has_the_ring_role() {
        assert_eq!(framework(&triangle()).role(s(1)), Some(Role::Ring));
    }

    #[test]
    fn linker_atom_has_the_linker_role() {
        assert_eq!(framework(&linked_rings()).role(s(4)), Some(Role::Linker));
    }

    #[test]
    fn stripped_atom_has_the_side_chain_role() {
        assert_eq!(
            framework(&ring_with_tail()).role(s(4)),
            Some(Role::SideChain)
        );
    }

    #[test]
    fn absent_site_has_no_role() {
        assert_eq!(framework(&triangle()).role(s(99)), None);
    }

    #[test]
    fn a_single_ring_has_no_linkers() {
        assert_eq!(framework(&triangle()).linkers().count(), 0);
    }

    #[test]
    fn a_ring_is_entirely_framework() {
        let fw = framework(&triangle());
        assert_eq!(fw.sites().collect::<Vec<_>>(), vec![s(1), s(2), s(3)]);
        assert_eq!(fw.side_chains().count(), 0);
    }

    #[test]
    fn an_acyclic_molecule_is_entirely_side_chains() {
        let fw = framework(&chain());
        assert_eq!(fw.sites().count(), 0);
        assert_eq!(fw.side_chains().collect::<Vec<_>>(), vec![s(1), s(2), s(3)]);
    }

    #[test]
    fn classifies_a_scaffold_with_every_role() {
        let fw = framework(&decorated());
        assert_eq!(fw.role(s(1)), Some(Role::Ring));
        assert_eq!(fw.role(s(4)), Some(Role::Linker));
        assert_eq!(fw.role(s(8)), Some(Role::SideChain));
    }

    #[test]
    fn framework_is_the_ring_and_linker_atoms() {
        let sites: Vec<SiteId> = framework(&decorated()).sites().collect();
        assert_eq!(sites, vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7)]);
    }

    #[test]
    fn linkers_lists_only_the_linker_atoms() {
        assert_eq!(
            framework(&decorated()).linkers().collect::<Vec<_>>(),
            vec![s(4)],
        );
    }

    #[test]
    fn side_chains_lists_only_the_stripped_atoms() {
        assert_eq!(
            framework(&decorated()).side_chains().collect::<Vec<_>>(),
            vec![s(8)],
        );
    }

    #[test]
    fn framework_and_side_chains_partition_the_atoms() {
        let mol = decorated();
        let fw = framework(&mol);
        for site in mol.sites() {
            assert!(fw.role(site).is_some());
        }
        let mut all: Vec<SiteId> = fw.sites().chain(fw.side_chains()).collect();
        all.sort_unstable();
        assert_eq!(all, mol.sites().collect::<Vec<_>>());
    }

    #[test]
    fn framework_is_independent_of_input_order() {
        let shuffled = Mol {
            sites: vec![s(8), s(7), s(6), s(5), s(4), s(3), s(2), s(1)],
            bonds: vec![b(9), b(8), b(7), b(6), b(5), b(4), b(3), b(2), b(1)],
            endpoints: vec![
                (s(4), s(8)),
                (s(4), s(5)),
                (s(3), s(4)),
                (s(5), s(7)),
                (s(6), s(7)),
                (s(5), s(6)),
                (s(1), s(3)),
                (s(2), s(3)),
                (s(1), s(2)),
            ],
        };
        assert_eq!(framework(&decorated()), framework(&shuffled));
    }
}

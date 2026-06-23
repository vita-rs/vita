use std::collections::{HashMap, HashSet, VecDeque};

use vita_core::{HasSites, SiteId};

use crate::HasBonds;
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
/// Classifies every atom by its structural role: in a ring, on a linker between
/// rings, or on an acyclic side chain. The framework itself — the ring and
/// linker atoms, given by [`Self::sites`] — is what remains once the side chains
/// are stripped; an acyclic molecule has none.
///
/// Obtain via [`framework`].
pub struct Framework {
    roles: HashMap<SiteId, Role>,
}

impl Framework {
    /// Returns the structural role of `site`.
    ///
    /// Returns `None` if `site` is not present in the molecule.
    pub fn role(&self, site: SiteId) -> Option<Role> {
        self.roles.get(&site).copied()
    }

    /// Iterates the framework sites: the ring and linker atoms left once the
    /// side chains are stripped.
    pub fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.roles
            .iter()
            .filter(|&(_, &r)| r != Role::SideChain)
            .map(|(&s, _)| s)
    }

    /// Iterates the linker sites: framework atoms on a path between ring systems
    /// but in no ring.
    pub fn linkers(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.roles
            .iter()
            .filter(|&(_, &r)| r == Role::Linker)
            .map(|(&s, _)| s)
    }

    /// Iterates the side-chain sites: atoms on the acyclic branches stripped from
    /// the framework.
    pub fn side_chains(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.roles
            .iter()
            .filter(|&(_, &r)| r == Role::SideChain)
            .map(|(&s, _)| s)
    }
}

/// Bemis–Murcko framework of a molecule.
///
/// Peels away the acyclic side chains by repeatedly removing sites of degree
/// below two. The surviving sites form the framework, split into the ring atoms
/// and the linker atoms that join the rings.
///
/// # Complexity
///
/// O(V + E) time.
pub fn framework<M: HasBonds + HasSites>(mol: &M) -> Framework {
    let rings = membership(mol);

    let mut adj: HashMap<SiteId, Vec<SiteId>> = HashMap::new();
    for site in mol.sites() {
        adj.entry(site).or_default();
    }
    for bond in mol.bonds() {
        let (a, b) = mol.bond_endpoints(bond);
        adj.entry(a).or_default().push(b);
        adj.entry(b).or_default().push(a);
    }

    let mut degree: HashMap<SiteId, usize> = adj.iter().map(|(&s, nbrs)| (s, nbrs.len())).collect();
    let mut peeled: HashSet<SiteId> = HashSet::new();
    let mut queue: VecDeque<SiteId> = degree
        .iter()
        .filter(|&(_, &d)| d < 2)
        .map(|(&s, _)| s)
        .collect();
    while let Some(v) = queue.pop_front() {
        if !peeled.insert(v) {
            continue;
        }
        for &u in &adj[&v] {
            if peeled.contains(&u) {
                continue;
            }
            let d = degree.get_mut(&u).unwrap();
            *d -= 1;
            if *d == 1 {
                queue.push_back(u);
            }
        }
    }

    let mut roles: HashMap<SiteId, Role> = HashMap::new();
    for site in mol.sites() {
        let role = if peeled.contains(&site) {
            Role::SideChain
        } else if rings.site(site) {
            Role::Ring
        } else {
            Role::Linker
        };
        roles.insert(site, role);
    }

    Framework { roles }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BondId;
    use std::collections::HashSet;
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

    fn empty() -> Mol {
        mol(&[], &[])
    }

    fn single() -> Mol {
        mol(&[1], &[])
    }

    fn chain() -> Mol {
        mol(&[1, 2, 3], &[(1, 1, 2), (2, 2, 3)])
    }

    fn triangle() -> Mol {
        mol(&[1, 2, 3], &[(1, 1, 2), (2, 2, 3), (3, 1, 3)])
    }

    fn ring_with_tail() -> Mol {
        mol(&[1, 2, 3, 4], &[(1, 1, 2), (2, 2, 3), (3, 1, 3), (4, 1, 4)])
    }

    fn linked_rings() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6, 7],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 1, 3),
                (4, 4, 5),
                (5, 5, 6),
                (6, 4, 6),
                (7, 3, 7),
                (8, 7, 4),
            ],
        )
    }

    #[test]
    fn empty_has_no_sites() {
        let fw = framework(&empty());
        assert!(fw.role(s(1)).is_none());
        assert_eq!(fw.sites().count(), 0);
    }

    #[test]
    fn single_site_is_a_side_chain() {
        assert_eq!(framework(&single()).role(s(1)), Some(Role::SideChain));
    }

    #[test]
    fn acyclic_molecule_is_all_side_chains() {
        let fw = framework(&chain());
        assert_eq!(fw.sites().count(), 0);
        let mut side: Vec<SiteId> = fw.side_chains().collect();
        side.sort_unstable();
        assert_eq!(side, vec![s(1), s(2), s(3)]);
    }

    #[test]
    fn ring_is_all_framework() {
        let fw = framework(&triangle());
        let mut sites: Vec<SiteId> = fw.sites().collect();
        sites.sort_unstable();
        assert_eq!(sites, vec![s(1), s(2), s(3)]);
        assert_eq!(fw.side_chains().count(), 0);
    }

    #[test]
    fn tail_atom_is_a_side_chain() {
        assert_eq!(
            framework(&ring_with_tail()).role(s(4)),
            Some(Role::SideChain)
        );
    }

    #[test]
    fn ring_atoms_are_framework() {
        let fw = framework(&ring_with_tail());
        assert_eq!(fw.role(s(1)), Some(Role::Ring));
        let mut sites: Vec<SiteId> = fw.sites().collect();
        sites.sort_unstable();
        assert_eq!(sites, vec![s(1), s(2), s(3)]);
    }

    #[test]
    fn linker_atom_is_a_linker() {
        assert_eq!(framework(&linked_rings()).role(s(7)), Some(Role::Linker));
    }

    #[test]
    fn linkers_are_listed() {
        let linkers: Vec<SiteId> = framework(&linked_rings()).linkers().collect();
        assert_eq!(linkers, vec![s(7)]);
    }

    #[test]
    fn side_chains_are_listed() {
        let side: Vec<SiteId> = framework(&ring_with_tail()).side_chains().collect();
        assert_eq!(side, vec![s(4)]);
    }

    #[test]
    fn framework_is_rings_and_linkers() {
        let mut sites: Vec<SiteId> = framework(&linked_rings()).sites().collect();
        sites.sort_unstable();
        assert_eq!(sites, vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7)]);
    }

    #[test]
    fn role_of_unknown_site_is_none() {
        assert!(framework(&triangle()).role(s(99)).is_none());
    }

    #[test]
    fn roles_partition_all_sites() {
        let m = linked_rings();
        let fw = framework(&m);
        for site in m.sites() {
            assert!(fw.role(site).is_some());
        }
        let scaffold: HashSet<SiteId> = fw.sites().collect();
        let side: HashSet<SiteId> = fw.side_chains().collect();
        assert!(scaffold.is_disjoint(&side));
        let all: HashSet<SiteId> = m.sites().collect();
        assert_eq!(&scaffold | &side, all);
    }
}

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

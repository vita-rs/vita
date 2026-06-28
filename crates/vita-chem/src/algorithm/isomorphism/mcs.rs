use std::collections::HashSet;

use vita_core::SiteId;

use super::indexed::{Adjacency, index};
use crate::utils::embeddings;
use crate::{BondId, HasBonds};

/// The maximum common substructure of two molecules: the largest connected set of
/// bonds shared by both, with the atom and bond correspondence it induces.
///
/// Obtain via [`mcs`].
#[derive(Debug, PartialEq, Eq)]
pub struct CommonSubgraph {
    sites: Vec<(SiteId, SiteId)>,
    bonds: Vec<(BondId, BondId)>,
}

impl CommonSubgraph {
    /// Number of shared bonds — the size the common substructure maximises.
    pub fn len(&self) -> usize {
        self.bonds.len()
    }

    /// Returns `true` if the molecules share no bond.
    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }

    /// Iterates the matched `(first site, second site)` pairs, ordered by the
    /// first molecule's site.
    pub fn sites(&self) -> impl Iterator<Item = (SiteId, SiteId)> + '_ {
        self.sites.iter().copied()
    }

    /// Iterates the matched `(first bond, second bond)` pairs, ordered by the
    /// first molecule's bond.
    pub fn bonds(&self) -> impl Iterator<Item = (BondId, BondId)> + '_ {
        self.bonds.iter().copied()
    }
}

/// Finds a maximum common substructure of `a` and `b`.
///
/// Returns the largest connected set of bonds shared by both molecules — the
/// common scaffold a chemist compares structures by — with the atom and bond
/// correspondence between the two. The match is the caller's to define.
///
/// The substructure is connected and maximises shared bonds; molecules with no
/// compatible bond share an empty one. When several substructures are maximal, a
/// deterministic one is returned.
///
/// # Complexity
///
/// Exponential in the worst case, as the problem is NP-hard; the search grows
/// only connected fragments and abandons any that fail to embed, staying
/// tractable for the molecules chemistry poses.
pub fn mcs<A, B>(
    a: &A,
    b: &B,
    site_match: impl Fn(SiteId, SiteId) -> bool,
    bond_match: impl Fn(BondId, BondId) -> bool,
) -> CommonSubgraph
where
    A: HasBonds,
    B: HasBonds,
{
    let a_index = index(a);
    let b_index = index(b);

    let swapped = b_index.bonds.len() < a_index.bonds.len();
    let (small, large) = if swapped {
        (&b_index, &a_index)
    } else {
        (&a_index, &b_index)
    };
    let site_ok = |s: SiteId, l: SiteId| {
        if swapped {
            site_match(l, s)
        } else {
            site_match(s, l)
        }
    };
    let bond_ok = |s: BondId, l: BondId| {
        if swapped {
            bond_match(l, s)
        } else {
            bond_match(s, l)
        }
    };

    let mut ends = vec![(0usize, 0usize); small.bonds.len()];
    for (u, incident) in small.adjacency.iter().enumerate() {
        for &(v, e) in incident {
            ends[e] = (u, v);
        }
    }

    let embeds = |seed: &[usize]| -> Option<Vec<usize>> {
        let (atoms, pattern) = fragment(seed, &ends);
        embeddings(
            &pattern,
            large.adjacency.clone(),
            |l, lv| site_ok(small.sites[atoms[l]], large.sites[lv]),
            |e, lb| bond_ok(small.bonds[e], large.bonds[lb]),
        )
        .next()
    };

    let cap = small.bonds.len().min(large.bonds.len());
    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut stack: Vec<Vec<usize>> = Vec::new();
    let mut best: Vec<usize> = Vec::new();

    for e in 0..small.bonds.len() {
        let seed = vec![e];
        if seen.insert(seed.clone()) && embeds(&seed).is_some() {
            if seed.len() > best.len() {
                best = seed.clone();
            }
            stack.push(seed);
        }
    }

    while let Some(seed) = stack.pop() {
        if best.len() == cap {
            break;
        }
        let atoms = atoms_of(&seed, &ends);
        for next in frontier(&seed, &atoms, &small.adjacency) {
            let grown = with(&seed, next);
            if seen.insert(grown.clone()) && embeds(&grown).is_some() {
                if grown.len() > best.len() {
                    best = grown.clone();
                }
                stack.push(grown);
            }
        }
    }

    let embedding = embeds(&best).unwrap_or_default();
    let atoms = atoms_of(&best, &ends);

    let mut sites: Vec<(SiteId, SiteId)> = atoms
        .iter()
        .enumerate()
        .map(|(l, &atom)| {
            let (s, t) = (small.sites[atom], large.sites[embedding[l]]);
            if swapped { (t, s) } else { (s, t) }
        })
        .collect();
    sites.sort_unstable();

    let mut bonds: Vec<(BondId, BondId)> = best
        .iter()
        .map(|&e| {
            let (u, v) = ends[e];
            let (lu, lv) = (
                atoms.binary_search(&u).unwrap(),
                atoms.binary_search(&v).unwrap(),
            );
            let (xu, xv) = (embedding[lu], embedding[lv]);
            let image = large.adjacency[xu]
                .iter()
                .find(|&&(neighbour, _)| neighbour == xv)
                .map(|&(_, bond)| bond)
                .unwrap();
            let (s, t) = (small.bonds[e], large.bonds[image]);
            if swapped { (t, s) } else { (s, t) }
        })
        .collect();
    bonds.sort_unstable();

    CommonSubgraph { sites, bonds }
}

/// The endpoint atoms of a bond set, in ascending order, and the adjacency of the
/// fragment they span over those atoms — the form the matching engine consumes.
fn fragment(seed: &[usize], ends: &[(usize, usize)]) -> (Vec<usize>, Adjacency) {
    let atoms = atoms_of(seed, ends);
    let mut adjacency: Adjacency = vec![Vec::new(); atoms.len()];
    for &e in seed {
        let (u, v) = ends[e];
        let (lu, lv) = (
            atoms.binary_search(&u).unwrap(),
            atoms.binary_search(&v).unwrap(),
        );
        adjacency[lu].push((lv, e));
        adjacency[lv].push((lu, e));
    }
    (atoms, adjacency)
}

/// The atoms a bond set touches, in ascending order.
fn atoms_of(seed: &[usize], ends: &[(usize, usize)]) -> Vec<usize> {
    let mut atoms: Vec<usize> = seed.iter().flat_map(|&e| [ends[e].0, ends[e].1]).collect();
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

/// The bonds incident to a fragment's `atoms` that it does not already hold — the
/// bonds it can grow by.
fn frontier(seed: &[usize], atoms: &[usize], adjacency: &Adjacency) -> Vec<usize> {
    let mut bonds: Vec<usize> = atoms
        .iter()
        .flat_map(|&atom| adjacency[atom].iter().map(|&(_, bond)| bond))
        .collect();
    bonds.sort_unstable();
    bonds.dedup();
    bonds.retain(|bond| seed.binary_search(bond).is_err());
    bonds
}

/// A sorted bond set with one more bond.
fn with(seed: &[usize], bond: usize) -> Vec<usize> {
    let mut grown = seed.to_vec();
    let at = grown.binary_search(&bond).unwrap_or_else(|at| at);
    grown.insert(at, bond);
    grown
}

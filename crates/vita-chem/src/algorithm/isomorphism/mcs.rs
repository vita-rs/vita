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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondOrder, HasBondOrders};
    use vita_core::{Element, HasElements, HasSites};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasElements for Mol {
        fn element(&self, site: SiteId) -> Element {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.elements[i]
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

    impl HasBondOrders for Mol {
        fn bond_order(&self, bond: BondId) -> BondOrder {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.orders[i]
        }
    }

    fn common(a: &Mol, b: &Mol) -> CommonSubgraph {
        mcs(
            a,
            b,
            |x, y| a.element(x) == b.element(y),
            |p, q| a.bond_order(p) == b.bond_order(q),
        )
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            elements: m.elements.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
            orders: m.orders.iter().rev().copied().collect(),
        }
    }

    fn empty() -> Mol {
        Mol {
            sites: vec![],
            elements: vec![],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        }
    }

    fn ethane() -> Mol {
        Mol {
            sites: vec![s(1), s(2)],
            elements: vec![elem("C"), elem("C")],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
            orders: vec![BondOrder::Single],
        }
    }

    fn ethanol() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("C"), elem("C"), elem("O")],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    fn ethanamine() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("C"), elem("C"), elem("N")],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    fn propanol() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            elements: vec![elem("C"), elem("C"), elem("C"), elem("O")],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4))],
            orders: vec![BondOrder::Single, BondOrder::Single, BondOrder::Single],
        }
    }

    fn propene() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("C"), elem("C"), elem("C")],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Double, BondOrder::Single],
        }
    }

    fn propane() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("C"), elem("C"), elem("C")],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    fn two_ethanes() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            elements: vec![elem("C"); 4],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(3), s(4))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    #[test]
    fn empty_molecule_shares_nothing() {
        let c = common(&empty(), &ethanol());
        assert_eq!(c.len(), 0);
        assert!(c.is_empty());
    }

    #[test]
    fn identical_molecules_share_everything() {
        let c = common(&ethanol(), &ethanol());
        assert_eq!(c.len(), 2);
        assert_eq!(c.sites().count(), 3);
    }

    #[test]
    fn smaller_embeds_whole() {
        let c = common(&ethane(), &ethanol());
        assert_eq!(c.len(), 1);
        assert_eq!(c.bonds().collect::<Vec<_>>(), vec![(b(1), b(1))]);
    }

    #[test]
    fn shared_chain_is_found() {
        let c = common(&ethanol(), &propanol());
        assert_eq!(c.len(), 2);
        assert_eq!(c.sites().count(), 3);
    }

    #[test]
    fn element_constrains_the_match() {
        let c = common(&ethanol(), &ethanamine());
        assert_eq!(c.len(), 1);
        assert_eq!(
            c.sites().collect::<Vec<_>>(),
            vec![(s(1), s(1)), (s(2), s(2))]
        );
    }

    #[test]
    fn bond_order_constrains_the_match() {
        let c = common(&propene(), &propane());
        assert_eq!(c.len(), 1);
    }

    #[test]
    fn the_common_substructure_is_connected() {
        let c = common(&two_ethanes(), &propane());
        assert_eq!(c.len(), 1);
    }

    #[test]
    fn disjoint_molecules_share_nothing() {
        let mut methane = ethane();
        methane.sites = vec![s(1)];
        methane.elements = vec![elem("C")];
        methane.bonds = vec![];
        methane.endpoints = vec![];
        methane.orders = vec![];
        let mut water = ethane();
        water.sites = vec![s(1)];
        water.elements = vec![elem("O")];
        water.bonds = vec![];
        water.endpoints = vec![];
        water.orders = vec![];
        assert!(common(&methane, &water).is_empty());
    }

    #[test]
    fn correspondence_is_consistent() {
        let c = common(&ethanol(), &propanol());
        for (x, y) in c.sites() {
            assert_eq!(ethanol().element(x), propanol().element(y));
        }
        for (p, q) in c.bonds() {
            assert_eq!(ethanol().bond_order(p), propanol().bond_order(q));
        }
    }

    #[test]
    fn common_substructure_is_independent_of_input_order() {
        assert_eq!(
            common(&ethanol(), &propanol()),
            common(&reversed(&ethanol()), &reversed(&propanol()))
        );
    }
}

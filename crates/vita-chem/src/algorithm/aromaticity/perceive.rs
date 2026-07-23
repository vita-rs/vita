use vita_core::{HasElements, SiteId};

use crate::algorithm::utils::{DisjointSet, FxHashMap, FxHashSet, SortedMap, electronegativity};
use crate::capability::delegation::forward_capabilities;
use crate::topology::ring::{Ring, RingMembership, rings};
use crate::valence::lone_pairs;
use crate::{
    BondId, BondOrder, HasAromaticity, HasBondOrders, HasBonds, HasFormalCharges,
    HasRadicalElectrons,
};

/// The aromatic π systems perceived in a molecule.
///
/// A bond is aromatic when it lies in a ring that satisfies Hückel's rule; a
/// site is aromatic when one of its bonds is. In biphenyl the two sites joined
/// by the central bond are aromatic while that bond is not.
///
/// Obtain via [`perceive`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Aromaticity {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
}

impl Aromaticity {
    /// Returns `true` if `site` lies in an aromatic system.
    ///
    /// Returns `false` if `site` is absent from the molecule or is not
    /// aromatic.
    pub fn contains_site(&self, site: SiteId) -> bool {
        self.sites.binary_search(&site).is_ok()
    }

    /// Returns `true` if `bond` is aromatic.
    ///
    /// Returns `false` if `bond` is absent from the molecule or is not
    /// aromatic.
    pub fn contains_bond(&self, bond: BondId) -> bool {
        self.bonds.binary_search(&bond).is_ok()
    }

    /// Number of aromatic sites.
    pub fn site_count(&self) -> usize {
        self.sites.len()
    }

    /// Number of aromatic bonds.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Iterates the aromatic sites in ascending order.
    pub fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.sites.iter().copied()
    }

    /// Iterates the aromatic bonds in ascending order.
    pub fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bonds.iter().copied()
    }

    /// Returns `true` if the molecule has no aromatic system.
    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }

    /// Binds this perception to `mol`, yielding a view that implements
    /// [`HasAromaticity`].
    ///
    /// The view borrows both, so `mol` stays immutable while it is held — the
    /// perception cannot silently fall out of step with the molecule it
    /// describes. Use it to feed a perceived molecule to anything that reads
    /// the [`HasAromaticity`] capability.
    pub fn bind<'a, M: HasBonds>(&'a self, mol: &'a M) -> WithAromaticity<'a, M> {
        WithAromaticity {
            mol,
            aromaticity: self,
        }
    }
}

/// A molecule viewed together with its perceived [`Aromaticity`].
///
/// Answers aromaticity from the perception and forwards every other core and
/// chem capability to the molecule, so a computed result reads as the
/// [`HasAromaticity`] capability its consumers expect — at no cost beyond the
/// two references it holds.
///
/// Obtain via [`Aromaticity::bind`].
pub struct WithAromaticity<'a, M> {
    mol: &'a M,
    aromaticity: &'a Aromaticity,
}

impl<M> Copy for WithAromaticity<'_, M> {}

impl<M> Clone for WithAromaticity<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

forward_capabilities!(
    WithAromaticity,
    mol,
    HasAccelerations,
    HasElements,
    HasIsotopes,
    HasLattice,
    HasMasses,
    HasNetCharge,
    HasPositions,
    HasSites,
    HasVelocities,
    HasBondOrders,
    HasBonds,
    HasFormalCharges,
    HasHybridizations,
    HasPartialCharges,
    HasRadicalElectrons,
    HasStereoConfigurations,
);

impl<M: HasBonds> HasAromaticity for WithAromaticity<'_, M> {
    fn is_aromatic(&self, bond: BondId) -> bool {
        self.aromaticity.contains_bond(bond)
    }

    fn is_aromatic_site(&self, site: SiteId) -> bool {
        self.aromaticity.contains_site(site)
    }
}

/// Perceives the aromatic bonds and sites of a molecule.
///
/// A cycle is aromatic when it is conjugated all the way round and holds 4*n*+2
/// π electrons (Hückel's rule). Each ring atom donates a fixed number of
/// electrons to the perpendicular π system, read from its localized bonding:
///
/// - an endocyclic double bond donates one electron;
/// - an exocyclic double bond to a more electronegative atom donates none, the
///   pair drawn off the ring to leave an empty p orbital (the carbonyl carbon
///   of tropone); to an equal or less electronegative one it keeps the pair and
///   excludes the cycle (the methylene carbon of heptafulvene);
/// - an otherwise saturated atom with a lone pair donates the pair, two
///   electrons (the nitrogen of pyrrole, the oxygen of furan);
/// - an electron-deficient atom donates none through its empty p orbital (the
///   cation of tropylium, the boron of borole);
/// - a saturated sp³ atom, a π radical, or one of undefined valence cannot join
///   an aromatic system, and excludes any cycle running through it.
///
/// The cycles tested are the rings of the [minimum cycle basis](rings) and, for
/// each fused ring system, its outer perimeter — enough to recognize azulene,
/// whose aromaticity lives on the ten-membered rim and not on either ring.
///
/// Aromaticity is a property of the localized structure: a molecule already
/// carrying [`Aromatic`](BondOrder::Aromatic) bonds leaves the electron count
/// undefined and its atoms non-aromatic, just as [`valence`](crate::valence)
/// and [`lone_pairs`] report `None`. Kekulize first, then perceive.
///
/// # Complexity
///
/// O(V · E³ / w) time and O(V · E² / w) space, over the molecule's `V` sites and
/// `E` bonds for word width `w` = 64, dominated by the [minimum cycle
/// basis](rings).
pub fn perceive<M>(mol: &M) -> Aromaticity
where
    M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons,
{
    let rings = rings(mol);
    if rings.is_empty() {
        return Aromaticity {
            sites: Vec::new(),
            bonds: Vec::new(),
        };
    }

    // The π contribution of every ring atom, settled once from its bonding.
    let membership = rings.membership();
    let pi: SortedMap<SiteId, Option<u32>> = SortedMap::from_pairs(
        membership
            .sites()
            .map(|site| (site, contribution(mol, site, &membership))),
    );

    let basis: Vec<&Ring> = rings.iter().collect();
    let mut aromatic: FxHashSet<BondId> = FxHashSet::default();

    // Every basis ring is a candidate aromatic cycle.
    for ring in &basis {
        if huckel(ring.sites(), &pi) {
            aromatic.extend(ring.bonds().iter().copied());
        }
    }

    // So is the outer perimeter of each fused ring system. When the rim is
    // aromatic the whole system is, down to the bonds buried inside it: the rim
    // of azulene carries ten electrons, which makes its shared bond aromatic
    // although neither of its rings is on its own.
    for group in fused_systems(&basis) {
        if group.len() < 2 {
            continue;
        }
        if perimeter(mol, &basis, &group).is_some_and(|rim| huckel(&rim, &pi)) {
            for &ring in &group {
                aromatic.extend(basis[ring].bonds().iter().copied());
            }
        }
    }

    // A site is aromatic exactly when one of its bonds is.
    let mut sites: FxHashSet<SiteId> = FxHashSet::default();
    for &bond in &aromatic {
        let (a, b) = mol.bond_endpoints(bond);
        sites.insert(a);
        sites.insert(b);
    }

    let mut sites: Vec<SiteId> = sites.into_iter().collect();
    let mut bonds: Vec<BondId> = aromatic.into_iter().collect();
    sites.sort_unstable();
    bonds.sort_unstable();
    Aromaticity { sites, bonds }
}

/// Returns whether a cycle satisfies Hückel's rule: every atom donates to the
/// π system and the donated electrons number 4*n*+2.
///
/// An atom that cannot belong to an aromatic system (a `None` contribution)
/// disqualifies the whole cycle.
fn huckel(sites: &[SiteId], pi: &SortedMap<SiteId, Option<u32>>) -> bool {
    let mut electrons = 0;
    for site in sites {
        match pi.get(site).copied().flatten() {
            Some(donated) => electrons += donated,
            None => return false,
        }
    }
    electrons >= 2 && (electrons - 2) % 4 == 0
}

/// The π-electron contribution of a ring atom, or `None` when the atom cannot
/// belong to an aromatic system.
///
/// `rings` carries the ring-bond membership, used to tell an endocyclic double
/// bond (a ring bond) from an exocyclic one (a bridge, as in the carbonyl of
/// tropone or the methylene of heptafulvene).
fn contribution<M>(mol: &M, site: SiteId, rings: &RingMembership) -> Option<u32>
where
    M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons,
{
    let mut endocyclic = 0;
    let mut exocyclic = 0;
    let mut across = None;
    for (bond, neighbor) in mol.bonds_of(site) {
        let pi = match mol.bond_order(bond) {
            BondOrder::Single => 0,
            BondOrder::Double => 1,
            BondOrder::Triple => 2,
            // Delocalized input, or a metal–metal bond: not a localized π system.
            _ => return None,
        };
        if pi > 0 {
            if rings.contains_bond(bond) {
                endocyclic += pi;
            } else {
                exocyclic += pi;
                across = Some(neighbor);
            }
        }
    }

    // An aromatic atom carries one π bond at most; a second leaves it sp, with
    // no p orbital perpendicular to the ring.
    if endocyclic + exocyclic > 1 {
        return None;
    }
    if endocyclic == 1 {
        return Some(1); // one p electron, shared around the ring
    }
    if let Some(partner) = across {
        // The exocyclic π leaves an empty p orbital perpendicular to the ring
        // only when it is polarized toward a more electronegative atom, drawing
        // the pair off the ring (the carbonyl carbon of tropone); an unpolarized
        // bond keeps the pair localized within it (the methylene of heptafulvene).
        let polarized =
            electronegativity(mol.element(partner))? > electronegativity(mol.element(site))?;
        return polarized.then_some(0);
    }

    // No π bond: the perpendicular p orbital holds a lone pair, sits empty, or
    // is not there at all.
    if lone_pairs(mol, site)? >= 1 {
        return Some(2); // a lone pair joins the ring π system
    }
    if mol.radical_electron(site) >= 1 {
        return None; // an unpaired electron, not a closed-shell donor
    }
    if mol.degree(site) <= 3 {
        return Some(0); // an empty p orbital on an electron-deficient center
    }
    None // four σ bonds: saturated, with no p orbital to offer
}

/// Groups the basis rings into fused systems: rings sharing a site coalesce
/// into one group. Returns the ring indices of each group.
fn fused_systems(basis: &[&Ring]) -> Vec<Vec<usize>> {
    let mut sets = DisjointSet::new(basis.len());
    let mut owner: FxHashMap<SiteId, usize> = FxHashMap::default();
    for (i, ring) in basis.iter().enumerate() {
        for &site in ring.sites() {
            match owner.get(&site) {
                Some(&j) => {
                    sets.union(i, j);
                }
                None => {
                    owner.insert(site, i);
                }
            }
        }
    }

    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for i in 0..basis.len() {
        groups.entry(sets.find(i)).or_default().push(i);
    }
    groups.into_values().collect()
}

/// The sites on the outer perimeter of a fused ring system, in no particular
/// order — the rim drawn by the bonds lying in an odd number of its rings.
///
/// Returns `None` unless the rim is a single simple cycle: a spiro junction or
/// bridged cage has a boundary that branches or splits, no annulene to test.
fn perimeter<M: HasBondOrders>(mol: &M, basis: &[&Ring], group: &[usize]) -> Option<Vec<SiteId>> {
    // Symmetric difference of the rings' bonds: shared (internal) bonds cancel.
    let mut shared: FxHashMap<BondId, u32> = FxHashMap::default();
    for &ring in group {
        for &bond in basis[ring].bonds() {
            *shared.entry(bond).or_insert(0) += 1;
        }
    }
    let bonds: Vec<BondId> = shared
        .into_iter()
        .filter(|&(_, count)| count % 2 == 1)
        .map(|(bond, _)| bond)
        .collect();
    if bonds.is_empty() {
        return None;
    }

    // A single simple cycle: every site meets exactly two perimeter bonds, and
    // the bonds are all reachable from one another.
    let mut adjacency: FxHashMap<SiteId, Vec<SiteId>> = FxHashMap::default();
    for &bond in &bonds {
        let (u, v) = mol.bond_endpoints(bond);
        adjacency.entry(u).or_default().push(v);
        adjacency.entry(v).or_default().push(u);
    }
    if adjacency.values().any(|neighbors| neighbors.len() != 2) {
        return None;
    }

    let sites: Vec<SiteId> = adjacency.keys().copied().collect();
    let mut seen: FxHashSet<SiteId> = FxHashSet::default();
    let mut stack = vec![sites[0]];
    while let Some(site) = stack.pop() {
        if seen.insert(site) {
            stack.extend(adjacency[&site].iter().copied());
        }
    }
    if seen.len() != sites.len() {
        return None; // the boundary is several disjoint cycles
    }

    Some(sites)
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasSites};

    use crate::BondOrder::{Aromatic, Double, Single};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
        charges: Vec<i8>,
        radicals: Vec<u8>,
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

    impl HasFormalCharges for Mol {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.charges[i]
        }
    }

    impl HasRadicalElectrons for Mol {
        fn radical_electron(&self, site: SiteId) -> u8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.radicals[i]
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

    fn mol(atoms: &[(u32, &str, i8, u8)], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, ..)| s(id)).collect(),
            elements: atoms
                .iter()
                .map(|&(_, sym, ..)| Element::from_symbol(sym).unwrap())
                .collect(),
            charges: atoms.iter().map(|&(_, _, charge, _)| charge).collect(),
            radicals: atoms.iter().map(|&(_, _, _, radical)| radical).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c, _)| (s(a), s(c))).collect(),
            orders: bonds.iter().map(|&(_, _, _, order)| order).collect(),
        }
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            elements: m.elements.iter().rev().copied().collect(),
            charges: m.charges.iter().rev().copied().collect(),
            radicals: m.radicals.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
            orders: m.orders.iter().rev().copied().collect(),
        }
    }

    fn empty() -> Mol {
        mol(&[], &[])
    }

    fn ethane() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
                (5, 2, 6, Single),
                (6, 2, 7, Single),
                (7, 2, 8, Single),
            ],
        )
    }

    fn benzene() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Single),
                (3, 3, 4, Double),
                (4, 4, 5, Single),
                (5, 5, 6, Double),
                (6, 6, 1, Single),
                (7, 1, 7, Single),
                (8, 2, 8, Single),
                (9, 3, 9, Single),
                (10, 4, 10, Single),
                (11, 5, 11, Single),
                (12, 6, 12, Single),
            ],
        )
    }

    fn delocalized_benzene() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 6, Aromatic),
                (6, 6, 1, Aromatic),
                (7, 1, 7, Single),
                (8, 2, 8, Single),
                (9, 3, 9, Single),
                (10, 4, 10, Single),
                (11, 5, 11, Single),
                (12, 6, 12, Single),
            ],
        )
    }

    fn pyridine() -> Mol {
        mol(
            &[
                (1, "N", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Single),
                (3, 3, 4, Double),
                (4, 4, 5, Single),
                (5, 5, 6, Double),
                (6, 6, 1, Single),
                (7, 2, 7, Single),
                (8, 3, 8, Single),
                (9, 4, 9, Single),
                (10, 5, 10, Single),
                (11, 6, 11, Single),
            ],
        )
    }

    fn pyrrole() -> Mol {
        mol(
            &[
                (1, "N", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Double),
                (3, 3, 4, Single),
                (4, 4, 5, Double),
                (5, 5, 1, Single),
                (6, 1, 6, Single),
                (7, 2, 7, Single),
                (8, 3, 8, Single),
                (9, 4, 9, Single),
                (10, 5, 10, Single),
            ],
        )
    }

    fn imidazole() -> Mol {
        mol(
            &[
                (1, "N", 0, 0),
                (2, "C", 0, 0),
                (3, "N", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Double),
                (3, 3, 4, Single),
                (4, 4, 5, Double),
                (5, 5, 1, Single),
                (6, 1, 6, Single),
                (7, 2, 7, Single),
                (8, 4, 8, Single),
                (9, 5, 9, Single),
            ],
        )
    }

    fn tropylium_cation() -> Mol {
        mol(
            &[
                (1, "C", 1, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "C", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Double),
                (3, 3, 4, Single),
                (4, 4, 5, Double),
                (5, 5, 6, Single),
                (6, 6, 7, Double),
                (7, 7, 1, Single),
                (8, 1, 8, Single),
                (9, 2, 9, Single),
                (10, 3, 10, Single),
                (11, 4, 11, Single),
                (12, 5, 12, Single),
                (13, 6, 13, Single),
                (14, 7, 14, Single),
            ],
        )
    }

    fn cyclopentadienyl_anion() -> Mol {
        mol(
            &[
                (1, "C", -1, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Double),
                (3, 3, 4, Single),
                (4, 4, 5, Double),
                (5, 5, 1, Single),
                (6, 1, 6, Single),
                (7, 2, 7, Single),
                (8, 3, 8, Single),
                (9, 4, 9, Single),
                (10, 5, 10, Single),
            ],
        )
    }

    fn tropone() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "C", 0, 0),
                (8, "O", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Double),
                (3, 3, 4, Single),
                (4, 4, 5, Double),
                (5, 5, 6, Single),
                (6, 6, 7, Double),
                (7, 7, 1, Single),
                (8, 1, 8, Double),
                (9, 2, 9, Single),
                (10, 3, 10, Single),
                (11, 4, 11, Single),
                (12, 5, 12, Single),
                (13, 6, 13, Single),
                (14, 7, 14, Single),
            ],
        )
    }

    fn cyclohexane() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
                (15, "H", 0, 0),
                (16, "H", 0, 0),
                (17, "H", 0, 0),
                (18, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Single),
                (3, 3, 4, Single),
                (4, 4, 5, Single),
                (5, 5, 6, Single),
                (6, 6, 1, Single),
                (7, 1, 7, Single),
                (8, 1, 8, Single),
                (9, 2, 9, Single),
                (10, 2, 10, Single),
                (11, 3, 11, Single),
                (12, 3, 12, Single),
                (13, 4, 13, Single),
                (14, 4, 14, Single),
                (15, 5, 15, Single),
                (16, 5, 16, Single),
                (17, 6, 17, Single),
                (18, 6, 18, Single),
            ],
        )
    }

    fn cyclobutadiene() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Single),
                (3, 3, 4, Double),
                (4, 4, 1, Single),
                (5, 1, 5, Single),
                (6, 2, 6, Single),
                (7, 3, 7, Single),
                (8, 4, 8, Single),
            ],
        )
    }

    fn cyclopentadiene() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Double),
                (3, 3, 4, Single),
                (4, 4, 5, Double),
                (5, 5, 1, Single),
                (6, 1, 6, Single),
                (7, 1, 7, Single),
                (8, 2, 8, Single),
                (9, 3, 9, Single),
                (10, 4, 10, Single),
                (11, 5, 11, Single),
            ],
        )
    }

    fn cycloheptatrienyl_radical() -> Mol {
        mol(
            &[
                (1, "C", 0, 1),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "C", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Double),
                (3, 3, 4, Single),
                (4, 4, 5, Double),
                (5, 5, 6, Single),
                (6, 6, 7, Double),
                (7, 7, 1, Single),
                (8, 1, 8, Single),
                (9, 2, 9, Single),
                (10, 3, 10, Single),
                (11, 4, 11, Single),
                (12, 5, 12, Single),
                (13, 6, 13, Single),
                (14, 7, 14, Single),
            ],
        )
    }

    fn naphthalene() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "C", 0, 0),
                (8, "C", 0, 0),
                (9, "C", 0, 0),
                (10, "C", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
                (15, "H", 0, 0),
                (16, "H", 0, 0),
                (17, "H", 0, 0),
                (18, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Single),
                (3, 3, 4, Double),
                (4, 4, 9, Single),
                (5, 9, 10, Double),
                (6, 10, 1, Single),
                (7, 5, 6, Double),
                (8, 6, 7, Single),
                (9, 7, 8, Double),
                (10, 8, 9, Single),
                (11, 10, 5, Single),
                (12, 1, 11, Single),
                (13, 2, 12, Single),
                (14, 3, 13, Single),
                (15, 4, 14, Single),
                (16, 5, 15, Single),
                (17, 6, 16, Single),
                (18, 7, 17, Single),
                (19, 8, 18, Single),
            ],
        )
    }

    fn azulene() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "C", 0, 0),
                (8, "C", 0, 0),
                (9, "C", 0, 0),
                (10, "C", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
                (15, "H", 0, 0),
                (16, "H", 0, 0),
                (17, "H", 0, 0),
                (18, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Double),
                (3, 3, 4, Single),
                (4, 4, 5, Double),
                (5, 5, 1, Single),
                (6, 2, 6, Single),
                (7, 6, 7, Double),
                (8, 7, 8, Single),
                (9, 8, 9, Double),
                (10, 9, 10, Single),
                (11, 10, 1, Double),
                (12, 3, 11, Single),
                (13, 4, 12, Single),
                (14, 5, 13, Single),
                (15, 6, 14, Single),
                (16, 7, 15, Single),
                (17, 8, 16, Single),
                (18, 9, 17, Single),
                (19, 10, 18, Single),
            ],
        )
    }

    fn biphenyl() -> Mol {
        mol(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "C", 0, 0),
                (8, "C", 0, 0),
                (9, "C", 0, 0),
                (10, "C", 0, 0),
                (11, "C", 0, 0),
                (12, "C", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
                (15, "H", 0, 0),
                (16, "H", 0, 0),
                (17, "H", 0, 0),
                (18, "H", 0, 0),
                (19, "H", 0, 0),
                (20, "H", 0, 0),
                (21, "H", 0, 0),
                (22, "H", 0, 0),
            ],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Single),
                (3, 3, 4, Double),
                (4, 4, 5, Single),
                (5, 5, 6, Double),
                (6, 6, 1, Single),
                (7, 7, 8, Double),
                (8, 8, 9, Single),
                (9, 9, 10, Double),
                (10, 10, 11, Single),
                (11, 11, 12, Double),
                (12, 12, 7, Single),
                (13, 1, 7, Single),
                (14, 2, 13, Single),
                (15, 3, 14, Single),
                (16, 4, 15, Single),
                (17, 5, 16, Single),
                (18, 6, 17, Single),
                (19, 8, 18, Single),
                (20, 9, 19, Single),
                (21, 10, 20, Single),
                (22, 11, 21, Single),
                (23, 12, 22, Single),
            ],
        )
    }

    #[test]
    fn empty_molecule_is_not_aromatic() {
        assert!(perceive(&empty()).is_empty());
    }

    #[test]
    fn acyclic_molecule_is_not_aromatic() {
        assert!(perceive(&ethane()).is_empty());
    }

    #[test]
    fn benzene_is_aromatic() {
        let aromaticity = perceive(&benzene());
        assert!(!aromaticity.is_empty());
        assert!((1..=6).all(|i| aromaticity.contains_bond(b(i))));
        assert!((1..=6).all(|i| aromaticity.contains_site(s(i))));
    }

    #[test]
    fn pyridine_is_aromatic() {
        assert!((1..=6).all(|i| perceive(&pyridine()).contains_bond(b(i))));
    }

    #[test]
    fn pyrrole_is_aromatic_through_its_nitrogen_lone_pair() {
        let aromaticity = perceive(&pyrrole());
        assert!((1..=5).all(|i| aromaticity.contains_bond(b(i))));
        assert!(aromaticity.contains_site(s(1)));
    }

    #[test]
    fn tropylium_cation_is_aromatic_through_its_empty_orbital() {
        let aromaticity = perceive(&tropylium_cation());
        assert!((1..=7).all(|i| aromaticity.contains_bond(b(i))));
        assert!(aromaticity.contains_site(s(1)));
    }

    #[test]
    fn cyclopentadienyl_anion_is_aromatic_through_its_carbanion() {
        let aromaticity = perceive(&cyclopentadienyl_anion());
        assert!((1..=5).all(|i| aromaticity.contains_bond(b(i))));
        assert!(aromaticity.contains_site(s(1)));
    }

    #[test]
    fn tropone_is_aromatic_through_its_exocyclic_carbonyl() {
        let aromaticity = perceive(&tropone());
        assert!((1..=7).all(|i| aromaticity.contains_bond(b(i))));
        assert!(!aromaticity.contains_bond(b(8)));
    }

    #[test]
    fn saturated_ring_is_not_aromatic() {
        assert!(perceive(&cyclohexane()).is_empty());
    }

    #[test]
    fn four_electron_ring_is_not_aromatic() {
        assert!(perceive(&cyclobutadiene()).is_empty());
    }

    #[test]
    fn interrupted_conjugation_is_not_aromatic() {
        assert!(perceive(&cyclopentadiene()).is_empty());
    }

    #[test]
    fn a_ring_radical_excludes_the_ring() {
        assert!(perceive(&cycloheptatrienyl_radical()).is_empty());
    }

    #[test]
    fn delocalized_input_is_not_perceived() {
        assert!(perceive(&delocalized_benzene()).is_empty());
    }

    #[test]
    fn non_aromatic_bond_is_not_reported() {
        assert!(!perceive(&benzene()).contains_bond(b(7)));
    }

    #[test]
    fn unknown_bond_is_not_aromatic() {
        assert!(!perceive(&benzene()).contains_bond(b(99)));
    }

    #[test]
    fn unknown_site_is_not_aromatic() {
        assert!(!perceive(&benzene()).contains_site(s(99)));
    }

    #[test]
    fn azulene_is_aromatic_on_its_perimeter() {
        let aromaticity = perceive(&azulene());
        assert!(!aromaticity.is_empty());
        assert!(aromaticity.contains_bond(b(1)));
    }

    #[test]
    fn naphthalene_is_aromatic_across_both_rings() {
        assert!((1..=11).all(|i| perceive(&naphthalene()).contains_bond(b(i))));
    }

    #[test]
    fn biphenyl_rings_are_aromatic_but_the_link_is_not() {
        let aromaticity = perceive(&biphenyl());
        assert!(aromaticity.contains_site(s(1)));
        assert!(aromaticity.contains_site(s(7)));
        assert!(!aromaticity.contains_bond(b(13)));
    }

    #[test]
    fn imidazole_is_aromatic_with_both_nitrogen_kinds() {
        assert!((1..=5).all(|i| perceive(&imidazole()).contains_bond(b(i))));
    }

    #[test]
    fn sites_and_bonds_are_listed_in_ascending_order() {
        let aromaticity = perceive(&benzene());
        let sites: Vec<SiteId> = aromaticity.sites().collect();
        let bonds: Vec<BondId> = aromaticity.bonds().collect();
        assert!(sites.windows(2).all(|w| w[0] < w[1]));
        assert!(bonds.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn counts_report_the_number_of_aromatic_sites_and_bonds() {
        let aromaticity = perceive(&benzene());
        assert_eq!(aromaticity.site_count(), 6);
        assert_eq!(aromaticity.bond_count(), 6);
    }

    #[test]
    fn bound_view_answers_the_aromaticity_capability() {
        let mol = benzene();
        let aromaticity = perceive(&mol);
        let view = aromaticity.bind(&mol);
        assert!(view.is_aromatic(b(1)));
        assert!(!view.is_aromatic(b(7)));
        assert!(view.is_aromatic_site(s(1)));
        assert!(!view.is_aromatic_site(s(7)));
    }

    #[test]
    fn bound_view_forwards_the_skeleton() {
        let mol = benzene();
        let aromaticity = perceive(&mol);
        let view = aromaticity.bind(&mol);
        assert_eq!(view.element(s(1)), Element::from_symbol("C").unwrap());
        assert_eq!(view.bond_endpoints(b(1)), mol.bond_endpoints(b(1)));
        assert_eq!(view.bond_count(), mol.bond_count());
    }

    #[test]
    fn perception_is_independent_of_input_order() {
        assert_eq!(perceive(&benzene()), perceive(&reversed(&benzene())));
    }
}

use std::collections::{HashMap, HashSet};

use vita_core::{HasElements, SiteId};

use crate::capability::delegation::forward_capabilities;
use crate::topology::ring::{Ring, RingMembership, rings};
use crate::utils::electronegativity;
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
pub struct Aromaticity {
    sites: HashSet<SiteId>,
    bonds: HashSet<BondId>,
}

impl Aromaticity {
    /// Returns `true` if `site` lies in an aromatic system.
    ///
    /// Returns `false` if `site` is absent from the molecule or is not
    /// aromatic.
    pub fn site(&self, site: SiteId) -> bool {
        self.sites.contains(&site)
    }

    /// Returns `true` if `bond` is aromatic.
    ///
    /// Returns `false` if `bond` is absent from the molecule or is not
    /// aromatic.
    pub fn bond(&self, bond: BondId) -> bool {
        self.bonds.contains(&bond)
    }

    /// Iterates the sites that lie in an aromatic system.
    pub fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.sites.iter().copied()
    }

    /// Iterates the aromatic bonds.
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
);

impl<M: HasBonds> HasAromaticity for WithAromaticity<'_, M> {
    fn is_aromatic(&self, bond: BondId) -> bool {
        self.aromaticity.bond(bond)
    }

    fn is_aromatic_site(&self, site: SiteId) -> bool {
        self.aromaticity.site(site)
    }
}

/// Perceives the aromatic bonds and sites of a molecule.
///
/// A cycle is aromatic when it is conjugated all the way round and holds 4*n*+2
/// π electrons (Hückel's rule). Each ring atom donates a fixed number of
/// electrons to the perpendicular π system, read from its localised bonding:
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
/// each fused ring system, its outer perimeter — enough to recognise azulene,
/// whose aromaticity lives on the ten-membered rim and not on either ring.
///
/// Aromaticity is a property of the localised structure: a molecule already
/// carrying [`Aromatic`](BondOrder::Aromatic) bonds leaves the electron count
/// undefined and its atoms non-aromatic, just as [`valence`](crate::valence)
/// and [`lone_pairs`] report `None`. Kekulise first, then perceive.
///
/// # Complexity
///
/// O(V² · E), dominated by the minimum cycle basis.
pub fn perceive<M>(mol: &M) -> Aromaticity
where
    M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons,
{
    let rings = rings(mol);
    if rings.is_empty() {
        return Aromaticity {
            sites: HashSet::new(),
            bonds: HashSet::new(),
        };
    }

    // The π contribution of every ring atom, settled once from its bonding.
    let membership = rings.membership();
    let pi: HashMap<SiteId, Option<u32>> = membership
        .sites()
        .map(|site| (site, contribution(mol, site, &membership)))
        .collect();

    let basis: Vec<&Ring> = rings.iter().collect();
    let mut aromatic: HashSet<BondId> = HashSet::new();

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
    let mut sites: HashSet<SiteId> = HashSet::new();
    for &bond in &aromatic {
        let (a, b) = mol.bond_endpoints(bond);
        sites.insert(a);
        sites.insert(b);
    }

    Aromaticity {
        sites,
        bonds: aromatic,
    }
}

/// Returns whether a cycle satisfies Hückel's rule: every atom donates to the
/// π system and the donated electrons number 4*n*+2.
///
/// An atom that cannot belong to an aromatic system (a `None` contribution)
/// disqualifies the whole cycle.
fn huckel(sites: &[SiteId], pi: &HashMap<SiteId, Option<u32>>) -> bool {
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
    for (bond, neighbour) in mol.bonds_of(site) {
        let pi = match mol.bond_order(bond) {
            BondOrder::Single => 0,
            BondOrder::Double => 1,
            BondOrder::Triple => 2,
            // Delocalised input, or a metal–metal bond: not a localised π system.
            _ => return None,
        };
        if pi > 0 {
            if rings.bond(bond) {
                endocyclic += pi;
            } else {
                exocyclic += pi;
                across = Some(neighbour);
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
        // only when it is polarised toward a more electronegative atom, drawing
        // the pair off the ring (the carbonyl carbon of tropone); an unpolarised
        // bond keeps the pair localised within it (the methylene of heptafulvene).
        let polarised =
            electronegativity(mol.element(partner))? > electronegativity(mol.element(site))?;
        return polarised.then_some(0);
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
        return Some(0); // an empty p orbital on an electron-deficient centre
    }
    None // four σ bonds: saturated, with no p orbital to offer
}

/// Groups the basis rings into fused systems: rings sharing a site coalesce
/// into one group. Returns the ring indices of each group.
fn fused_systems(basis: &[&Ring]) -> Vec<Vec<usize>> {
    let k = basis.len();
    let mut parent: Vec<usize> = (0..k).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    let mut owner: HashMap<SiteId, usize> = HashMap::new();
    for (i, ring) in basis.iter().enumerate() {
        for &site in ring.sites() {
            match owner.get(&site) {
                Some(&j) => {
                    let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                    if a != b {
                        parent[a] = b;
                    }
                }
                None => {
                    owner.insert(site, i);
                }
            }
        }
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..k {
        let root = find(&mut parent, i);
        groups.entry(root).or_default().push(i);
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
    let mut shared: HashMap<BondId, u32> = HashMap::new();
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
    let mut adjacency: HashMap<SiteId, Vec<SiteId>> = HashMap::new();
    for &bond in &bonds {
        let (u, v) = mol.bond_endpoints(bond);
        adjacency.entry(u).or_default().push(v);
        adjacency.entry(v).or_default().push(u);
    }
    if adjacency.values().any(|neighbours| neighbours.len() != 2) {
        return None;
    }

    let sites: Vec<SiteId> = adjacency.keys().copied().collect();
    let mut seen: HashSet<SiteId> = HashSet::new();
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
    use crate::{BondId, HasAromaticity, HasBonds};
    use vita_core::{Element, HasSites};

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
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.elements[i]
        }
    }

    impl HasFormalCharges for Mol {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.charges[i]
        }
    }

    impl HasRadicalElectrons for Mol {
        fn radical_electron(&self, site: SiteId) -> u8 {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
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

    fn empty() -> Mol {
        Mol {
            sites: vec![],
            elements: vec![],
            charges: vec![],
            radicals: vec![],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        }
    }

    fn ethane() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8)],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7)],
            endpoints: vec![
                (s(1), s(2)),
                (s(1), s(3)),
                (s(1), s(4)),
                (s(1), s(5)),
                (s(2), s(6)),
                (s(2), s(7)),
                (s(2), s(8)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn benzene() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
                (s(1), s(7)),
                (s(2), s(8)),
                (s(3), s(9)),
                (s(4), s(10)),
                (s(5), s(11)),
                (s(6), s(12)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn pyridine() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
            ],
            elements: vec![
                elem("N"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
                (s(2), s(7)),
                (s(3), s(8)),
                (s(4), s(9)),
                (s(5), s(10)),
                (s(6), s(11)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn pyrimidine() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9), s(10)],
            elements: vec![
                elem("N"),
                elem("C"),
                elem("N"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), b(10)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
                (s(2), s(7)),
                (s(4), s(8)),
                (s(5), s(9)),
                (s(6), s(10)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn pyrrole() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9), s(10)],
            elements: vec![
                elem("N"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), b(10)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(1)),
                (s(1), s(6)),
                (s(2), s(7)),
                (s(3), s(8)),
                (s(4), s(9)),
                (s(5), s(10)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn furan() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9)],
            elements: vec![
                elem("O"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(1)),
                (s(2), s(6)),
                (s(3), s(7)),
                (s(4), s(8)),
                (s(5), s(9)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn thiophene() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9)],
            elements: vec![
                elem("S"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(1)),
                (s(2), s(6)),
                (s(3), s(7)),
                (s(4), s(8)),
                (s(5), s(9)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn imidazole() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9)],
            elements: vec![
                elem("N"),
                elem("C"),
                elem("N"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(1)),
                (s(1), s(6)),
                (s(2), s(7)),
                (s(4), s(8)),
                (s(5), s(9)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn cyclopentadienyl_anion() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9), s(10)],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![-1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), b(10)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(1)),
                (s(1), s(6)),
                (s(2), s(7)),
                (s(3), s(8)),
                (s(4), s(9)),
                (s(5), s(10)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn tropylium_cation() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(1)),
                (s(1), s(8)),
                (s(2), s(9)),
                (s(3), s(10)),
                (s(4), s(11)),
                (s(5), s(12)),
                (s(6), s(13)),
                (s(7), s(14)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn phenyl_radical() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
                (s(2), s(7)),
                (s(3), s(8)),
                (s(4), s(9)),
                (s(5), s(10)),
                (s(6), s(11)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn tropone() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("O"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(1)),
                (s(1), s(8)),
                (s(2), s(9)),
                (s(3), s(10)),
                (s(4), s(11)),
                (s(5), s(12)),
                (s(6), s(13)),
                (s(7), s(14)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn naphthalene() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
                s(15),
                s(16),
                s(17),
                s(18),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
                b(15),
                b(16),
                b(17),
                b(18),
                b(19),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(9)),
                (s(9), s(10)),
                (s(10), s(1)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(8)),
                (s(8), s(9)),
                (s(10), s(5)),
                (s(1), s(11)),
                (s(2), s(12)),
                (s(3), s(13)),
                (s(4), s(14)),
                (s(5), s(15)),
                (s(6), s(16)),
                (s(7), s(17)),
                (s(8), s(18)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn indole() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
                s(15),
                s(16),
            ],
            elements: vec![
                elem("N"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
                b(15),
                b(16),
                b(17),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(9)),
                (s(9), s(1)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(8)),
                (s(8), s(9)),
                (s(1), s(10)),
                (s(2), s(11)),
                (s(3), s(12)),
                (s(5), s(13)),
                (s(6), s(14)),
                (s(7), s(15)),
                (s(8), s(16)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn azulene() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
                s(15),
                s(16),
                s(17),
                s(18),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
                b(15),
                b(16),
                b(17),
                b(18),
                b(19),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(1)),
                (s(2), s(6)),
                (s(6), s(7)),
                (s(7), s(8)),
                (s(8), s(9)),
                (s(9), s(10)),
                (s(10), s(1)),
                (s(3), s(11)),
                (s(4), s(12)),
                (s(5), s(13)),
                (s(6), s(14)),
                (s(7), s(15)),
                (s(8), s(16)),
                (s(9), s(17)),
                (s(10), s(18)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn biphenyl() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
                s(15),
                s(16),
                s(17),
                s(18),
                s(19),
                s(20),
                s(21),
                s(22),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ],
            radicals: vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
                b(15),
                b(16),
                b(17),
                b(18),
                b(19),
                b(20),
                b(21),
                b(22),
                b(23),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
                (s(7), s(8)),
                (s(8), s(9)),
                (s(9), s(10)),
                (s(10), s(11)),
                (s(11), s(12)),
                (s(12), s(7)),
                (s(1), s(7)),
                (s(2), s(13)),
                (s(3), s(14)),
                (s(4), s(15)),
                (s(5), s(16)),
                (s(6), s(17)),
                (s(8), s(18)),
                (s(9), s(19)),
                (s(10), s(20)),
                (s(11), s(21)),
                (s(12), s(22)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn cyclohexane() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
                s(15),
                s(16),
                s(17),
                s(18),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
                b(15),
                b(16),
                b(17),
                b(18),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
                (s(1), s(7)),
                (s(1), s(8)),
                (s(2), s(9)),
                (s(2), s(10)),
                (s(3), s(11)),
                (s(3), s(12)),
                (s(4), s(13)),
                (s(4), s(14)),
                (s(5), s(15)),
                (s(5), s(16)),
                (s(6), s(17)),
                (s(6), s(18)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn cyclopentadiene() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(1)),
                (s(1), s(6)),
                (s(1), s(7)),
                (s(2), s(8)),
                (s(3), s(9)),
                (s(4), s(10)),
                (s(5), s(11)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn cyclobutadiene() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8)],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(1)),
                (s(1), s(5)),
                (s(2), s(6)),
                (s(3), s(7)),
                (s(4), s(8)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn cyclooctatetraene() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
                s(15),
                s(16),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
                b(15),
                b(16),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(8)),
                (s(8), s(1)),
                (s(1), s(9)),
                (s(2), s(10)),
                (s(3), s(11)),
                (s(4), s(12)),
                (s(5), s(13)),
                (s(6), s(14)),
                (s(7), s(15)),
                (s(8), s(16)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn cycloheptatrienyl_radical() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(1)),
                (s(1), s(8)),
                (s(2), s(9)),
                (s(3), s(10)),
                (s(4), s(11)),
                (s(5), s(12)),
                (s(6), s(13)),
                (s(7), s(14)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn fulvene() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(1)),
                (s(1), s(6)),
                (s(2), s(7)),
                (s(3), s(8)),
                (s(4), s(9)),
                (s(5), s(10)),
                (s(6), s(11)),
                (s(6), s(12)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn heptafulvene() -> Mol {
        Mol {
            sites: vec![
                s(1),
                s(2),
                s(3),
                s(4),
                s(5),
                s(6),
                s(7),
                s(8),
                s(9),
                s(10),
                s(11),
                s(12),
                s(13),
                s(14),
                s(15),
                s(16),
            ],
            elements: vec![
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("C"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
                elem("H"),
            ],
            charges: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            radicals: vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            bonds: vec![
                b(1),
                b(2),
                b(3),
                b(4),
                b(5),
                b(6),
                b(7),
                b(8),
                b(9),
                b(10),
                b(11),
                b(12),
                b(13),
                b(14),
                b(15),
                b(16),
            ],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(1)),
                (s(1), s(8)),
                (s(2), s(9)),
                (s(3), s(10)),
                (s(4), s(11)),
                (s(5), s(12)),
                (s(6), s(13)),
                (s(7), s(14)),
                (s(8), s(15)),
                (s(8), s(16)),
            ],
            orders: vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn aromatic(mol: &Mol) -> Vec<BondId> {
        let mut bonds: Vec<BondId> = perceive(mol).bonds().collect();
        bonds.sort();
        bonds
    }

    fn ring_bonds(count: u32) -> Vec<BondId> {
        (1..=count).map(b).collect()
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
        assert_eq!(aromatic(&benzene()), ring_bonds(6));
    }

    #[test]
    fn pyridine_is_aromatic() {
        assert_eq!(aromatic(&pyridine()), ring_bonds(6));
    }

    #[test]
    fn pyrimidine_is_aromatic() {
        assert_eq!(aromatic(&pyrimidine()), ring_bonds(6));
    }

    #[test]
    fn pyrrole_is_aromatic() {
        assert_eq!(aromatic(&pyrrole()), ring_bonds(5));
    }

    #[test]
    fn furan_is_aromatic() {
        assert_eq!(aromatic(&furan()), ring_bonds(5));
    }

    #[test]
    fn thiophene_is_aromatic() {
        assert_eq!(aromatic(&thiophene()), ring_bonds(5));
    }

    #[test]
    fn imidazole_is_aromatic() {
        assert_eq!(aromatic(&imidazole()), ring_bonds(5));
    }

    #[test]
    fn cyclopentadienyl_anion_is_aromatic() {
        assert_eq!(aromatic(&cyclopentadienyl_anion()), ring_bonds(5));
    }

    #[test]
    fn tropylium_cation_is_aromatic() {
        assert_eq!(aromatic(&tropylium_cation()), ring_bonds(7));
    }

    #[test]
    fn phenyl_radical_is_aromatic() {
        assert_eq!(aromatic(&phenyl_radical()), ring_bonds(6));
    }

    #[test]
    fn tropone_is_aromatic_despite_its_exocyclic_carbonyl() {
        assert_eq!(aromatic(&tropone()), ring_bonds(7));
    }

    #[test]
    fn naphthalene_is_aromatic() {
        assert_eq!(aromatic(&naphthalene()), ring_bonds(11));
    }

    #[test]
    fn indole_is_aromatic() {
        assert_eq!(aromatic(&indole()), ring_bonds(10));
    }

    #[test]
    fn azulene_is_aromatic_around_its_perimeter() {
        assert_eq!(aromatic(&azulene()), ring_bonds(11));
    }

    #[test]
    fn biphenyl_rings_are_aromatic() {
        assert_eq!(aromatic(&biphenyl()), ring_bonds(12));
    }

    #[test]
    fn cyclohexane_is_not_aromatic() {
        assert!(perceive(&cyclohexane()).is_empty());
    }

    #[test]
    fn cyclopentadiene_is_not_aromatic() {
        assert!(perceive(&cyclopentadiene()).is_empty());
    }

    #[test]
    fn cyclobutadiene_is_not_aromatic() {
        assert!(perceive(&cyclobutadiene()).is_empty());
    }

    #[test]
    fn cyclooctatetraene_is_not_aromatic() {
        assert!(perceive(&cyclooctatetraene()).is_empty());
    }

    #[test]
    fn cycloheptatrienyl_radical_is_not_aromatic() {
        assert!(perceive(&cycloheptatrienyl_radical()).is_empty());
    }

    #[test]
    fn fulvene_is_not_aromatic() {
        assert!(perceive(&fulvene()).is_empty());
    }

    #[test]
    fn heptafulvene_is_not_aromatic() {
        assert!(perceive(&heptafulvene()).is_empty());
    }

    #[test]
    fn bond_is_aromatic_only_inside_the_ring() {
        let a = perceive(&benzene());
        assert!(a.bond(b(1)));
        assert!(!a.bond(b(7)));
    }

    #[test]
    fn site_is_aromatic_only_on_the_ring() {
        let a = perceive(&benzene());
        assert!(a.site(s(1)));
        assert!(!a.site(s(7)));
    }

    #[test]
    fn aromatic_sites_are_the_ring_atoms() {
        let mut sites: Vec<SiteId> = perceive(&benzene()).sites().collect();
        sites.sort();
        assert_eq!(sites, (1..=6).map(s).collect::<Vec<_>>());
    }

    #[test]
    fn biphenyl_central_bond_is_not_aromatic() {
        let a = perceive(&biphenyl());
        assert!(!a.bond(b(13)));
        assert!(a.site(s(1)));
        assert!(a.site(s(7)));
    }

    #[test]
    fn unknown_site_is_not_aromatic() {
        assert!(!perceive(&benzene()).site(s(99)));
    }

    #[test]
    fn unknown_bond_is_not_aromatic() {
        assert!(!perceive(&benzene()).bond(b(99)));
    }

    #[test]
    fn bound_view_answers_the_aromaticity_capability() {
        let mol = benzene();
        let perceived = perceive(&mol);
        let view = perceived.bind(&mol);
        assert!(view.is_aromatic(b(1)));
        assert!(!view.is_aromatic(b(7)));
        assert!(view.is_aromatic_site(s(1)));
        assert!(!view.is_aromatic_site(s(7)));
    }

    #[test]
    fn bound_view_forwards_the_skeleton() {
        let mol = pyridine();
        let perceived = perceive(&mol);
        let view = perceived.bind(&mol);

        let mut view_sites: Vec<SiteId> = view.sites().collect();
        let mut mol_sites: Vec<SiteId> = mol.sites().collect();
        view_sites.sort();
        mol_sites.sort();
        assert_eq!(view_sites, mol_sites);
        assert_eq!(view.site_count(), mol.site_count());
        assert!(view.contains_site(s(1)));

        let mut view_bonds: Vec<BondId> = view.bonds().collect();
        let mut mol_bonds: Vec<BondId> = mol.bonds().collect();
        view_bonds.sort();
        mol_bonds.sort();
        assert_eq!(view_bonds, mol_bonds);
        assert_eq!(view.bond_count(), mol.bond_count());
        assert!(view.contains_bond(b(1)));
        assert_eq!(view.bond_endpoints(b(1)), mol.bond_endpoints(b(1)));
        assert_eq!(view.bond_between(s(1), s(2)), mol.bond_between(s(1), s(2)));

        assert_eq!(view.degree(s(1)), mol.degree(s(1)));
        let mut view_incident: Vec<(BondId, SiteId)> = view.bonds_of(s(1)).collect();
        let mut mol_incident: Vec<(BondId, SiteId)> = mol.bonds_of(s(1)).collect();
        view_incident.sort();
        mol_incident.sort();
        assert_eq!(view_incident, mol_incident);
        let mut view_neighbours: Vec<SiteId> = view.neighbors(s(1)).collect();
        let mut mol_neighbours: Vec<SiteId> = mol.neighbors(s(1)).collect();
        view_neighbours.sort();
        mol_neighbours.sort();
        assert_eq!(view_neighbours, mol_neighbours);
    }
}

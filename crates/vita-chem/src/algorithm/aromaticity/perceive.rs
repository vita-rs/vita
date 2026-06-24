use std::collections::{HashMap, HashSet};

use vita_core::{HasElements, HasSites, SiteId};

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
/// Implements [`HasAromaticity`] by forwarding the bond and site skeleton to
/// the molecule and answering aromaticity from the perception, so a computed
/// result reads as the capability that consumers expect. The view costs nothing
/// beyond the two references it holds.
///
/// Obtain via [`Aromaticity::bind`].
pub struct WithAromaticity<'a, M> {
    mol: &'a M,
    aromaticity: &'a Aromaticity,
}

impl<M: HasBonds> HasSites for WithAromaticity<'_, M> {
    fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.mol.sites()
    }

    fn site_count(&self) -> usize {
        self.mol.site_count()
    }

    fn contains_site(&self, site: SiteId) -> bool {
        self.mol.contains_site(site)
    }
}

impl<M: HasBonds> HasBonds for WithAromaticity<'_, M> {
    fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.mol.bonds()
    }

    fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
        self.mol.bond_endpoints(bond)
    }

    fn bond_count(&self) -> usize {
        self.mol.bond_count()
    }

    fn contains_bond(&self, bond: BondId) -> bool {
        self.mol.contains_bond(bond)
    }

    fn bond_between(&self, a: SiteId, b: SiteId) -> Option<BondId> {
        self.mol.bond_between(a, b)
    }

    fn bonds_of(&self, site: SiteId) -> impl Iterator<Item = (BondId, SiteId)> + '_ {
        self.mol.bonds_of(site)
    }

    fn neighbors(&self, site: SiteId) -> impl Iterator<Item = SiteId> + '_ {
        self.mol.neighbors(site)
    }

    fn degree(&self, site: SiteId) -> usize {
        self.mol.degree(site)
    }
}

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
/// - an exocyclic double bond donates none — the p orbital is spent on it (the
///   carbonyl carbon of tropone, the methylene carbon of fulvene);
/// - an otherwise saturated atom with a lone pair donates the pair, two
///   electrons (the nitrogen of pyrrole, the oxygen of furan);
/// - an electron-deficient atom donates none through its empty p orbital (the
///   cation of tropylium, the boron of borole);
/// - a saturated sp³ atom, a radical centre, or one of undefined valence cannot
///   join an aromatic system, and excludes any cycle running through it.
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
/// tropone or the methylene of fulvene).
fn contribution<M>(mol: &M, site: SiteId, rings: &RingMembership) -> Option<u32>
where
    M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons,
{
    let mut endocyclic = 0;
    let mut exocyclic = 0;
    for (bond, _) in mol.bonds_of(site) {
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
    if exocyclic == 1 {
        return Some(0); // the p orbital is spent on the exocyclic π bond
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

use std::collections::{HashMap, HashSet};

use vita_core::{HasElements, SiteId};

use crate::capability::delegation::forward_capabilities;
use crate::utils::{maximum_matching, valence_electrons};
use crate::{BondId, BondOrder, HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// The Kekulé structure resolved for a molecule's aromatic bonds.
///
/// Maps each formerly [`Aromatic`](BondOrder::Aromatic) bond to the localised
/// [`Single`](BondOrder::Single) or [`Double`](BondOrder::Double) order that
/// stands in for it. Bonds that were already localised are absent and unchanged.
///
/// Obtain via [`kekulize`].
pub struct Kekule {
    orders: HashMap<BondId, BondOrder>,
}

impl Kekule {
    /// Returns the localised order resolved for `bond`.
    ///
    /// Returns `None` if `bond` was not aromatic, or is absent from the
    /// molecule.
    pub fn order(&self, bond: BondId) -> Option<BondOrder> {
        self.orders.get(&bond).copied()
    }

    /// Iterates the resolved bonds and the localised order each took.
    pub fn orders(&self) -> impl Iterator<Item = (BondId, BondOrder)> + '_ {
        self.orders.iter().map(|(&bond, &order)| (bond, order))
    }

    /// Returns `true` if the molecule had no aromatic bonds to resolve.
    pub fn is_empty(&self) -> bool {
        self.orders.is_empty()
    }

    /// Binds this structure to `mol`, yielding a view that implements
    /// [`HasBondOrders`].
    ///
    /// The view borrows both, so `mol` stays immutable while it is held — the
    /// structure cannot silently fall out of step with the molecule it
    /// localises. Use it to feed a kekulised molecule to anything that reads the
    /// [`HasBondOrders`] capability.
    pub fn bind<'a, M: HasBondOrders>(&'a self, mol: &'a M) -> WithKekule<'a, M> {
        WithKekule { mol, kekule: self }
    }
}

/// A molecule viewed together with its resolved [`Kekule`] structure.
///
/// Answers bond orders from the resolution — every aromatic bond localised to
/// single or double — and forwards every other core and chem capability to the
/// molecule, so a kekulised result reads as the [`HasBondOrders`] capability its
/// consumers expect, at no cost beyond the two references it holds.
///
/// Obtain via [`Kekule::bind`].
pub struct WithKekule<'a, M> {
    mol: &'a M,
    kekule: &'a Kekule,
}

forward_capabilities!(
    WithKekule,
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
    HasAromaticity,
    HasBonds,
    HasFormalCharges,
    HasHybridizations,
    HasPartialCharges,
    HasRadicalElectrons,
);

impl<M: HasBondOrders> HasBondOrders for WithKekule<'_, M> {
    fn bond_order(&self, bond: BondId) -> BondOrder {
        self.kekule
            .order(bond)
            .unwrap_or_else(|| self.mol.bond_order(bond))
    }
}

/// Resolves a molecule's aromatic bonds into a Kekulé structure.
///
/// Each [`Aromatic`](BondOrder::Aromatic) bond is localised to a
/// [`Single`](BondOrder::Single) or [`Double`](BondOrder::Double) order so that
/// every atom's valence is met. The choice is a matching: an atom that brings
/// one electron to the π system — the carbon of benzene, the nitrogen of
/// pyridine — must take exactly one double bond, while a lone-pair donor (the
/// nitrogen of pyrrole) or an electron-deficient centre (the cation of
/// tropylium) takes none. Pairing each demanding atom with a neighbour across a
/// shared double bond is a perfect matching of the subgraph they span, blossoms
/// and all — enough to localise even the odd rings azulene and other
/// non-alternant systems carry.
///
/// Returns `None` when no Kekulé structure exists: an atom that cannot be sp²
/// aromatic, or a π system left with an odd, unpairable electron. A molecule
/// with no aromatic bonds is already localised and resolves to an empty
/// [`Kekule`].
///
/// Kekulisation is the inverse of [`perceive`](super::perceive())'s
/// representation: perceive reads the aromatic system a localised structure
/// stands for; kekulize writes a localised structure an aromatic system stands
/// for.
///
/// # Complexity
///
/// O(V³), dominated by the matching.
pub fn kekulize<M>(mol: &M) -> Option<Kekule>
where
    M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons,
{
    // Without aromatic bonds the molecule is already localised.
    let aromatic: Vec<BondId> = mol
        .bonds()
        .filter(|&bond| mol.bond_order(bond) == BondOrder::Aromatic)
        .collect();
    if aromatic.is_empty() {
        return Some(Kekule {
            orders: HashMap::new(),
        });
    }

    // The atoms that must each take one double bond are the vertices to match;
    // a `None` demand is an atom that cannot be sp² aromatic, so no Kekulé
    // structure exists.
    let mut vertices: Vec<SiteId> = Vec::new();
    let mut seen: HashSet<SiteId> = HashSet::new();
    for &bond in &aromatic {
        let (a, b) = mol.bond_endpoints(bond);
        for site in [a, b] {
            if seen.insert(site) && needs_double(mol, site)? {
                vertices.push(site);
            }
        }
    }
    vertices.sort();

    // Aromatic bonds between two demanding atoms are the matchable edges.
    let index: HashMap<SiteId, usize> = vertices.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertices.len()];
    for &bond in &aromatic {
        let (a, b) = mol.bond_endpoints(bond);
        if let (Some(&i), Some(&j)) = (index.get(&a), index.get(&b)) {
            adjacency[i].push(j);
            adjacency[j].push(i);
        }
    }

    // Every demanding atom must be paired; an unmatched one is an unpaired π
    // electron — a radical, not a closed-shell Kekulé structure.
    let matching = maximum_matching(&adjacency);
    if matching.iter().any(Option::is_none) {
        return None;
    }

    // The matched aromatic bonds become double, every other aromatic bond single.
    let mut orders: HashMap<BondId, BondOrder> = HashMap::new();
    for &bond in &aromatic {
        let (a, b) = mol.bond_endpoints(bond);
        let double = match (index.get(&a), index.get(&b)) {
            (Some(&i), Some(&j)) => matching[i] == Some(j),
            _ => false,
        };
        orders.insert(
            bond,
            if double {
                BondOrder::Double
            } else {
                BondOrder::Single
            },
        );
    }

    Some(Kekule { orders })
}

/// Whether an aromatic atom must take exactly one double bond.
///
/// With its aromatic bonds counted as single, the electrons left over the σ
/// skeleton — `valence electrons − formal charge − bond-order sum − radicals` —
/// fall into the perpendicular p orbital. An odd number leaves one unpaired
/// electron that a single double bond pairs off (the carbon of benzene, the
/// nitrogen of pyridine); an even number is a lone pair donated whole (the
/// nitrogen of pyrrole) or an empty orbital (the cation of tropylium), wanting
/// no double bond.
///
/// Returns `None` when the skeleton is overfull: the atom cannot be sp²
/// aromatic, so no Kekulé structure exists.
fn needs_double<M>(mol: &M, site: SiteId) -> Option<bool>
where
    M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons,
{
    let mut bonding = 0i32;
    for (bond, _) in mol.bonds_of(site) {
        bonding += match mol.bond_order(bond) {
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => return None,
        };
    }
    let electrons = valence_electrons(mol.element(site))? as i32;
    let free =
        electrons - mol.formal_charge(site) as i32 - bonding - mol.radical_electron(site) as i32;
    if free < 0 {
        return None;
    }
    Some(free % 2 == 1)
}

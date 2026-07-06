use vita_core::{HasElements, SiteId};

use crate::algorithm::utils::{FxHashSet, SortedMap, maximum_matching, valence_electrons};
use crate::capability::delegation::forward_capabilities;
use crate::{BondId, BondOrder, HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// The Kekulé structure resolved for a molecule's aromatic bonds.
///
/// Maps each formerly [`Aromatic`](BondOrder::Aromatic) bond to the localised
/// [`Single`](BondOrder::Single) or [`Double`](BondOrder::Double) order that
/// stands in for it. Bonds that were already localised are absent and unchanged.
///
/// Obtain via [`kekulize`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Kekule {
    orders: SortedMap<BondId, BondOrder>,
}

impl Kekule {
    /// Number of bonds the Kekulé structure resolves.
    pub fn len(&self) -> usize {
        self.orders.len()
    }

    /// Returns `true` if the molecule had no aromatic bonds to resolve.
    pub fn is_empty(&self) -> bool {
        self.orders.is_empty()
    }

    /// Returns the localised order resolved for `bond`.
    ///
    /// Returns `None` if `bond` was not aromatic, or is absent from the
    /// molecule.
    pub fn order(&self, bond: BondId) -> Option<BondOrder> {
        self.orders.get(&bond).copied()
    }

    /// Iterates the resolved bonds and the localised order each took, in
    /// ascending bond order.
    pub fn orders(&self) -> impl Iterator<Item = (BondId, BondOrder)> + '_ {
        self.orders.iter().map(|(&bond, &order)| (bond, order))
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

impl<M> Copy for WithKekule<'_, M> {}

impl<M> Clone for WithKekule<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
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
/// O(V³) time and O(V + E) space, over the molecule's `V` sites and `E` bonds,
/// dominated by the matching.
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
            orders: SortedMap::from_pairs(Vec::new()),
        });
    }

    // The atoms that must each take one double bond are the vertices to match;
    // a `None` demand is an atom that cannot be sp² aromatic, so no Kekulé
    // structure exists. Sorted, the vertices are their own index into the
    // matching by binary search.
    let mut vertices: Vec<SiteId> = Vec::new();
    let mut seen: FxHashSet<SiteId> = FxHashSet::default();
    for &bond in &aromatic {
        let (a, b) = mol.bond_endpoints(bond);
        for site in [a, b] {
            if seen.insert(site) && needs_double(mol, site)? {
                vertices.push(site);
            }
        }
    }
    vertices.sort_unstable();

    // Aromatic bonds between two demanding atoms are the matchable edges.
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertices.len()];
    for &bond in &aromatic {
        let (a, b) = mol.bond_endpoints(bond);
        if let (Ok(i), Ok(j)) = (vertices.binary_search(&a), vertices.binary_search(&b)) {
            adjacency[i].push(j);
            adjacency[j].push(i);
        }
    }

    // Every demanding atom must be paired; an unmatched one is an unpaired π
    // electron — a radical, not a closed-shell Kekulé structure.
    let matching = maximum_matching(&adjacency);
    if !matching.is_perfect() {
        return None;
    }

    // The matched aromatic bonds become double, every other aromatic bond single.
    let orders = aromatic.iter().map(|&bond| {
        let (a, b) = mol.bond_endpoints(bond);
        let double = match (vertices.binary_search(&a), vertices.binary_search(&b)) {
            (Ok(i), Ok(j)) => matching.mate(i) == Some(j),
            _ => false,
        };
        let order = if double {
            BondOrder::Double
        } else {
            BondOrder::Single
        };
        (bond, order)
    });

    Some(Kekule {
        orders: SortedMap::from_pairs(orders),
    })
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

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasSites};

    use crate::BondOrder::{Aromatic, Double, Single};
    use crate::{BondId, HasBonds};

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

    fn doubles(kekule: &Kekule) -> usize {
        kekule
            .orders()
            .filter(|&(_, order)| order == Double)
            .count()
    }

    fn incident_doubles(mol: &Mol, kekule: &Kekule, site: SiteId) -> usize {
        mol.bonds_of(site)
            .filter(|&(bond, _)| kekule.order(bond) == Some(Double))
            .count()
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 6, Aromatic),
                (6, 6, 1, Aromatic),
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 1, Aromatic),
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 1, Aromatic),
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 6, Aromatic),
                (6, 6, 7, Aromatic),
                (7, 7, 1, Aromatic),
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 1, Aromatic),
                (6, 1, 6, Single),
                (7, 2, 7, Single),
                (8, 3, 8, Single),
                (9, 4, 9, Single),
                (10, 5, 10, Single),
            ],
        )
    }

    fn cyclopentadienyl_ring() -> Mol {
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
            ],
            &[
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 1, Aromatic),
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 6, Aromatic),
                (6, 6, 7, Aromatic),
                (7, 7, 1, Aromatic),
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 9, Aromatic),
                (5, 9, 10, Aromatic),
                (6, 10, 1, Aromatic),
                (7, 5, 6, Aromatic),
                (8, 6, 7, Aromatic),
                (9, 7, 8, Aromatic),
                (10, 8, 9, Aromatic),
                (11, 10, 5, Aromatic),
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 1, Aromatic),
                (6, 2, 6, Aromatic),
                (7, 6, 7, Aromatic),
                (8, 7, 8, Aromatic),
                (9, 8, 9, Aromatic),
                (10, 9, 10, Aromatic),
                (11, 10, 1, Aromatic),
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
                (1, 1, 2, Aromatic),
                (2, 2, 3, Aromatic),
                (3, 3, 4, Aromatic),
                (4, 4, 5, Aromatic),
                (5, 5, 6, Aromatic),
                (6, 6, 1, Aromatic),
                (7, 7, 8, Aromatic),
                (8, 8, 9, Aromatic),
                (9, 9, 10, Aromatic),
                (10, 10, 11, Aromatic),
                (11, 11, 12, Aromatic),
                (12, 12, 7, Aromatic),
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
    fn already_localised_molecule_resolves_to_empty() {
        assert!(kekulize(&ethane()).unwrap().is_empty());
    }

    #[test]
    fn benzene_takes_three_alternating_double_bonds() {
        let mol = benzene();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 3);
        assert!((1..=6).all(|i| incident_doubles(&mol, &kekule, s(i)) == 1));
    }

    #[test]
    fn pyridine_takes_three_double_bonds() {
        assert_eq!(doubles(&kekulize(&pyridine()).unwrap()), 3);
    }

    #[test]
    fn pyrrole_leaves_its_nitrogen_single() {
        let mol = pyrrole();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 2);
        assert_eq!(incident_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn tropylium_cation_leaves_its_cation_single() {
        let mol = tropylium_cation();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 3);
        assert_eq!(incident_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn cyclopentadienyl_anion_leaves_its_carbanion_single() {
        let mol = cyclopentadienyl_anion();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 2);
        assert_eq!(incident_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn tropone_leaves_its_carbonyl_carbon_single() {
        let mol = tropone();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 3);
        assert_eq!(incident_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn closed_shell_odd_ring_has_no_kekule_structure() {
        assert!(kekulize(&cyclopentadienyl_ring()).is_none());
    }

    #[test]
    fn order_is_none_for_a_non_aromatic_bond() {
        let kekule = kekulize(&benzene()).unwrap();
        assert!(kekule.order(b(1)).is_some());
        assert_eq!(kekule.order(b(7)), None);
        assert_eq!(kekule.order(b(99)), None);
    }

    #[test]
    fn azulene_kekulises_its_non_alternant_rings() {
        assert_eq!(doubles(&kekulize(&azulene()).unwrap()), 5);
    }

    #[test]
    fn naphthalene_takes_five_double_bonds() {
        assert_eq!(doubles(&kekulize(&naphthalene()).unwrap()), 5);
    }

    #[test]
    fn biphenyl_kekulises_each_ring_leaving_the_link() {
        let kekule = kekulize(&biphenyl()).unwrap();
        assert_eq!(doubles(&kekule), 6);
        assert_eq!(kekule.order(b(13)), None);
    }

    #[test]
    fn imidazole_leaves_its_pyrrole_nitrogen_single() {
        let mol = imidazole();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 2);
        assert_eq!(incident_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn every_aromatic_bond_is_localised() {
        let kekule = kekulize(&benzene()).unwrap();
        assert_eq!(kekule.orders().count(), 6);
        assert!(kekule.orders().all(|(_, order)| order != Aromatic));
    }

    #[test]
    fn len_counts_the_resolved_bonds() {
        assert_eq!(kekulize(&benzene()).unwrap().len(), 6);
        assert_eq!(kekulize(&ethane()).unwrap().len(), 0);
    }

    #[test]
    fn orders_are_listed_in_ascending_bond_order() {
        let kekule = kekulize(&benzene()).unwrap();
        let bonds: Vec<BondId> = kekule.orders().map(|(bond, _)| bond).collect();
        assert!(bonds.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn resolution_is_independent_of_input_order() {
        let resolved = |m: &Mol| -> Vec<BondId> {
            kekulize(m)
                .unwrap()
                .orders()
                .map(|(bond, _)| bond)
                .collect()
        };
        assert_eq!(resolved(&benzene()), resolved(&reversed(&benzene())));
    }

    #[test]
    fn bound_view_answers_the_bond_order_capability() {
        let mol = benzene();
        let kekule = kekulize(&mol).unwrap();
        let view = kekule.bind(&mol);
        assert!((1..=6).all(|i| view.bond_order(b(i)) != Aromatic));
        assert_eq!(view.bond_order(b(7)), Single);
    }

    #[test]
    fn bound_view_forwards_the_skeleton() {
        let mol = pyridine();
        let kekule = kekulize(&mol).unwrap();
        let view = kekule.bind(&mol);
        assert_eq!(view.bond_count(), mol.bond_count());
        assert_eq!(view.bond_endpoints(b(1)), mol.bond_endpoints(b(1)));
    }
}

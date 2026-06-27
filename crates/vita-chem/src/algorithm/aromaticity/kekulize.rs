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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondId, HasBonds};
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
            orders: vec![BondOrder::Single; 7],
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
            charges: vec![0; 12],
            radicals: vec![0; 12],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 11],
            radicals: vec![0; 11],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 10],
            radicals: vec![0; 10],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 9],
            radicals: vec![0; 9],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 9],
            radicals: vec![0; 9],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            radicals: vec![0; 10],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            radicals: vec![0; 14],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 14],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 14],
            radicals: vec![0; 14],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 18],
            radicals: vec![0; 18],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 16],
            radicals: vec![0; 16],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 18],
            radicals: vec![0; 18],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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
            charges: vec![0; 22],
            radicals: vec![0; 22],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
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

    fn cyclopentadienyl_ring() -> Mol {
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
            charges: vec![0; 10],
            radicals: vec![0; 10],
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
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        }
    }

    fn doubles(kekule: &Kekule) -> usize {
        kekule
            .orders()
            .filter(|(_, order)| *order == BondOrder::Double)
            .count()
    }

    fn incident_ring_doubles(mol: &Mol, kekule: &Kekule, site: SiteId) -> usize {
        mol.bonds_of(site)
            .filter(|&(bond, _)| kekule.order(bond) == Some(BondOrder::Double))
            .count()
    }

    #[test]
    fn already_localised_molecule_resolves_to_empty() {
        let kekule = kekulize(&ethane()).unwrap();
        assert!(kekule.is_empty());
    }

    #[test]
    fn benzene_alternates_three_double_bonds() {
        let mol = benzene();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 3);
        assert!((1..=6).all(|i| incident_ring_doubles(&mol, &kekule, s(i)) == 1));
    }

    #[test]
    fn pyridine_alternates_three_double_bonds() {
        assert_eq!(doubles(&kekulize(&pyridine()).unwrap()), 3);
    }

    #[test]
    fn pyrrole_leaves_its_nitrogen_single() {
        let mol = pyrrole();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 2);
        assert_eq!(incident_ring_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn furan_leaves_its_oxygen_single() {
        let mol = furan();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 2);
        assert_eq!(incident_ring_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn imidazole_leaves_its_pyrrole_nitrogen_single() {
        let mol = imidazole();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 2);
        assert_eq!(incident_ring_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn cyclopentadienyl_anion_takes_two_double_bonds() {
        let mol = cyclopentadienyl_anion();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 2);
        assert_eq!(incident_ring_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn tropylium_cation_takes_three_double_bonds() {
        let mol = tropylium_cation();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 3);
        assert_eq!(incident_ring_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn cycloheptatrienyl_radical_leaves_its_radical_carbon_single() {
        let mol = cycloheptatrienyl_radical();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 3);
        assert_eq!(incident_ring_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn tropone_leaves_its_carbonyl_carbon_single() {
        let mol = tropone();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 3);
        assert_eq!(incident_ring_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn naphthalene_takes_five_double_bonds() {
        assert_eq!(doubles(&kekulize(&naphthalene()).unwrap()), 5);
    }

    #[test]
    fn indole_leaves_its_nitrogen_single() {
        let mol = indole();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 4);
        assert_eq!(incident_ring_doubles(&mol, &kekule, s(1)), 0);
    }

    #[test]
    fn azulene_takes_five_double_bonds() {
        let mol = azulene();
        let kekule = kekulize(&mol).unwrap();
        assert_eq!(doubles(&kekule), 5);
        assert!((1..=10).all(|i| incident_ring_doubles(&mol, &kekule, s(i)) == 1));
    }

    #[test]
    fn biphenyl_kekulises_each_ring_leaving_the_central_bond() {
        let kekule = kekulize(&biphenyl()).unwrap();
        assert_eq!(doubles(&kekule), 6);
        assert_eq!(kekule.order(b(13)), None);
    }

    #[test]
    fn closed_shell_odd_ring_has_no_kekule_structure() {
        assert!(kekulize(&cyclopentadienyl_ring()).is_none());
    }

    #[test]
    fn resolved_bonds_localise_every_aromatic_bond() {
        let kekule = kekulize(&benzene()).unwrap();
        let resolved: Vec<BondId> = kekule.orders().map(|(bond, _)| bond).collect();
        assert_eq!(resolved.len(), 6);
        assert!(
            kekule
                .orders()
                .all(|(_, order)| order != BondOrder::Aromatic)
        );
    }

    #[test]
    fn order_is_none_for_a_non_aromatic_bond() {
        let kekule = kekulize(&benzene()).unwrap();
        assert!(kekule.order(b(1)).is_some());
        assert_eq!(kekule.order(b(7)), None);
        assert_eq!(kekule.order(b(99)), None);
    }

    #[test]
    fn bound_view_answers_the_bond_order_capability() {
        let mol = benzene();
        let kekule = kekulize(&mol).unwrap();
        let view = kekule.bind(&mol);

        for bond in (1..=6).map(b) {
            assert_ne!(view.bond_order(bond), BondOrder::Aromatic);
        }
        let doubled = (1..=6)
            .map(b)
            .filter(|&bond| view.bond_order(bond) == BondOrder::Double)
            .count();
        assert_eq!(doubled, 3);
        assert_eq!(view.bond_order(b(7)), BondOrder::Single);
    }

    #[test]
    fn bound_view_forwards_the_skeleton() {
        let mol = pyridine();
        let kekule = kekulize(&mol).unwrap();
        let view = kekule.bind(&mol);

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

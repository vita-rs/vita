use vita_core::{HasElements, SiteId};

use super::lone_pairs::lone_pairs;
use crate::capability::delegation::forward_capabilities;
use crate::{
    HasBondOrders, HasFormalCharges, HasHybridizations, HasRadicalElectrons, Hybridization,
};

/// Hybridization of `site` from its electron-domain count.
///
/// Counts electron domains — bonded neighbours plus [`lone_pairs`] — and names
/// the geometry: one or none is [`S`](Hybridization::S), two
/// [`Sp`](Hybridization::Sp), three [`Sp2`](Hybridization::Sp2), four
/// [`Sp3`](Hybridization::Sp3), five [`Sp3d`](Hybridization::Sp3d), six
/// [`Sp3d2`](Hybridization::Sp3d2), seven [`Sp3d3`](Hybridization::Sp3d3), and
/// eight or more [`Other`](Hybridization::Other).
///
/// This is the steric (VSEPR) count, blind to conjugation and to coordinates:
/// an amide nitrogen is [`Sp3`](Hybridization::Sp3) rather than `Sp2`, and
/// square-planar [`Sp2d`](Hybridization::Sp2d) is not told apart from
/// [`Sp3`](Hybridization::Sp3).
///
/// Returns `None` exactly when [`lone_pairs`] does — a d-/f-block element, an
/// aromatic (delocalised) bond, or an impossible valence — leaving the domain
/// count undefined.
///
/// # Complexity
///
/// O(degree) time.
pub fn hybridization<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>(
    mol: &M,
    site: SiteId,
) -> Option<Hybridization> {
    let domains = mol.degree(site) as u32 + lone_pairs(mol, site)?;
    Some(match domains {
        0 | 1 => Hybridization::S,
        2 => Hybridization::Sp,
        3 => Hybridization::Sp2,
        4 => Hybridization::Sp3,
        5 => Hybridization::Sp3d,
        6 => Hybridization::Sp3d2,
        7 => Hybridization::Sp3d3,
        _ => Hybridization::Other,
    })
}

/// A molecule viewed together with its computed hybridizations.
///
/// Answers hybridization from [`hybridization`] — the catch-all
/// [`Other`](Hybridization::Other) where it is undefined — and forwards every
/// other core and chem capability to the molecule, so a computed result reads
/// as the [`HasHybridizations`] capability its consumers expect, at no cost
/// beyond the single reference it holds.
///
/// Obtain via [`new`](Self::new).
pub struct WithHybridizations<'a, M> {
    mol: &'a M,
}

impl<'a, M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>
    WithHybridizations<'a, M>
{
    /// Views `mol` together with its computed hybridizations, yielding a value
    /// that implements [`HasHybridizations`].
    ///
    /// The view borrows `mol` and computes each hybridization on demand, holding
    /// no buffer; it cannot fall out of step with the molecule. Feed it to
    /// anything that reads the [`HasHybridizations`] capability.
    pub fn new(mol: &'a M) -> Self {
        Self { mol }
    }
}

forward_capabilities!(
    WithHybridizations,
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
    HasBondOrders,
    HasBonds,
    HasFormalCharges,
    HasPartialCharges,
    HasRadicalElectrons,
);

impl<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons> HasHybridizations
    for WithHybridizations<'_, M>
{
    fn hybridization(&self, site: SiteId) -> Hybridization {
        hybridization(self.mol, site).unwrap_or(Hybridization::Other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondId, BondOrder, HasBonds};
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

    fn atom(symbol: &str, charge: i8, radicals: u8, bonds: &[BondOrder]) -> Mol {
        let mut mol = Mol {
            sites: vec![s(1)],
            elements: vec![elem(symbol)],
            charges: vec![charge],
            radicals: vec![radicals],
            bonds: Vec::new(),
            endpoints: Vec::new(),
            orders: Vec::new(),
        };
        for (i, &order) in bonds.iter().enumerate() {
            let neighbour = s(i as u32 + 2);
            mol.sites.push(neighbour);
            mol.elements.push(elem("H"));
            mol.charges.push(0);
            mol.radicals.push(0);
            mol.bonds.push(b(i as u32 + 1));
            mol.endpoints.push((s(1), neighbour));
            mol.orders.push(order);
        }
        mol
    }

    #[test]
    fn carbon_geometries() {
        assert_eq!(
            hybridization(&atom("C", 0, 0, &[BondOrder::Double; 2]), s(1)),
            Some(Hybridization::Sp),
        );
        assert_eq!(
            hybridization(&atom("C", 1, 0, &[BondOrder::Single; 3]), s(1)),
            Some(Hybridization::Sp2),
        );
        assert_eq!(
            hybridization(&atom("C", 0, 0, &[BondOrder::Single; 4]), s(1)),
            Some(Hybridization::Sp3),
        );
    }

    #[test]
    fn lone_pairs_raise_the_count() {
        assert_eq!(
            hybridization(&atom("O", 0, 0, &[BondOrder::Single; 2]), s(1)),
            Some(Hybridization::Sp3),
        );
        assert_eq!(
            hybridization(&atom("N", 0, 0, &[BondOrder::Single; 3]), s(1)),
            Some(Hybridization::Sp3),
        );
        assert_eq!(
            hybridization(&atom("O", 0, 0, &[BondOrder::Double]), s(1)),
            Some(Hybridization::Sp2),
        );
        assert_eq!(
            hybridization(&atom("N", 0, 0, &[BondOrder::Triple]), s(1)),
            Some(Hybridization::Sp),
        );
    }

    #[test]
    fn hypervalent_geometries() {
        assert_eq!(
            hybridization(&atom("P", 0, 0, &[BondOrder::Single; 5]), s(1)),
            Some(Hybridization::Sp3d),
        );
        assert_eq!(
            hybridization(&atom("S", 0, 0, &[BondOrder::Single; 6]), s(1)),
            Some(Hybridization::Sp3d2),
        );
        assert_eq!(
            hybridization(&atom("I", 0, 0, &[BondOrder::Single; 7]), s(1)),
            Some(Hybridization::Sp3d3),
        );
    }

    #[test]
    fn zero_or_one_domain_is_s() {
        assert_eq!(
            hybridization(&atom("Na", 1, 0, &[]), s(1)),
            Some(Hybridization::S),
        );
        assert_eq!(
            hybridization(&atom("H", 0, 0, &[BondOrder::Single]), s(1)),
            Some(Hybridization::S),
        );
    }

    #[test]
    fn eight_domains_is_other() {
        assert_eq!(
            hybridization(&atom("Xe", 0, 0, &[BondOrder::Single; 8]), s(1)),
            Some(Hybridization::Other),
        );
    }

    #[test]
    fn undefined_when_lone_pairs_undefined() {
        assert_eq!(
            hybridization(&atom("C", 0, 0, &[BondOrder::Aromatic; 2]), s(1)),
            None,
        );
        assert_eq!(
            hybridization(&atom("Fe", 0, 0, &[BondOrder::Single; 2]), s(1)),
            None
        );
    }

    #[test]
    fn bound_view_answers_the_hybridization_capability() {
        let mol = atom("C", 0, 0, &[BondOrder::Single; 4]);
        let view = WithHybridizations::new(&mol);
        assert_eq!(view.hybridization(s(1)), Hybridization::Sp3);
    }

    #[test]
    fn bound_view_takes_other_where_undefined() {
        let mol = atom("Fe", 0, 0, &[BondOrder::Single; 2]);
        let view = WithHybridizations::new(&mol);
        assert_eq!(view.hybridization(s(1)), Hybridization::Other);
    }

    #[test]
    fn bound_view_forwards_the_skeleton() {
        let mol = atom("C", 0, 0, &[BondOrder::Single; 2]);
        let view = WithHybridizations::new(&mol);

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

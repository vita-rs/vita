use vita_core::{HasElements, SiteId};

use super::lone_pairs;
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
/// O(d) time and O(1) space, where `d` is the degree of `site`, assuming
/// [`degree`](crate::HasBonds::degree) and [`bonds_of`](crate::HasBonds::bonds_of)
/// run in O(degree).
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

    use vita_core::{Element, HasSites};

    use crate::{BondId, BondOrder, HasBonds};

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
        formal_charges: Vec<i8>,
        radicals: Vec<u8>,
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

    impl HasFormalCharges for Mol {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.formal_charges[i]
        }
    }

    impl HasRadicalElectrons for Mol {
        fn radical_electron(&self, site: SiteId) -> u8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.radicals[i]
        }
    }

    fn atom(symbol: &str, charge: i8, radicals: u8, orders: &[BondOrder]) -> Mol {
        let n = orders.len() as u32;
        let mut sites = vec![s(1)];
        let mut elements = vec![elem(symbol)];
        let mut formal_charges = vec![charge];
        let mut radical_counts = vec![radicals];
        for i in 2..=n + 1 {
            sites.push(s(i));
            elements.push(elem("H"));
            formal_charges.push(0);
            radical_counts.push(0);
        }
        Mol {
            sites,
            elements,
            bonds: (1..=n).map(b).collect(),
            endpoints: (2..=n + 1).map(|i| (s(1), s(i))).collect(),
            orders: orders.to_vec(),
            formal_charges,
            radicals: radical_counts,
        }
    }

    #[test]
    fn no_or_one_domain_is_s() {
        assert_eq!(
            hybridization(&atom("H", 0, 0, &[]), s(1)),
            Some(Hybridization::S)
        );
        let bonded = atom("H", 0, 0, &[BondOrder::Single]);
        assert_eq!(hybridization(&bonded, s(1)), Some(Hybridization::S));
    }

    #[test]
    fn two_domains_are_sp() {
        let carbon_dioxide = atom("C", 0, 0, &[BondOrder::Double, BondOrder::Double]);
        assert_eq!(
            hybridization(&carbon_dioxide, s(1)),
            Some(Hybridization::Sp)
        );
    }

    #[test]
    fn three_domains_are_sp2() {
        let borane = atom("B", 0, 0, &[BondOrder::Single; 3]);
        assert_eq!(hybridization(&borane, s(1)), Some(Hybridization::Sp2));
    }

    #[test]
    fn four_domains_are_sp3() {
        let methane = atom("C", 0, 0, &[BondOrder::Single; 4]);
        assert_eq!(hybridization(&methane, s(1)), Some(Hybridization::Sp3));
    }

    #[test]
    fn five_domains_are_sp3d() {
        let phosphorus_pentachloride = atom("P", 0, 0, &[BondOrder::Single; 5]);
        assert_eq!(
            hybridization(&phosphorus_pentachloride, s(1)),
            Some(Hybridization::Sp3d),
        );
    }

    #[test]
    fn six_domains_are_sp3d2() {
        let sulfur_hexafluoride = atom("S", 0, 0, &[BondOrder::Single; 6]);
        assert_eq!(
            hybridization(&sulfur_hexafluoride, s(1)),
            Some(Hybridization::Sp3d2),
        );
    }

    #[test]
    fn seven_domains_are_sp3d3() {
        let iodine_heptafluoride = atom("I", 0, 0, &[BondOrder::Single; 7]);
        assert_eq!(
            hybridization(&iodine_heptafluoride, s(1)),
            Some(Hybridization::Sp3d3),
        );
    }

    #[test]
    fn eight_or_more_domains_are_other() {
        let over_coordinated = atom("Xe", 0, 0, &[BondOrder::Single; 8]);
        assert_eq!(
            hybridization(&over_coordinated, s(1)),
            Some(Hybridization::Other),
        );
    }

    #[test]
    fn undefined_when_the_domain_count_is_undefined() {
        assert_eq!(hybridization(&atom("Fe", 0, 0, &[]), s(1)), None);
        let aromatic = atom("C", 0, 0, &[BondOrder::Aromatic]);
        assert_eq!(hybridization(&aromatic, s(1)), None);
    }

    #[test]
    fn lone_pairs_count_as_domains() {
        let water = atom("O", 0, 0, &[BondOrder::Single, BondOrder::Single]);
        assert_eq!(hybridization(&water, s(1)), Some(Hybridization::Sp3));
    }

    #[test]
    fn view_reports_the_computed_hybridization() {
        let methane = atom("C", 0, 0, &[BondOrder::Single; 4]);
        let view = WithHybridizations::new(&methane);
        assert_eq!(view.hybridization(s(1)), Hybridization::Sp3);
    }

    #[test]
    fn view_falls_back_to_other_where_undefined() {
        let iron = atom("Fe", 0, 0, &[]);
        let view = WithHybridizations::new(&iron);
        assert_eq!(view.hybridization(s(1)), Hybridization::Other);
    }

    #[test]
    fn view_forwards_other_capabilities() {
        let mol = atom("C", 0, 0, &[BondOrder::Double]);
        let view = WithHybridizations::new(&mol);
        assert_eq!(view.element(s(1)), elem("C"));
        assert_eq!(view.bond_order(b(1)), BondOrder::Double);
    }
}

use vita_core::{HasElements, HasSites, SiteId};

use crate::algorithm::conjugation::systems;
use crate::algorithm::utils::{FxHashSet, SortedMap};
use crate::algorithm::valence::lone_pairs;
use crate::capability::delegation::forward_capabilities;
use crate::{
    BondOrder, HasBondOrders, HasFormalCharges, HasHybridizations, HasRadicalElectrons,
    Hybridization,
};

/// The hybridization of every site, perceived from the molecular graph.
///
/// Sites whose count the model cannot fix — d- and f-block elements,
/// impossible structures — carry no label.
///
/// Obtain via [`perceive`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hybridizations {
    labels: SortedMap<SiteId, Hybridization>,
}

impl Hybridizations {
    /// Number of labeled sites.
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// Returns `true` if no site carries a label.
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// Returns the hybridization of `site`.
    ///
    /// Returns `None` if `site` is absent from the molecule or carries no
    /// label.
    pub fn hybridization(&self, site: SiteId) -> Option<Hybridization> {
        self.labels.get(&site).copied()
    }

    /// Iterates the labeled sites with their hybridizations, in ascending
    /// site order.
    pub fn hybridizations(&self) -> impl Iterator<Item = (SiteId, Hybridization)> + '_ {
        self.labels.iter().map(|(&site, &label)| (site, label))
    }

    /// Binds this perception to `mol`, yielding a view that implements
    /// [`HasHybridizations`].
    ///
    /// The view borrows both, so `mol` stays immutable while it is held — the
    /// perception cannot silently fall out of step with the molecule it
    /// describes. Use it to feed a perceived molecule to anything that reads
    /// the [`HasHybridizations`] capability.
    pub fn bind<'a, M: HasSites>(&'a self, mol: &'a M) -> WithHybridizations<'a, M> {
        WithHybridizations {
            mol,
            hybridizations: self,
        }
    }
}

/// A molecule viewed together with its perceived [`Hybridizations`].
///
/// Answers hybridization from the perception — the catch-all
/// [`Other`](Hybridization::Other) for unlabeled sites — and forwards every
/// other core and chem capability to the molecule, so a computed result reads
/// as the [`HasHybridizations`] capability its consumers expect, at no cost
/// beyond the two references it holds.
///
/// Obtain via [`Hybridizations::bind`].
pub struct WithHybridizations<'a, M> {
    mol: &'a M,
    hybridizations: &'a Hybridizations,
}

impl<M> Copy for WithHybridizations<'_, M> {}

impl<M> Clone for WithHybridizations<'_, M> {
    fn clone(&self) -> Self {
        *self
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
    HasStereoConfigurations,
);

impl<M: HasSites> HasHybridizations for WithHybridizations<'_, M> {
    fn hybridization(&self, site: SiteId) -> Hybridization {
        self.hybridizations
            .hybridization(site)
            .unwrap_or(Hybridization::Other)
    }
}

/// Perceives the hybridization of every site from the molecular graph.
///
/// A site hybridizes the orbitals its electrons keep in the σ-frame: the
/// electron domains — bonded neighbors plus lone pairs — minus the lone
/// pairs it donates into a conjugated system, which relocate to pure
/// p-orbitals. An amide nitrogen counts four domains, donates one pair, and
/// comes out [`Sp2`](Hybridization::Sp2); an azide is [`Sp`](Hybridization::Sp)
/// end to end. One domain or none names [`S`](Hybridization::S), two
/// [`Sp`](Hybridization::Sp), three [`Sp2`](Hybridization::Sp2), four
/// [`Sp3`](Hybridization::Sp3), five [`Sp3d`](Hybridization::Sp3d), six
/// [`Sp3d2`](Hybridization::Sp3d2), seven [`Sp3d3`](Hybridization::Sp3d3),
/// and more [`Other`](Hybridization::Other). Sites on
/// [`Aromatic`](BondOrder::Aromatic) bonds are [`Sp2`](Hybridization::Sp2) by
/// declaration.
///
/// The perception is topological and inherits the conjugation model (see
/// [`systems`]): it reports the idealised, maximally planar answer — square
/// planar [`Sp2d`](Hybridization::Sp2d) is only reachable from geometry — and
/// it is independent of the drawn resonance form, the donation exactly
/// offsetting the shifted bonds and charges. A site is left unlabeled when
/// no exact count exists: a d- or f-block element, or arithmetic describing
/// an impossible structure.
///
/// # Complexity
///
/// O((V + E) · log (V + E)) time and O(V + E) space, over the molecule's `V`
/// sites and `E` bonds, assuming [`degree`](crate::HasBonds::degree) and
/// [`bonds_of`](crate::HasBonds::bonds_of) run in O(degree); the conjugated
/// systems are perceived internally.
pub fn perceive<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>(
    mol: &M,
) -> Hybridizations {
    let conjugated = systems(mol);

    let mut aromatic: FxHashSet<SiteId> = FxHashSet::default();
    for bond in mol.bonds() {
        if mol.bond_order(bond) == BondOrder::Aromatic {
            let (a, b) = mol.bond_endpoints(bond);
            aromatic.insert(a);
            aromatic.insert(b);
        }
    }

    let labels = mol.sites().filter_map(|site| {
        if aromatic.contains(&site) {
            return Some((site, Hybridization::Sp2));
        }
        let domains = mol.degree(site) as u32 + lone_pairs(mol, site)?;
        let donated: u32 = conjugated
            .of_site(site)
            .map(|system| system.donated_pairs(site))
            .sum();
        let label = match domains - donated {
            0 | 1 => Hybridization::S,
            2 => Hybridization::Sp,
            3 => Hybridization::Sp2,
            4 => Hybridization::Sp3,
            5 => Hybridization::Sp3d,
            6 => Hybridization::Sp3d2,
            7 => Hybridization::Sp3d3,
            _ => Hybridization::Other,
        };
        Some((site, label))
    });

    Hybridizations {
        labels: SortedMap::from_pairs(labels),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::iter::once;

    use vita_core::Element;

    use crate::{BondId, HasBonds};

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

    fn molecule(atoms: &[(u32, &str, i8, u8)], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, ..)| s(id)).collect(),
            elements: atoms.iter().map(|&(_, symbol, ..)| elem(symbol)).collect(),
            formal_charges: atoms.iter().map(|&(_, _, charge, _)| charge).collect(),
            radicals: atoms.iter().map(|&(.., radicals)| radicals).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, u, v, _)| (s(u), s(v))).collect(),
            orders: bonds.iter().map(|&(.., order)| order).collect(),
        }
    }

    fn empty() -> Mol {
        molecule(&[], &[])
    }

    fn ethane() -> Mol {
        molecule(
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
                (1, 1, 2, BondOrder::Single),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 2, 6, BondOrder::Single),
                (6, 2, 7, BondOrder::Single),
                (7, 2, 8, BondOrder::Single),
            ],
        )
    }

    fn acetylene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Triple),
                (2, 1, 3, BondOrder::Single),
                (3, 2, 4, BondOrder::Single),
            ],
        )
    }

    fn ethylene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 2, 5, BondOrder::Single),
                (5, 2, 6, BondOrder::Single),
            ],
        )
    }

    fn water() -> Mol {
        molecule(
            &[(1, "O", 0, 0), (2, "H", 0, 0), (3, "H", 0, 0)],
            &[(1, 1, 2, BondOrder::Single), (2, 1, 3, BondOrder::Single)],
        )
    }

    fn formamide() -> Mol {
        molecule(
            &[
                (1, "O", 0, 0),
                (2, "C", 0, 0),
                (3, "N", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 2, 4, BondOrder::Single),
                (4, 3, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
            ],
        )
    }

    fn kekule_furan() -> Mol {
        molecule(
            &[
                (1, "O", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Double),
                (3, 3, 4, BondOrder::Single),
                (4, 4, 5, BondOrder::Double),
                (5, 5, 1, BondOrder::Single),
                (6, 2, 6, BondOrder::Single),
                (7, 3, 7, BondOrder::Single),
                (8, 4, 8, BondOrder::Single),
                (9, 5, 9, BondOrder::Single),
            ],
        )
    }

    fn aniline() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "N", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Aromatic),
                (2, 2, 3, BondOrder::Aromatic),
                (3, 3, 4, BondOrder::Aromatic),
                (4, 4, 5, BondOrder::Aromatic),
                (5, 5, 6, BondOrder::Aromatic),
                (6, 6, 1, BondOrder::Aromatic),
                (7, 1, 7, BondOrder::Single),
                (8, 7, 8, BondOrder::Single),
                (9, 7, 9, BondOrder::Single),
                (10, 2, 10, BondOrder::Single),
                (11, 3, 11, BondOrder::Single),
                (12, 4, 12, BondOrder::Single),
                (13, 5, 13, BondOrder::Single),
                (14, 6, 14, BondOrder::Single),
            ],
        )
    }

    fn chloroacetylene() -> Mol {
        molecule(
            &[
                (1, "Cl", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Triple),
                (3, 3, 4, BondOrder::Single),
            ],
        )
    }

    fn azide() -> Mol {
        molecule(
            &[(1, "N", -1, 0), (2, "N", 1, 0), (3, "N", -1, 0)],
            &[(1, 1, 2, BondOrder::Double), (2, 2, 3, BondOrder::Double)],
        )
    }

    fn dimethyl_sulfone() -> Mol {
        molecule(
            &[
                (1, "S", 0, 0),
                (2, "O", 0, 0),
                (3, "O", 0, 0),
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
                (1, 1, 2, BondOrder::Double),
                (2, 1, 3, BondOrder::Double),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 4, 6, BondOrder::Single),
                (6, 4, 7, BondOrder::Single),
                (7, 4, 8, BondOrder::Single),
                (8, 5, 9, BondOrder::Single),
                (9, 5, 10, BondOrder::Single),
                (10, 5, 11, BondOrder::Single),
            ],
        )
    }

    fn allyl(charge: i8, radicals: u8) -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", charge, radicals),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 2, 6, BondOrder::Single),
                (6, 3, 7, BondOrder::Single),
                (7, 3, 8, BondOrder::Single),
            ],
        )
    }

    fn butadiene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 3, 4, BondOrder::Double),
                (4, 1, 5, BondOrder::Single),
                (5, 1, 6, BondOrder::Single),
                (6, 2, 7, BondOrder::Single),
                (7, 3, 8, BondOrder::Single),
                (8, 4, 9, BondOrder::Single),
                (9, 4, 10, BondOrder::Single),
            ],
        )
    }

    fn iron_vinyl() -> Mol {
        molecule(
            &[
                (1, "Fe", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Double),
                (3, 2, 4, BondOrder::Single),
                (4, 3, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
            ],
        )
    }

    fn bridged_hydrogen() -> Mol {
        molecule(
            &[(1, "H", 0, 0), (2, "C", 0, 0), (3, "C", 0, 0)],
            &[(1, 1, 2, BondOrder::Single), (2, 1, 3, BondOrder::Single)],
        )
    }

    fn formaldehyde() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "O", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
            ],
        )
    }

    fn formamide_polar() -> Mol {
        molecule(
            &[
                (1, "O", -1, 0),
                (2, "C", 0, 0),
                (3, "N", 1, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Double),
                (3, 2, 4, BondOrder::Single),
                (4, 3, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
            ],
        )
    }

    #[test]
    fn an_empty_molecule_has_no_labels() {
        let perceived = perceive(&empty());
        assert_eq!(perceived.len(), 0);
        assert!(perceived.is_empty());
    }

    #[test]
    fn an_isolated_site_is_s() {
        let helium = molecule(&[(1, "He", 0, 0)], &[]);
        assert_eq!(
            perceive(&helium).hybridization(s(1)),
            Some(Hybridization::S)
        );
    }

    #[test]
    fn a_hydrogen_is_s() {
        assert_eq!(
            perceive(&ethane()).hybridization(s(3)),
            Some(Hybridization::S)
        );
    }

    #[test]
    fn an_alkyne_carbon_is_sp() {
        assert_eq!(
            perceive(&acetylene()).hybridization(s(1)),
            Some(Hybridization::Sp)
        );
    }

    #[test]
    fn an_alkene_carbon_is_sp2() {
        assert_eq!(
            perceive(&ethylene()).hybridization(s(1)),
            Some(Hybridization::Sp2)
        );
    }

    #[test]
    fn a_saturated_carbon_is_sp3() {
        let perceived = perceive(&ethane());
        assert_eq!(perceived.hybridization(s(1)), Some(Hybridization::Sp3));
        assert_eq!(perceived.hybridization(s(2)), Some(Hybridization::Sp3));
    }

    #[test]
    fn a_water_oxygen_is_sp3() {
        assert_eq!(
            perceive(&water()).hybridization(s(1)),
            Some(Hybridization::Sp3)
        );
    }

    #[test]
    fn expanded_octets_map_to_d_labels() {
        for (center, neighbors, label) in [
            ("P", 5, Hybridization::Sp3d),
            ("S", 6, Hybridization::Sp3d2),
            ("I", 7, Hybridization::Sp3d3),
            ("Xe", 8, Hybridization::Other),
        ] {
            let atoms: Vec<(u32, &str, i8, u8)> = once((1, center, 0, 0))
                .chain((2..=neighbors + 1).map(|id| (id, "F", 0, 0)))
                .collect();
            let bonds: Vec<(u32, u32, u32, BondOrder)> = (2..=neighbors + 1)
                .map(|id| (id - 1, 1, id, BondOrder::Single))
                .collect();
            let fluoride = molecule(&atoms, &bonds);
            assert_eq!(perceive(&fluoride).hybridization(s(1)), Some(label));
        }
    }

    #[test]
    fn an_amide_nitrogen_is_sp2() {
        let perceived = perceive(&formamide());
        assert_eq!(perceived.hybridization(s(3)), Some(Hybridization::Sp2));
        assert_eq!(perceived.hybridization(s(2)), Some(Hybridization::Sp2));
        assert_eq!(perceived.hybridization(s(1)), Some(Hybridization::Sp2));
    }

    #[test]
    fn a_kekule_heteroring_donor_is_sp2() {
        assert_eq!(
            perceive(&kekule_furan()).hybridization(s(1)),
            Some(Hybridization::Sp2)
        );
    }

    #[test]
    fn an_aniline_nitrogen_is_sp2() {
        assert_eq!(
            perceive(&aniline()).hybridization(s(7)),
            Some(Hybridization::Sp2)
        );
    }

    #[test]
    fn a_declared_aromatic_site_is_sp2() {
        assert_eq!(
            perceive(&aniline()).hybridization(s(1)),
            Some(Hybridization::Sp2)
        );
    }

    #[test]
    fn a_donor_into_both_planes_is_sp() {
        assert_eq!(
            perceive(&chloroacetylene()).hybridization(s(1)),
            Some(Hybridization::Sp)
        );
    }

    #[test]
    fn an_azide_is_linear_throughout() {
        let perceived = perceive(&azide());
        for site in [s(1), s(2), s(3)] {
            assert_eq!(perceived.hybridization(site), Some(Hybridization::Sp));
        }
    }

    #[test]
    fn a_hypervalent_center_keeps_its_domain_count() {
        let perceived = perceive(&dimethyl_sulfone());
        assert_eq!(perceived.hybridization(s(1)), Some(Hybridization::Sp3));
        assert_eq!(perceived.hybridization(s(2)), Some(Hybridization::Sp2));
    }

    #[test]
    fn a_radical_or_carbenium_center_is_sp2() {
        let radical = perceive(&allyl(0, 1));
        let cation = perceive(&allyl(1, 0));
        assert_eq!(radical.hybridization(s(3)), Some(Hybridization::Sp2));
        assert_eq!(cation.hybridization(s(3)), Some(Hybridization::Sp2));
    }

    #[test]
    fn a_donation_free_molecule_keeps_plain_domain_labels() {
        let perceived = perceive(&butadiene());
        for site in [s(1), s(2), s(3), s(4)] {
            assert_eq!(perceived.hybridization(site), Some(Hybridization::Sp2));
        }
        for site in (5..=10).map(s) {
            assert_eq!(perceived.hybridization(site), Some(Hybridization::S));
        }
    }

    #[test]
    fn a_d_block_site_is_unlabeled() {
        let perceived = perceive(&iron_vinyl());
        assert_eq!(perceived.hybridization(s(1)), None);
        assert_eq!(perceived.hybridization(s(2)), Some(Hybridization::Sp2));
    }

    #[test]
    fn an_overbonded_hydrogen_is_unlabeled() {
        assert_eq!(perceive(&bridged_hydrogen()).hybridization(s(1)), None);
    }

    #[test]
    fn the_table_reports_and_iterates_its_labels() {
        let perceived = perceive(&formaldehyde());
        assert_eq!(perceived.len(), 4);
        assert!(!perceived.is_empty());
        assert_eq!(perceived.hybridization(s(99)), None);
        assert_eq!(
            perceived.hybridizations().collect::<Vec<_>>(),
            vec![
                (s(1), Hybridization::Sp2),
                (s(2), Hybridization::Sp2),
                (s(3), Hybridization::S),
                (s(4), Hybridization::S),
            ]
        );
    }

    #[test]
    fn the_table_holds_only_the_labeled_sites() {
        let mol = iron_vinyl();
        let perceived = perceive(&mol);
        assert_eq!(perceived.len(), mol.sites().count() - 1);
        assert!(perceived.hybridizations().all(|(site, _)| site != s(1)));
    }

    #[test]
    fn the_bound_view_answers_the_capability() {
        let mol = iron_vinyl();
        let perceived = perceive(&mol);
        let viewed = perceived.bind(&mol);
        assert_eq!(viewed.hybridization(s(2)), Hybridization::Sp2);
        assert_eq!(viewed.hybridization(s(1)), Hybridization::Other);
    }

    #[test]
    fn the_bound_view_forwards_the_skeleton() {
        let mol = iron_vinyl();
        let perceived = perceive(&mol);
        let viewed = perceived.bind(&mol);
        assert_eq!(viewed.element(s(1)), elem("Fe"));
        assert_eq!(viewed.bond_order(b(2)), BondOrder::Double);
    }

    #[test]
    fn the_labels_are_independent_of_the_resonance_form() {
        assert_eq!(perceive(&formamide()), perceive(&formamide_polar()));
    }

    #[test]
    fn the_labels_are_independent_of_input_order() {
        let shuffled = molecule(
            &[
                (8, "H", 0, 0),
                (3, "C", 0, 0),
                (10, "H", 0, 0),
                (1, "C", 0, 0),
                (6, "H", 0, 0),
                (4, "C", 0, 0),
                (9, "H", 0, 0),
                (2, "C", 0, 0),
                (5, "H", 0, 0),
                (7, "H", 0, 0),
            ],
            &[
                (7, 3, 8, BondOrder::Single),
                (3, 3, 4, BondOrder::Double),
                (9, 4, 10, BondOrder::Single),
                (1, 1, 2, BondOrder::Double),
                (5, 1, 6, BondOrder::Single),
                (8, 4, 9, BondOrder::Single),
                (2, 2, 3, BondOrder::Single),
                (6, 2, 7, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
            ],
        );
        assert_eq!(perceive(&butadiene()), perceive(&shuffled));
    }
}

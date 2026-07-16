use std::collections::BTreeSet;

use vita_core::SiteId;

use super::stereocenters;
use crate::{BondId, HasBondOrders, HasStereoConfigurations, StereoKind, StereoLocus};

/// The disagreements between a molecule's stereocentres and its declared
/// configurations.
///
/// A stereogenic locus with no configuration is *unspecified*; a configuration at a
/// locus that is not stereogenic is *overspecified*.
///
/// Obtain via [`consistency`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoConsistency {
    unspecified: Vec<StereoLocus>,
    overspecified: Vec<StereoLocus>,
}

impl StereoConsistency {
    /// Returns `true` if `locus` is stereogenic yet carries no configuration.
    pub fn contains_unspecified(&self, locus: StereoLocus) -> bool {
        self.unspecified.binary_search(&locus).is_ok()
    }

    /// Returns `true` if `locus` carries a configuration yet is not stereogenic.
    pub fn contains_overspecified(&self, locus: StereoLocus) -> bool {
        self.overspecified.binary_search(&locus).is_ok()
    }

    /// Number of stereogenic loci left without a configuration.
    pub fn unspecified_count(&self) -> usize {
        self.unspecified.len()
    }

    /// Number of loci carrying a configuration yet not stereogenic.
    pub fn overspecified_count(&self) -> usize {
        self.overspecified.len()
    }

    /// Iterates the stereogenic loci left without a configuration, ascending.
    pub fn unspecified(&self) -> impl Iterator<Item = StereoLocus> + '_ {
        self.unspecified.iter().copied()
    }

    /// Iterates the loci that carry a configuration yet are not stereogenic,
    /// ascending.
    pub fn overspecified(&self) -> impl Iterator<Item = StereoLocus> + '_ {
        self.overspecified.iter().copied()
    }

    /// Returns `true` if every stereocentre is specified and no configuration is
    /// spurious.
    pub fn is_consistent(&self) -> bool {
        self.unspecified.is_empty() && self.overspecified.is_empty()
    }
}

/// Reconciles a molecule's stereocentres with the configurations it declares.
///
/// Detects the stereogenic loci under the caller's coloring and `candidate` (as
/// [`stereocenters`]), then compares them with the loci [`stereo_configurations`]
/// speak for: a stereocentre without a configuration is reported unspecified, a
/// configuration off a stereocentre overspecified.
///
/// [`stereo_configurations`]: HasStereoConfigurations::stereo_configurations
///
/// # Complexity
///
/// As [`stereocenters`]: O((V + L) · V · (V + E) · log V) time and O(V + E) space, over
/// the molecule's `V` sites and `E` bonds, for `L` candidate loci.
pub fn consistency<M, VK, EK>(
    mol: &M,
    site_key: impl Fn(SiteId) -> VK,
    bond_key: impl Fn(BondId) -> EK,
    candidate: impl Fn(StereoLocus) -> Option<StereoKind>,
) -> StereoConsistency
where
    M: HasStereoConfigurations + HasBondOrders,
    VK: Ord,
    EK: Ord,
{
    let detected: BTreeSet<StereoLocus> = stereocenters(mol, site_key, bond_key, candidate)
        .iter()
        .collect();
    let declared: BTreeSet<StereoLocus> = mol
        .stereo_configurations()
        .map(|config| config.locus())
        .collect();

    StereoConsistency {
        unspecified: detected.difference(&declared).copied().collect(),
        overspecified: declared.difference(&detected).copied().collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondOrder::Single;
    use crate::{BondOrder, HasBonds, StereoConfiguration};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn centers_at(sites: &'static [u32]) -> impl Fn(StereoLocus) -> Option<StereoKind> {
        move |locus| match locus {
            StereoLocus::Site(site) if sites.contains(&site.get()) => Some(StereoKind::Tetrahedral),
            _ => None,
        }
    }

    fn config(site: u32, order: [u32; 4]) -> StereoConfiguration {
        StereoConfiguration::new(
            StereoLocus::Site(s(site)),
            StereoKind::Tetrahedral,
            order.map(s),
        )
        .unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        colors: Vec<u32>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
        configs: Vec<StereoConfiguration>,
    }

    impl Mol {
        fn color(&self, site: SiteId) -> u32 {
            self.colors[self.sites.iter().position(|&x| x == site).unwrap()]
        }
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
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

    impl HasStereoConfigurations for Mol {
        fn stereo_configurations(&self) -> impl Iterator<Item = StereoConfiguration> + '_ {
            self.configs.iter().cloned()
        }
    }

    fn reconcile(
        mol: &Mol,
        candidate: impl Fn(StereoLocus) -> Option<StereoKind>,
    ) -> StereoConsistency {
        consistency(mol, |site| mol.color(site), |_| 0u32, candidate)
    }

    fn mol(
        atoms: &[(u32, u32)],
        bonds: &[(u32, u32, u32)],
        configs: Vec<StereoConfiguration>,
    ) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, _)| s(id)).collect(),
            colors: atoms.iter().map(|&(_, color)| color).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c)| (s(a), s(c))).collect(),
            orders: bonds.iter().map(|_| Single).collect(),
            configs,
        }
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            colors: m.colors.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
            orders: m.orders.iter().rev().copied().collect(),
            configs: m.configs.iter().rev().cloned().collect(),
        }
    }

    fn empty() -> Mol {
        mol(&[], &[], Vec::new())
    }

    fn specified() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
            vec![config(1, [2, 3, 4, 5])],
        )
    }

    fn unspecified() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
            Vec::new(),
        )
    }

    fn overspecified() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 1), (4, 2), (5, 3)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
            vec![config(1, [2, 3, 4, 5])],
        )
    }

    fn mixed() -> Mol {
        mol(
            &[
                (1, 0),
                (2, 1),
                (3, 2),
                (4, 3),
                (5, 4),
                (6, 0),
                (7, 5),
                (8, 5),
                (9, 6),
                (10, 7),
            ],
            &[
                (1, 1, 2),
                (2, 1, 3),
                (3, 1, 4),
                (4, 1, 5),
                (5, 6, 7),
                (6, 6, 8),
                (7, 6, 9),
                (8, 6, 10),
            ],
            vec![config(6, [7, 8, 9, 10])],
        )
    }

    #[test]
    fn empty_molecule_is_consistent() {
        let report = reconcile(&empty(), centers_at(&[1]));
        assert!(report.is_consistent());
        assert_eq!(report.unspecified_count(), 0);
        assert_eq!(report.overspecified_count(), 0);
    }

    #[test]
    fn a_configured_stereocenter_is_consistent() {
        assert!(reconcile(&specified(), centers_at(&[1])).is_consistent());
    }

    #[test]
    fn a_stereocenter_without_a_configuration_is_unspecified() {
        let report = reconcile(&unspecified(), centers_at(&[1]));
        assert!(report.contains_unspecified(StereoLocus::Site(s(1))));
        assert!(!report.contains_overspecified(StereoLocus::Site(s(1))));
        assert_eq!(report.unspecified_count(), 1);
        assert!(!report.is_consistent());
    }

    #[test]
    fn a_configuration_off_a_stereocenter_is_overspecified() {
        let report = reconcile(&overspecified(), centers_at(&[1]));
        assert!(report.contains_overspecified(StereoLocus::Site(s(1))));
        assert!(!report.contains_unspecified(StereoLocus::Site(s(1))));
        assert_eq!(report.overspecified_count(), 1);
        assert!(!report.is_consistent());
    }

    #[test]
    fn unspecified_lists_the_stereocenters_without_a_configuration() {
        let report = reconcile(&unspecified(), centers_at(&[1]));
        let loci: Vec<StereoLocus> = report.unspecified().collect();
        assert_eq!(loci, vec![StereoLocus::Site(s(1))]);
    }

    #[test]
    fn overspecified_lists_the_configurations_off_a_stereocenter() {
        let report = reconcile(&overspecified(), centers_at(&[1]));
        let loci: Vec<StereoLocus> = report.overspecified().collect();
        assert_eq!(loci, vec![StereoLocus::Site(s(1))]);
    }

    #[test]
    fn both_kinds_of_disagreement_are_reported_together() {
        let report = reconcile(&mixed(), centers_at(&[1, 6]));
        assert!(report.contains_unspecified(StereoLocus::Site(s(1))));
        assert!(report.contains_overspecified(StereoLocus::Site(s(6))));
        assert!(!report.is_consistent());
    }

    #[test]
    fn reconciliation_is_independent_of_input_order() {
        let molecule = mixed();
        assert_eq!(
            reconcile(&molecule, centers_at(&[1, 6])),
            reconcile(&reversed(&molecule), centers_at(&[1, 6])),
        );
    }
}

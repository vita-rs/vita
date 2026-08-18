use vita_core::SiteId;

use super::{candidate_loci, realizable, refined_key, settle};
use crate::{BondId, HasBondOrders, HasStereoConfigurations, StereoKind, StereoLocus};

/// The stereogenic units of a molecule: the loci whose substituents symmetry cannot
/// interchange, so more than one configuration is realizable.
///
/// Membership is a function of the graph, the caller's coloring, and the
/// configurations the molecule declares — the last because a pseudo-asymmetric
/// center is stereogenic only for a particular pairing of its neighbors' handedness.
/// A locus survives when the substituent classes its geometry admits realize more
/// than one distinct arrangement, so repeated substituents are weighed correctly: an
/// octahedral MA₄B₂ center is stereogenic (cis and trans), an MA₅B one is not.
///
/// Obtain via [`stereocenters`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Stereocenters {
    loci: Vec<StereoLocus>,
}

impl Stereocenters {
    /// The number of stereogenic units.
    pub fn len(&self) -> usize {
        self.loci.len()
    }

    /// Returns `true` if the molecule has no stereogenic units.
    pub fn is_empty(&self) -> bool {
        self.loci.is_empty()
    }

    /// Returns `true` if `locus` is stereogenic.
    ///
    /// Returns `false` if `locus` is absent from the molecule or is not
    /// stereogenic.
    pub fn contains(&self, locus: StereoLocus) -> bool {
        self.loci.binary_search(&locus).is_ok()
    }

    /// Iterates the stereogenic loci in ascending order: sites, then bonds, then axes.
    pub fn iter(&self) -> impl Iterator<Item = StereoLocus> + '_ {
        self.loci.iter().copied()
    }
}

/// The stereogenic loci of a molecule.
///
/// `candidate` is the caller's: it reports which loci could bear stereochemistry and
/// of which [`StereoKind`] — a center's coordination geometry — which the library
/// cannot know. The library then keeps only those a symmetry does *not* neutralise:
/// `site_key` and `bond_key` color the graph as for [`canonicalize`](crate::canonical::canonicalize),
/// the declared configurations refine that coloring to a fixpoint — so a center whose
/// stereogenicity turns on its neighbors' handedness is resolved — and a candidate
/// survives when, once its frame is individualised, its substituent classes admit
/// more than one arrangement under the geometry's group.
///
/// A molecule that declares no configurations is colored by its constitution alone,
/// and detection is then purely topological. The substituents are taken from the
/// bonds present; a stereocenter bearing an implicit hydrogen must give it
/// explicitly, or fold the coordination into `site_key`. A bond's or axis's
/// substituents hang off the termini of its rigid double-bond chain, so a plain
/// double bond and a cumulene resolve alike.
///
/// Detection speaks for the molecule as declared. Which loci could be stereogenic
/// under *some* assignment is a different, intrinsically enumerative question —
/// answered by the loci [`stereoisomers`](super::stereoisomers) assigns.
///
/// # Complexity
///
/// O((V + L) · V · (V + E) · log V) time and O(V + E) space, over the molecule's `V`
/// sites and `E` bonds, for `L` candidate loci — the fixpoint refinement and one
/// [`canonicalize`](crate::canonical::canonicalize) per candidate.
pub fn stereocenters<M, VK, EK>(
    mol: &M,
    site_key: impl Fn(SiteId) -> VK,
    bond_key: impl Fn(BondId) -> EK,
    candidate: impl Fn(StereoLocus) -> Option<StereoKind>,
) -> Stereocenters
where
    M: HasStereoConfigurations + HasBondOrders,
    VK: Ord,
    EK: Ord,
{
    let signal = settle(mol, &site_key, &bond_key);
    let key = refined_key(&signal, &site_key);

    let surviving = |locus| {
        let kind = candidate(locus)?;
        (realizable(mol, locus, kind, &key, &bond_key).len() > 1).then_some(locus)
    };
    let mut loci: Vec<StereoLocus> = candidate_loci(mol).filter_map(surviving).collect();
    loci.sort_unstable();
    Stereocenters { loci }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondOrder::{Double, Single};
    use crate::CoordinationGeometry::*;
    use crate::{
        BondOrder, CoordinationGeometry, HasBonds, StereoConfiguration, StereogenicGeometry,
    };

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn center(geometry: CoordinationGeometry) -> StereoKind {
        StereoKind::Center(StereogenicGeometry::new(geometry).expect("the geometry is stereogenic"))
    }

    fn only(target: StereoLocus, kind: StereoKind) -> impl Fn(StereoLocus) -> Option<StereoKind> {
        move |locus| (locus == target).then_some(kind)
    }

    fn centers_at(sites: &'static [u32]) -> impl Fn(StereoLocus) -> Option<StereoKind> {
        move |locus| match locus {
            StereoLocus::Site(site) if sites.contains(&site.get()) => Some(center(Tetrahedral)),
            _ => None,
        }
    }

    fn config(site: u32, order: [u32; 4]) -> StereoConfiguration {
        StereoConfiguration::new(
            StereoLocus::Site(s(site)),
            center(Tetrahedral),
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

        fn with(mut self, configs: impl IntoIterator<Item = StereoConfiguration>) -> Self {
            self.configs = configs.into_iter().collect();
            self
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

    fn detect(mol: &Mol, candidate: impl Fn(StereoLocus) -> Option<StereoKind>) -> Stereocenters {
        stereocenters(mol, |site| mol.color(site), |_| 0u32, candidate)
    }

    fn mol(atoms: &[(u32, u32)], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, _)| s(id)).collect(),
            colors: atoms.iter().map(|&(_, color)| color).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c, _)| (s(a), s(c))).collect(),
            orders: bonds.iter().map(|&(_, _, _, order)| order).collect(),
            configs: Vec::new(),
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
        mol(&[], &[])
    }

    fn tetrahedral() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
            ],
        )
    }

    fn two_like_substituents() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 1), (4, 2), (5, 3)],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
            ],
        )
    }

    fn three_substituents() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 2), (4, 3)],
            &[(1, 1, 2, Single), (2, 1, 3, Single), (3, 1, 4, Single)],
        )
    }

    fn alkene() -> Mol {
        mol(
            &[(1, 0), (2, 0), (3, 1), (4, 2), (5, 3), (6, 4)],
            &[
                (1, 1, 2, Double),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 2, 5, Single),
                (5, 2, 6, Single),
            ],
        )
    }

    fn symmetric_alkene() -> Mol {
        mol(
            &[(1, 0), (2, 0), (3, 1), (4, 1), (5, 2), (6, 3)],
            &[
                (1, 1, 2, Double),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 2, 5, Single),
                (5, 2, 6, Single),
            ],
        )
    }

    fn allene() -> Mol {
        mol(
            &[(1, 0), (2, 5), (3, 0), (4, 1), (5, 2), (6, 3), (7, 4)],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Double),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
                (5, 3, 6, Single),
                (6, 3, 7, Single),
            ],
        )
    }

    fn octahedral_m_a4_b2() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 1), (4, 1), (5, 1), (6, 2), (7, 2)],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
                (5, 1, 6, Single),
                (6, 1, 7, Single),
            ],
        )
    }

    fn octahedral_m_a5_b() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (7, 2)],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
                (5, 1, 6, Single),
                (6, 1, 7, Single),
            ],
        )
    }

    fn paired_centers() -> Mol {
        mol(
            &[
                (1, 0),
                (2, 0),
                (3, 1),
                (4, 2),
                (5, 3),
                (6, 1),
                (7, 2),
                (8, 3),
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

    fn trihydroxyglutaric(arms: [[u32; 4]; 2]) -> Mol {
        mol(
            &[
                (1, 0),
                (2, 0),
                (3, 0),
                (4, 1),
                (5, 2),
                (6, 1),
                (7, 2),
                (8, 1),
                (9, 2),
                (10, 3),
                (11, 3),
            ],
            &[
                (1, 1, 2, Single),
                (2, 2, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
                (5, 1, 10, Single),
                (6, 2, 6, Single),
                (7, 2, 7, Single),
                (8, 3, 8, Single),
                (9, 3, 9, Single),
                (10, 3, 11, Single),
            ],
        )
        .with([
            config(1, arms[0]),
            config(2, [1, 6, 7, 3]),
            config(3, arms[1]),
        ])
    }

    #[test]
    fn empty_molecule_has_no_stereocenters() {
        let centers = detect(&empty(), |_| Some(center(Tetrahedral)));
        assert_eq!(centers.len(), 0);
        assert!(centers.is_empty());
    }

    #[test]
    fn a_molecule_the_candidate_rejects_has_no_stereocenters() {
        assert!(detect(&tetrahedral(), |_| None).is_empty());
    }

    #[test]
    fn a_carbon_with_four_distinct_substituents_is_a_stereocenter() {
        let centers = detect(
            &tetrahedral(),
            only(StereoLocus::Site(s(1)), center(Tetrahedral)),
        );
        assert!(centers.contains(StereoLocus::Site(s(1))));
    }

    #[test]
    fn a_double_bond_with_distinct_ends_is_a_stereocenter() {
        let centers = detect(&alkene(), only(StereoLocus::Bond(b(1)), StereoKind::Bond));
        assert!(centers.contains(StereoLocus::Bond(b(1))));
    }

    #[test]
    fn an_allene_with_distinct_termini_is_a_stereocenter() {
        let centers = detect(&allene(), only(StereoLocus::Axis(s(2)), StereoKind::Axis));
        assert!(centers.contains(StereoLocus::Axis(s(2))));
    }

    #[test]
    fn an_octahedral_m_a4_b2_center_is_a_stereocenter() {
        let centers = detect(
            &octahedral_m_a4_b2(),
            only(StereoLocus::Site(s(1)), center(Octahedral)),
        );
        assert!(centers.contains(StereoLocus::Site(s(1))));
    }

    #[test]
    fn a_carbon_with_two_like_substituents_is_not_a_stereocenter() {
        let centers = detect(
            &two_like_substituents(),
            only(StereoLocus::Site(s(1)), center(Tetrahedral)),
        );
        assert!(!centers.contains(StereoLocus::Site(s(1))));
    }

    #[test]
    fn a_double_bond_with_a_symmetric_end_is_not_a_stereocenter() {
        let centers = detect(
            &symmetric_alkene(),
            only(StereoLocus::Bond(b(1)), StereoKind::Bond),
        );
        assert!(!centers.contains(StereoLocus::Bond(b(1))));
    }

    #[test]
    fn an_octahedral_m_a5_b_center_is_not_a_stereocenter() {
        let centers = detect(
            &octahedral_m_a5_b(),
            only(StereoLocus::Site(s(1)), center(Octahedral)),
        );
        assert!(!centers.contains(StereoLocus::Site(s(1))));
    }

    #[test]
    fn a_candidate_off_its_anchor_is_not_a_stereocenter() {
        let centers = detect(
            &alkene(),
            only(StereoLocus::Bond(b(1)), center(Tetrahedral)),
        );
        assert!(centers.is_empty());
    }

    #[test]
    fn a_center_short_of_its_substituents_is_not_a_stereocenter() {
        let centers = detect(
            &three_substituents(),
            only(StereoLocus::Site(s(1)), center(Tetrahedral)),
        );
        assert!(!centers.contains(StereoLocus::Site(s(1))));
    }

    #[test]
    fn symmetry_equivalent_centers_are_each_detected() {
        let centers = detect(&paired_centers(), centers_at(&[1, 2]));
        assert!(centers.contains(StereoLocus::Site(s(1))));
        assert!(centers.contains(StereoLocus::Site(s(2))));
    }

    #[test]
    fn a_pseudo_asymmetric_center_is_detected() {
        let mol = trihydroxyglutaric([[2, 4, 5, 10], [8, 2, 9, 11]]);
        let centers = detect(&mol, only(StereoLocus::Site(s(2)), center(Tetrahedral)));
        assert!(centers.contains(StereoLocus::Site(s(2))));
    }

    #[test]
    fn a_dormant_pseudo_center_is_not_detected() {
        let mol = trihydroxyglutaric([[2, 4, 5, 10], [2, 8, 9, 11]]);
        let centers = detect(&mol, only(StereoLocus::Site(s(2)), center(Tetrahedral)));
        assert!(!centers.contains(StereoLocus::Site(s(2))));
    }

    #[test]
    fn count_reports_the_number_of_stereocenters() {
        let centers = detect(&paired_centers(), centers_at(&[1, 2]));
        assert_eq!(centers.len(), 2);
        assert!(!centers.is_empty());
    }

    #[test]
    fn contains_is_false_for_an_absent_locus() {
        let centers = detect(
            &tetrahedral(),
            only(StereoLocus::Site(s(1)), center(Tetrahedral)),
        );
        assert!(!centers.contains(StereoLocus::Site(s(99))));
    }

    #[test]
    fn iter_yields_the_stereocenters_in_ascending_order() {
        let centers = detect(&paired_centers(), centers_at(&[1, 2]));
        let loci: Vec<StereoLocus> = centers.iter().collect();
        assert_eq!(loci, vec![StereoLocus::Site(s(1)), StereoLocus::Site(s(2))]);
    }

    #[test]
    fn detection_is_independent_of_input_order() {
        let candidate = centers_at(&[1, 2]);
        let molecule = paired_centers();
        assert_eq!(
            detect(&molecule, &candidate),
            detect(&reversed(&molecule), &candidate),
        );
    }
}

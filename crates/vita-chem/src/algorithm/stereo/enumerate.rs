use vita_core::SiteId;

use super::{
    StereoConfigurations, StereoForm, candidate_loci, form, realizable, refined_key, settle,
};
use crate::{BondId, HasBondOrders, StereoConfiguration, StereoKind, StereoLocus};

/// The stereoisomers of a molecule's constitution — every assignment of a
/// configuration to its stereogenic units that no symmetry equates.
///
/// Each isomer is the set of [`StereoConfigurations`] it bears, over one shared
/// constitution; bind it to the molecule with [`StereoConfigurations::bind`] to read
/// it as a whole. Enantiomers appear as distinct isomers, a meso form once, and a
/// pseudo-asymmetric center — stereogenic only for a particular pairing of its
/// neighbors' handedness — is enumerated exactly where it is realized.
///
/// Obtain via [`stereoisomers`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Stereoisomers {
    isomers: Vec<StereoConfigurations>,
}

impl Stereoisomers {
    /// The number of distinct stereoisomers.
    pub fn len(&self) -> usize {
        self.isomers.len()
    }

    /// Returns `true` if the constitution has no stereoisomers — it never does; a
    /// molecule with no stereocenters is its own sole isomer.
    pub fn is_empty(&self) -> bool {
        self.isomers.is_empty()
    }

    /// Iterates the stereoisomers in a canonical order, each the configurations it
    /// bears.
    pub fn iter(&self) -> impl Iterator<Item = &StereoConfigurations> + '_ {
        self.isomers.iter()
    }
}

/// Enumerates the stereoisomers a molecule's constitution admits.
///
/// `candidate`, `site_key`, and `bond_key` are the caller's, exactly as for
/// [`stereocenters`](super::stereocenters): the first reports which loci could bear
/// stereochemistry and of which [`StereoKind`], the latter two color the graph. The
/// stereogenic units are found under the configurations assigned so far and one is
/// given each of its realizable configurations in turn, recursively — so a
/// pseudo-asymmetric center surfaces once the neighbors it depends on are fixed.
/// Assignments the molecule's symmetry equates, and the two halves of a meso form,
/// collapse to one isomer by their [`StereoForm`].
///
/// The molecule's own declared configurations, if any, are ignored: the enumeration
/// speaks for the constitution alone. Every locus that can ever be stereogenic is
/// assigned in some isomer, so the assigned loci, unioned over the isomers, are the
/// constitution's complete potential stereocenters — pseudo-asymmetric ones included.
///
/// # Complexity
///
/// O(N · (V + L) · V · (V + E) · log V) time and O(N · k + V + E) space, over the
/// molecule's `V` sites and `E` bonds, for `L` candidate loci and `N` assignments of
/// configurations to the `k` stereocenters — up to the product of their configuration
/// counts, so exponential in `k`. Each search node runs a stereogenic detection as
/// costly as [`stereocenters`](super::stereocenters); every assignment is held to fold
/// symmetric and meso duplicates by their [`StereoForm`].
pub fn stereoisomers<M, VK, EK>(
    mol: &M,
    site_key: impl Fn(SiteId) -> VK,
    bond_key: impl Fn(BondId) -> EK,
    candidate: impl Fn(StereoLocus) -> Option<StereoKind>,
) -> Stereoisomers
where
    M: HasBondOrders,
    VK: Ord,
    EK: Ord,
{
    let mut assignments: Vec<Vec<StereoConfiguration>> = Vec::new();
    extend(
        mol,
        Vec::new(),
        &site_key,
        &bond_key,
        &candidate,
        &mut assignments,
    );

    let mut isomers: Vec<(StereoForm<VK, EK>, StereoConfigurations)> = assignments
        .into_iter()
        .map(|assignment| {
            let isomer = StereoConfigurations::from_configurations(assignment);
            let identity = form(&isomer.bind(mol), &site_key, &bond_key);
            (identity, isomer)
        })
        .collect();
    isomers.sort_by(|a, b| a.0.cmp(&b.0));
    isomers.dedup_by(|a, b| a.0 == b.0);

    Stereoisomers {
        isomers: isomers.into_iter().map(|(_, isomer)| isomer).collect(),
    }
}

/// Extends a partial assignment: gives the smallest stereogenic unit it leaves open
/// each of its realizable configurations in turn, or records the assignment once none
/// remains.
fn extend<M, VK, EK>(
    mol: &M,
    assigned: Vec<StereoConfiguration>,
    site_key: &impl Fn(SiteId) -> VK,
    bond_key: &impl Fn(BondId) -> EK,
    candidate: &impl Fn(StereoLocus) -> Option<StereoKind>,
    out: &mut Vec<Vec<StereoConfiguration>>,
) where
    M: HasBondOrders,
    VK: Ord,
    EK: Ord,
{
    let assigned_so_far = StereoConfigurations::from_configurations(assigned.clone());
    let view = assigned_so_far.bind(mol);
    let signal = settle(&view, site_key, bond_key);
    let key = refined_key(&signal, site_key);

    let open = |locus: StereoLocus| -> Option<Vec<StereoConfiguration>> {
        if assigned.iter().any(|config| config.locus() == locus) {
            return None;
        }
        let configs = realizable(&view, locus, candidate(locus)?, &key, bond_key);
        (configs.len() > 1).then_some(configs)
    };

    let next = candidate_loci(&view)
        .filter_map(|locus| Some((locus, open(locus)?)))
        .min_by_key(|(locus, _)| *locus);

    match next {
        None => out.push(assigned),
        Some((_, configs)) => {
            for config in configs {
                let mut branch = assigned.clone();
                branch.push(config);
                extend(mol, branch, site_key, bond_key, candidate, out);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondOrder::Single;
    use crate::CoordinationGeometry::*;
    use crate::{
        BondOrder, CoordinationGeometry, HasBonds, HasStereoConfigurations, StereogenicGeometry,
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

    fn centers_at(sites: &'static [u32]) -> impl Fn(StereoLocus) -> Option<StereoKind> {
        move |locus| match locus {
            StereoLocus::Site(site) if sites.contains(&site.get()) => Some(center(Tetrahedral)),
            _ => None,
        }
    }

    struct Mol {
        sites: Vec<SiteId>,
        colors: Vec<u32>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
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
        fn bond_order(&self, _: BondId) -> BondOrder {
            Single
        }
    }

    fn enumerate(
        mol: &Mol,
        candidate: impl Fn(StereoLocus) -> Option<StereoKind>,
    ) -> Stereoisomers {
        stereoisomers(mol, |site| mol.color(site), |_| 0u32, candidate)
    }

    fn mol(atoms: &[(u32, u32)], bonds: &[(u32, u32, u32)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, _)| s(id)).collect(),
            colors: atoms.iter().map(|&(_, color)| color).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c)| (s(a), s(c))).collect(),
        }
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            colors: m.colors.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
        }
    }

    fn tetrahedral() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
        )
    }

    fn tartaric() -> Mol {
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
                (1, 1, 2),
                (2, 1, 3),
                (3, 1, 4),
                (4, 1, 5),
                (5, 2, 6),
                (6, 2, 7),
                (7, 2, 8),
            ],
        )
    }

    fn trihydroxyglutaric() -> Mol {
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
                (1, 1, 2),
                (2, 2, 3),
                (3, 1, 4),
                (4, 1, 5),
                (5, 1, 10),
                (6, 2, 6),
                (7, 2, 7),
                (8, 3, 8),
                (9, 3, 9),
                (10, 3, 11),
            ],
        )
    }

    #[test]
    fn a_constitution_without_stereocenters_has_a_single_isomer() {
        assert_eq!(enumerate(&tetrahedral(), |_| None).len(), 1);
    }

    #[test]
    fn a_single_stereocenter_has_two_isomers() {
        assert_eq!(enumerate(&tetrahedral(), centers_at(&[1])).len(), 2);
    }

    #[test]
    fn two_equivalent_centers_have_three_stereoisomers() {
        assert_eq!(enumerate(&tartaric(), centers_at(&[1, 2])).len(), 3);
    }

    #[test]
    fn a_pseudo_asymmetric_constitution_has_four_stereoisomers() {
        assert_eq!(
            enumerate(&trihydroxyglutaric(), centers_at(&[1, 2, 3])).len(),
            4,
        );
    }

    #[test]
    fn a_lone_isomer_is_not_empty() {
        assert!(!enumerate(&tetrahedral(), |_| None).is_empty());
    }

    #[test]
    fn each_isomer_binds_a_configuration_to_the_stereocenter() {
        let molecule = tetrahedral();
        let isomers = enumerate(&molecule, centers_at(&[1]));
        for isomer in isomers.iter() {
            assert_eq!(isomer.bind(&molecule).stereo_configuration_count(), 1);
        }
    }

    #[test]
    fn enumeration_is_independent_of_input_order() {
        let molecule = trihydroxyglutaric();
        assert_eq!(
            enumerate(&molecule, centers_at(&[1, 2, 3])),
            enumerate(&reversed(&molecule), centers_at(&[1, 2, 3])),
        );
    }
}

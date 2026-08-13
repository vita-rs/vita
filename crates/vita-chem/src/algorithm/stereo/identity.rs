use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

use vita_core::SiteId;

use super::{Configured, refined};
use crate::algorithm::canonical::{Canonical, canonicalize};
use crate::{BondId, HasStereoConfigurations};

/// How two molecules relate stereochemically.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StereoRelationship {
    /// The same stereoisomer.
    Identical,
    /// Non-superimposable mirror images.
    Enantiomers,
    /// Stereoisomers that are not mirror images.
    Diastereomers,
    /// Not stereoisomers — their constitutions differ.
    Unrelated,
}

/// The stereo-aware identity of a molecule: its constitution and the
/// stereochemistry laid over it, canonicalized as one.
///
/// Two molecules a coloring makes isomorphic share a form exactly when they are
/// the same stereoisomer — pseudo-asymmetric centers and all — so a `StereoForm`
/// compares, orders, and hashes as a portable identity, and answers the questions
/// the constitution alone cannot: whether the molecule is chiral, and how it
/// relates to another.
///
/// Obtain via [`form`].
#[derive(Clone, Debug)]
pub struct StereoForm<VK, EK> {
    constitution: Canonical<VK, EK>,
    // The stereo-aware canonical form: the identity molecules of the same
    // stereoisomer share, and the basis of `Eq`, `Ord`, and `Hash`. The
    // constitution and mirror are functions of it, kept only for `constitution`,
    // `is_chiral`, and `relate`.
    configured: Configured<VK, EK>,
    mirror: Configured<VK, EK>,
}

impl<VK, EK> StereoForm<VK, EK> {
    /// Returns `true` if the molecule is chiral: not superimposable on its mirror
    /// image. A molecule with no stereocenters, or a meso form whose centers
    /// cancel, is achiral.
    pub fn is_chiral(&self) -> bool
    where
        VK: PartialEq,
        EK: PartialEq,
    {
        self.configured != self.mirror
    }

    /// The canonical constitution underlying this form — its identity blind to
    /// stereochemistry, which enantiomers share.
    pub fn constitution(&self) -> &Canonical<VK, EK> {
        &self.constitution
    }

    /// How this molecule relates to `other`: the same stereoisomer, its
    /// enantiomer, a diastereomer, or — if their constitutions differ — unrelated.
    pub fn relate(&self, other: &StereoForm<VK, EK>) -> StereoRelationship
    where
        VK: PartialEq,
        EK: PartialEq,
    {
        if self.constitution != other.constitution {
            StereoRelationship::Unrelated
        } else if self.configured == other.configured {
            StereoRelationship::Identical
        } else if self.configured == other.mirror {
            StereoRelationship::Enantiomers
        } else {
            StereoRelationship::Diastereomers
        }
    }
}

impl<VK: PartialEq, EK: PartialEq> PartialEq for StereoForm<VK, EK> {
    fn eq(&self, other: &Self) -> bool {
        self.configured == other.configured
    }
}

impl<VK: Eq, EK: Eq> Eq for StereoForm<VK, EK> {}

impl<VK: Ord, EK: Ord> PartialOrd for StereoForm<VK, EK> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<VK: Ord, EK: Ord> Ord for StereoForm<VK, EK> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.configured.cmp(&other.configured)
    }
}

impl<VK: Hash, EK: Hash> Hash for StereoForm<VK, EK> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.configured.hash(state);
    }
}

/// The stereo-aware identity of a molecule under the caller's coloring.
///
/// Canonically labels the constitution (as [`canonicalize`]) and, over it, the
/// declared configurations — refined to a fixpoint so stereochemistry that breaks
/// a constitutional symmetry is resolved. `site_key` and `bond_key` define
/// identity exactly as they do for the constitution.
///
/// # Complexity
///
/// A bounded number of canonical labelings, each O(V · (V + E) · log V) time and
/// O(V + E) space, over the molecule's `V` sites and `E` bonds — one refinement
/// per constitutional symmetry stereochemistry breaks, at most `V`.
pub fn form<M, VK, EK>(
    mol: &M,
    site_key: impl Fn(SiteId) -> VK,
    bond_key: impl Fn(BondId) -> EK,
) -> StereoForm<VK, EK>
where
    M: HasStereoConfigurations,
    VK: Ord,
    EK: Ord,
{
    StereoForm {
        constitution: canonicalize(mol, &site_key, &bond_key),
        configured: refined(mol, &site_key, &bond_key, false),
        mirror: refined(mol, &site_key, &bond_key, true),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::hash_map::DefaultHasher;

    use vita_core::HasSites;

    use crate::CoordinationGeometry::*;
    use crate::{
        CoordinationGeometry, HasBonds, StereoConfiguration, StereoKind, StereoLocus,
        StereogenicGeometry,
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

    impl HasStereoConfigurations for Mol {
        fn stereo_configurations(&self) -> impl Iterator<Item = StereoConfiguration> + '_ {
            self.configs.iter().cloned()
        }
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
            configs,
        }
    }

    fn stereo_form(mol: &Mol) -> StereoForm<u32, u32> {
        form(mol, |site| mol.color(site), |_| 0u32)
    }

    fn hash_of(form: &StereoForm<u32, u32>) -> u64 {
        let mut hasher = DefaultHasher::new();
        form.hash(&mut hasher);
        hasher.finish()
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            colors: m.colors.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
            configs: m.configs.iter().rev().cloned().collect(),
        }
    }

    fn single(order: [u32; 4]) -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
            vec![config(1, order)],
        )
    }

    fn recolored(order: [u32; 4]) -> Mol {
        mol(
            &[(1, 0), (2, 5), (3, 6), (4, 7), (5, 8)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
            vec![config(1, order)],
        )
    }

    fn tartaric(first: [u32; 4], second: [u32; 4]) -> Mol {
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
            vec![config(1, first), config(2, second)],
        )
    }

    fn trihydroxyglutaric(first: [u32; 4], middle: [u32; 4], third: [u32; 4]) -> Mol {
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
            vec![config(1, first), config(2, middle), config(3, third)],
        )
    }

    fn cis_trans(bond: u32, order: [u32; 4]) -> StereoConfiguration {
        StereoConfiguration::new(StereoLocus::Bond(b(bond)), StereoKind::Bond, order.map(s))
            .unwrap()
    }

    fn butene(order: [u32; 4]) -> Mol {
        mol(
            &[(1, 1), (2, 0), (3, 0), (4, 1), (5, 2), (6, 2)],
            &[(1, 2, 3), (2, 1, 2), (3, 3, 4), (4, 2, 5), (5, 3, 6)],
            vec![cis_trans(1, order)],
        )
    }

    #[test]
    fn a_molecule_without_configurations_is_achiral() {
        let plain = mol(&[(1, 0), (2, 1)], &[(1, 1, 2)], Vec::new());
        assert!(!stereo_form(&plain).is_chiral());
    }

    #[test]
    fn a_single_stereocenter_is_chiral() {
        assert!(stereo_form(&single([2, 3, 4, 5])).is_chiral());
    }

    #[test]
    fn a_meso_form_is_achiral() {
        assert!(!stereo_form(&tartaric([2, 3, 4, 5], [6, 1, 7, 8])).is_chiral());
    }

    #[test]
    fn a_homochiral_form_is_chiral() {
        assert!(stereo_form(&tartaric([2, 3, 4, 5], [1, 6, 7, 8])).is_chiral());
    }

    #[test]
    fn a_double_bond_geometry_is_achiral() {
        assert!(!stereo_form(&butene([1, 5, 4, 6])).is_chiral());
    }

    #[test]
    fn a_pseudo_asymmetric_center_distinguishes_its_diastereomers() {
        let r = stereo_form(&trihydroxyglutaric(
            [2, 4, 5, 10],
            [1, 6, 7, 3],
            [8, 2, 9, 11],
        ));
        let t = stereo_form(&trihydroxyglutaric(
            [2, 4, 5, 10],
            [3, 6, 7, 1],
            [8, 2, 9, 11],
        ));
        assert_ne!(r, t);
    }

    #[test]
    fn a_dormant_pseudo_center_does_not_split_its_diastereomers() {
        let a = stereo_form(&trihydroxyglutaric(
            [2, 4, 5, 10],
            [1, 6, 7, 3],
            [2, 8, 9, 11],
        ));
        let c = stereo_form(&trihydroxyglutaric(
            [2, 4, 5, 10],
            [3, 6, 7, 1],
            [2, 8, 9, 11],
        ));
        assert_eq!(a, c);
    }

    #[test]
    fn a_form_relates_to_itself_as_identical() {
        let form = stereo_form(&single([2, 3, 4, 5]));
        assert_eq!(form.relate(&form), StereoRelationship::Identical);
    }

    #[test]
    fn opposite_configurations_are_enantiomers() {
        let right = stereo_form(&single([2, 3, 4, 5]));
        let left = stereo_form(&single([3, 2, 4, 5]));
        assert_eq!(right.relate(&left), StereoRelationship::Enantiomers);
    }

    #[test]
    fn the_homochiral_and_meso_forms_are_diastereomers() {
        let homochiral = stereo_form(&tartaric([2, 3, 4, 5], [1, 6, 7, 8]));
        let meso = stereo_form(&tartaric([2, 3, 4, 5], [6, 1, 7, 8]));
        assert_eq!(homochiral.relate(&meso), StereoRelationship::Diastereomers);
    }

    #[test]
    fn the_double_bond_isomers_are_diastereomers() {
        let cis = stereo_form(&butene([1, 5, 4, 6]));
        let trans = stereo_form(&butene([1, 5, 6, 4]));
        assert_eq!(cis.relate(&trans), StereoRelationship::Diastereomers);
    }

    #[test]
    fn molecules_of_different_constitution_are_unrelated() {
        let a = stereo_form(&single([2, 3, 4, 5]));
        let c = stereo_form(&recolored([2, 3, 4, 5]));
        assert_eq!(a.relate(&c), StereoRelationship::Unrelated);
    }

    #[test]
    fn enantiomers_share_a_constitution() {
        let right = stereo_form(&single([2, 3, 4, 5]));
        let left = stereo_form(&single([3, 2, 4, 5]));
        assert_eq!(right.constitution(), left.constitution());
    }

    #[test]
    fn the_same_stereoisomer_has_equal_forms() {
        assert_eq!(
            stereo_form(&single([2, 3, 4, 5])),
            stereo_form(&single([2, 3, 4, 5])),
        );
    }

    #[test]
    fn enantiomers_have_distinct_forms() {
        assert_ne!(
            stereo_form(&single([2, 3, 4, 5])),
            stereo_form(&single([3, 2, 4, 5])),
        );
    }

    #[test]
    fn equal_forms_hash_equally() {
        assert_eq!(
            hash_of(&stereo_form(&single([2, 3, 4, 5]))),
            hash_of(&stereo_form(&single([2, 3, 4, 5]))),
        );
    }

    #[test]
    fn ordering_agrees_with_equality() {
        let right = stereo_form(&single([2, 3, 4, 5]));
        let left = stereo_form(&single([3, 2, 4, 5]));
        assert_eq!(right.cmp(&right), Ordering::Equal);
        assert_ne!(right.cmp(&left), Ordering::Equal);
    }

    #[test]
    fn the_form_is_independent_of_input_order() {
        let molecule = tartaric([2, 3, 4, 5], [1, 6, 7, 8]);
        assert_eq!(stereo_form(&molecule), stereo_form(&reversed(&molecule)));
    }

    #[test]
    fn the_pseudo_asymmetric_form_is_independent_of_input_order() {
        let molecule = trihydroxyglutaric([2, 4, 5, 10], [1, 6, 7, 3], [8, 2, 9, 11]);
        assert_eq!(stereo_form(&molecule), stereo_form(&reversed(&molecule)));
    }
}

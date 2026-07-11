use std::cmp::Ordering;

use vita_core::SiteId;

use super::{Token, geometry};
use crate::algorithm::canonical::{Canonical, canonicalize};
use crate::{BondId, HasStereoConfigurations, StereoKind, StereoLocus};

/// A configuration reduced to an arrangement-invariant fact: the orbit its locus
/// anchors, its kind, and its class token.
type Entry = (OrbitKey, StereoKind, Token);

/// The canonical ranks of the orbits a locus anchors — the key under which
/// symmetry-equivalent configurations are pooled, so which of two equivalent centres
/// takes which rank cannot change the identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct OrbitKey(u64);

/// The stereochemistry of a molecule, reduced against a canonical labeling.
///
/// Each configuration becomes an [`Entry`] keyed on the orbit it anchors; pooling
/// symmetry-equivalent loci into a sorted set makes the layer invariant to the atom
/// order and to the arrangement of an orbit's shared ranks — which is what lets a
/// meso form register as its own mirror image.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct StereoLayer {
    entries: Vec<Entry>,
}

impl StereoLayer {
    /// Reduces `mol`'s configurations against `canon`.
    fn of<M, VK, EK>(mol: &M, canon: &Canonical<VK, EK>) -> Self
    where
        M: HasStereoConfigurations,
    {
        let rank = |site: SiteId| {
            canon
                .rank(site)
                .expect("a configuration's site is in the molecule")
        };
        let orbit_rank = |site: SiteId| {
            canon
                .orbit(site)
                .expect("a configuration's site is in the molecule")
                .iter()
                .map(rank)
                .min()
                .expect("an orbit is non-empty") as u64
        };
        let key = |locus: StereoLocus| match locus {
            StereoLocus::Site(site) | StereoLocus::Axis(site) => OrbitKey(orbit_rank(site)),
            StereoLocus::Bond(bond) => {
                let (a, b) = mol.bond_endpoints(bond);
                let (lo, hi) = (
                    orbit_rank(a).min(orbit_rank(b)),
                    orbit_rank(a).max(orbit_rank(b)),
                );
                OrbitKey(lo << 32 | hi)
            }
        };

        let mut entries: Vec<Entry> = mol
            .stereo_configurations()
            .map(|config| {
                let class = geometry(config.kind()).token(config.neighbors(), rank);
                (key(config.locus()), config.kind(), class)
            })
            .collect();
        entries.sort_unstable();
        StereoLayer { entries }
    }

    /// The layer of the mirror image: every token reflected, then re-sorted.
    fn mirror(&self) -> Self {
        let mut entries: Vec<Entry> = self
            .entries
            .iter()
            .map(|&(key, kind, class)| (key, kind, geometry(kind).mirror(class)))
            .collect();
        entries.sort_unstable();
        StereoLayer { entries }
    }

    /// Returns `true` if the molecule is chiral.
    fn is_chiral(&self) -> bool {
        self.entries != self.mirror().entries
    }
}

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

/// The stereo-aware identity of a molecule: its canonical constitution together with
/// the stereochemistry laid over it.
///
/// Two molecules a coloring makes isomorphic share a form exactly when they are the
/// same stereoisomer, so a `StereoForm` compares, orders, and hashes as a portable
/// identity — key a registry by it — and answers the questions the constitution
/// alone cannot: whether the molecule is chiral, and how it relates to another.
///
/// Obtain via [`form`].
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereoForm<VK, EK> {
    canonical: Canonical<VK, EK>,
    layer: StereoLayer,
}

impl<VK: Ord, EK: Ord> PartialOrd for StereoForm<VK, EK> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<VK: Ord, EK: Ord> Ord for StereoForm<VK, EK> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.canonical
            .cmp(&other.canonical)
            .then_with(|| self.layer.cmp(&other.layer))
    }
}

impl<VK, EK> StereoForm<VK, EK> {
    /// Returns `true` if the molecule is chiral: not superimposable on its mirror
    /// image. A molecule with no stereocentres, or a meso form whose centres cancel,
    /// is achiral.
    pub fn is_chiral(&self) -> bool {
        self.layer.is_chiral()
    }

    /// The canonical constitution underlying this form.
    pub fn constitution(&self) -> &Canonical<VK, EK> {
        &self.canonical
    }
}

impl<VK: PartialEq, EK: PartialEq> StereoForm<VK, EK> {
    /// How this molecule relates to `other`: the same stereoisomer, its enantiomer, a
    /// diastereomer, or — if their constitutions differ — unrelated.
    pub fn relate(&self, other: &StereoForm<VK, EK>) -> StereoRelationship {
        if self.canonical != other.canonical {
            StereoRelationship::Unrelated
        } else if self.layer == other.layer {
            StereoRelationship::Identical
        } else if self.layer == other.layer.mirror() {
            StereoRelationship::Enantiomers
        } else {
            StereoRelationship::Diastereomers
        }
    }
}

/// The stereo-aware identity of a molecule under the caller's coloring.
///
/// Canonically labels the constitution (as [`canonicalize`]) and reduces the
/// declared configurations against it. `site_key` and `bond_key` define identity
/// exactly as they do for the constitution.
///
/// # Complexity
///
/// One canonical labeling, O(V · (V + E) · log V) time and O(V + E) space, over the
/// molecule's `V` sites and `E` bonds.
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
    let canonical = canonicalize(mol, site_key, bond_key);
    let layer = StereoLayer::of(mol, &canonical);
    StereoForm { canonical, layer }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    use vita_core::HasSites;

    use crate::{HasBonds, StereoConfiguration};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
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

    fn plain() -> Mol {
        mol(
            &[(1, 0), (2, 1), (3, 2), (4, 3), (5, 4)],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
            Vec::new(),
        )
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

    fn paired(first: [u32; 4], second: [u32; 4]) -> Mol {
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

    fn distinct_pair(first: [u32; 4], second: [u32; 4]) -> Mol {
        mol(
            &[
                (1, 0),
                (2, 1),
                (3, 2),
                (4, 3),
                (5, 0),
                (6, 4),
                (7, 5),
                (8, 6),
            ],
            &[
                (1, 1, 2),
                (2, 1, 3),
                (3, 1, 4),
                (4, 1, 5),
                (5, 5, 6),
                (6, 5, 7),
                (7, 5, 8),
            ],
            vec![config(1, first), config(5, second)],
        )
    }

    #[test]
    fn a_molecule_without_configurations_is_achiral() {
        assert!(!stereo_form(&plain()).is_chiral());
    }

    #[test]
    fn a_single_stereocenter_is_chiral() {
        assert!(stereo_form(&single([2, 3, 4, 5])).is_chiral());
    }

    #[test]
    fn a_meso_form_is_achiral() {
        assert!(!stereo_form(&paired([2, 3, 4, 5], [6, 1, 7, 8])).is_chiral());
    }

    #[test]
    fn a_form_of_matched_centers_is_chiral() {
        assert!(stereo_form(&paired([2, 3, 4, 5], [1, 6, 7, 8])).is_chiral());
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
    fn a_pair_flipped_at_one_center_are_diastereomers() {
        let a = stereo_form(&distinct_pair([2, 3, 4, 5], [1, 6, 7, 8]));
        let b = stereo_form(&distinct_pair([2, 3, 4, 5], [6, 1, 7, 8]));
        assert_eq!(a.relate(&b), StereoRelationship::Diastereomers);
    }

    #[test]
    fn a_pair_flipped_at_both_centers_are_enantiomers() {
        let a = stereo_form(&distinct_pair([2, 3, 4, 5], [1, 6, 7, 8]));
        let b = stereo_form(&distinct_pair([3, 2, 4, 5], [6, 1, 7, 8]));
        assert_eq!(a.relate(&b), StereoRelationship::Enantiomers);
    }

    #[test]
    fn molecules_of_different_constitution_are_unrelated() {
        let a = stereo_form(&single([2, 3, 4, 5]));
        let b = stereo_form(&recolored([2, 3, 4, 5]));
        assert_eq!(a.relate(&b), StereoRelationship::Unrelated);
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
        let molecule = single([2, 3, 4, 5]);
        assert_eq!(stereo_form(&molecule), stereo_form(&reversed(&molecule)));
    }
}

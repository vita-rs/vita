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

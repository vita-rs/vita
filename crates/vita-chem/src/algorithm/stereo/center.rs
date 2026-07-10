use vita_core::SiteId;

use super::{frame, geometry};
use crate::algorithm::canonical::{Canonical, canonicalize};
use crate::{BondId, HasBondOrders, StereoKind, StereoLocus};

/// The stereogenic units of a molecule: the loci whose substituents symmetry cannot
/// interchange, so more than one configuration is realisable.
///
/// Membership is topological — a function of the graph and the caller's coloring,
/// independent of any configuration a source may declare. A locus survives when the
/// substituent classes its geometry admits realise more than one distinct
/// arrangement, so repeated substituents are weighed correctly: an octahedral MA₄B₂
/// centre is stereogenic (cis and trans), an MA₅B one is not.
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

    /// Iterates the stereogenic loci ascending: sites, then bonds, then axes.
    pub fn iter(&self) -> impl Iterator<Item = StereoLocus> + '_ {
        self.loci.iter().copied()
    }
}

/// The stereogenic loci of a molecule.
///
/// `candidate` is the caller's: it reports which loci could bear stereochemistry and
/// of which [`StereoKind`] — a centre's coordination geometry — which the library
/// cannot know. The library then keeps only those a symmetry does *not* neutralise:
/// `site_key` and `bond_key` colour the graph exactly as for [`canonicalize`], and a
/// candidate survives when, once its frame is individualised, its substituent classes
/// admit more than one arrangement under the geometry's group.
///
/// The substituents are taken from the bonds present; a stereocentre bearing an
/// implicit hydrogen must give it explicitly, or fold the coordination into
/// `site_key`. A bond's or axis's substituents hang off the termini of its rigid
/// double-bond chain, so a plain double bond and a cumulene resolve alike.
///
/// # Complexity
///
/// O(L · V · (V + E) · log V) time and O(V + E) space, over the molecule's `V` sites
/// and `E` bonds, for `L` candidate loci — one [`canonicalize`] each.
pub fn stereocenters<M, VK, EK>(
    mol: &M,
    site_key: impl Fn(SiteId) -> VK,
    bond_key: impl Fn(BondId) -> EK,
    candidate: impl Fn(StereoLocus) -> Option<StereoKind>,
) -> Stereocenters
where
    M: HasBondOrders,
    VK: Ord,
    EK: Ord,
{
    let surviving = |locus| {
        let kind = candidate(locus)?;
        stereogenic(mol, locus, kind, &site_key, &bond_key).then_some(locus)
    };
    let sites = mol
        .sites()
        .filter_map(|site| surviving(StereoLocus::Site(site)));
    let axes = mol
        .sites()
        .filter_map(|site| surviving(StereoLocus::Axis(site)));
    let bonds = mol
        .bonds()
        .filter_map(|bond| surviving(StereoLocus::Bond(bond)));

    let mut loci: Vec<StereoLocus> = sites.chain(axes).chain(bonds).collect();
    loci.sort_unstable();
    Stereocenters { loci }
}

/// Whether the frame `locus` names realises more than one configuration of `kind`
/// once individualised.
fn stereogenic<M, VK, EK>(
    mol: &M,
    locus: StereoLocus,
    kind: StereoKind,
    site_key: &impl Fn(SiteId) -> VK,
    bond_key: &impl Fn(BondId) -> EK,
) -> bool
where
    M: HasBondOrders,
    VK: Ord,
    EK: Ord,
{
    if !locus.anchors(kind) {
        return false;
    }
    let Some(frame) = frame(mol, locus) else {
        return false;
    };
    if frame.substituents.len() != kind.slot_count() {
        return false;
    }
    let canon = canonicalize(
        mol,
        |other| (frame_mark(other, &frame.anchors), site_key(other)),
        bond_key,
    );
    geometry(kind).configuration_count(&classes(&canon, &frame.substituents)) > 1
}

/// The symmetry class of each site: the least canonical rank in its orbit, so
/// interchangeable substituents share a class.
fn classes<VK, EK>(canon: &Canonical<VK, EK>, sites: &[SiteId]) -> Vec<usize> {
    sites
        .iter()
        .map(|&site| {
            canon
                .orbit(site)
                .expect("a substituent is in the molecule")
                .iter()
                .map(|member| {
                    canon
                        .rank(member)
                        .expect("an orbit member is in the molecule")
                })
                .min()
                .expect("an orbit is non-empty")
        })
        .collect()
}

/// A distinct individualisation mark for each anchor of a frame — its position in
/// `anchors` plus one, or zero for a site outside it — so a canonical labeling holds
/// the frame fixed and tells its members apart. A site marks itself, an edge or axis
/// its two termini.
fn frame_mark(site: SiteId, anchors: &[SiteId]) -> u8 {
    anchors
        .iter()
        .position(|&anchor| anchor == site)
        .map_or(0, |index| index as u8 + 1)
}

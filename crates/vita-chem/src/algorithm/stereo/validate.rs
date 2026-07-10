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
/// As [`stereocenters`]: O(L · V · (V + E) · log V) time and O(V + E) space, over the
/// molecule's `V` sites and `E` bonds, for `L` candidate loci.
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

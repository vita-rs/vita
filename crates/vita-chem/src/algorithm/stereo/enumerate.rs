use vita_core::SiteId;

use super::{
    StereoConfigurations, StereoForm, candidate_loci, form, realisable, refined_key, settle,
};
use crate::{BondId, HasBondOrders, StereoConfiguration, StereoKind, StereoLocus};

/// The stereoisomers of a molecule's constitution — every assignment of a
/// configuration to its stereogenic units that no symmetry equates.
///
/// Each isomer is the set of [`StereoConfigurations`] it bears, over one shared
/// constitution; bind it to the molecule with [`StereoConfigurations::bind`] to read
/// it as a whole. Enantiomers appear as distinct isomers, a meso form once, and a
/// pseudo-asymmetric centre — stereogenic only for a particular pairing of its
/// neighbours' handedness — is enumerated exactly where it is realised.
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
    /// molecule with no stereocentres is its own sole isomer.
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
/// stereochemistry and of which [`StereoKind`], the latter two colour the graph. The
/// stereogenic units are found under the configurations assigned so far and one is
/// given each of its realisable configurations in turn, recursively — so a
/// pseudo-asymmetric centre surfaces once the neighbours it depends on are fixed.
/// Assignments the molecule's symmetry equates, and the two halves of a meso form,
/// collapse to one isomer by their [`StereoForm`].
///
/// The molecule's own declared configurations, if any, are ignored: the enumeration
/// speaks for the constitution alone. Every locus that can ever be stereogenic is
/// assigned in some isomer, so the assigned loci, unioned over the isomers, are the
/// constitution's complete potential stereocentres — pseudo-asymmetric ones included.
///
/// # Complexity
///
/// O(N · (V + L) · V · (V + E) · log V) time and O(N · k + V + E) space, over the
/// molecule's `V` sites and `E` bonds, for `L` candidate loci and `N` assignments of
/// configurations to the `k` stereocentres — up to the product of their configuration
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
/// each of its realisable configurations in turn, or records the assignment once none
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
        let configs = realisable(&view, locus, candidate(locus)?, &key, bond_key);
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

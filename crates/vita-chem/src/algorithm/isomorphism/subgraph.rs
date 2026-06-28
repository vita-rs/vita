use std::collections::HashMap;

use vita_core::SiteId;

use crate::utils::embeddings;
use crate::{BondId, HasBonds};

/// A subgraph match: an injection of a pattern's sites onto a target's under
/// which every pattern bond has a counterpart in the target.
///
/// Maps each pattern site to the target site it stands for.
///
/// Obtain via [`matches`](matches()).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mapping {
    pairs: Vec<(SiteId, SiteId)>,
}

impl Mapping {
    /// Number of matched sites.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Returns `true` if the match maps no sites.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Returns the target site `pattern_site` is matched to.
    ///
    /// Returns `None` if `pattern_site` is not in the pattern.
    pub fn get(&self, pattern_site: SiteId) -> Option<SiteId> {
        self.pairs
            .binary_search_by_key(&pattern_site, |&(pattern, _)| pattern)
            .ok()
            .map(|i| self.pairs[i].1)
    }

    /// Iterates the matched `(pattern site, target site)` pairs, ordered by
    /// pattern site.
    pub fn iter(&self) -> impl Iterator<Item = (SiteId, SiteId)> + '_ {
        self.pairs.iter().copied()
    }
}

/// Finds every way `pattern` occurs as a subgraph of `target`.
///
/// The match is the caller's to define: `site_match` decides which pattern site
/// may stand for which target site, and `bond_match` which pattern bond for which
/// target bond — pass element and bond-order equality to match by constitution,
/// loosen either to honour a query's own rules, and the library imposes no
/// default. Each match injects every pattern site onto a distinct target site so
/// that every pattern bond meets a `bond_match`-ing target bond; the target may
/// bear further bonds among the matched sites, so the match is a subgraph, not an
/// induced one.
///
/// # Complexity
///
/// Exponential in the worst case, as subgraph isomorphism is NP-complete;
/// near-linear in practice for the connected patterns chemistry poses.
pub fn matches<P, T>(
    pattern: &P,
    target: &T,
    site_match: impl Fn(SiteId, SiteId) -> bool,
    bond_match: impl Fn(BondId, BondId) -> bool,
) -> impl Iterator<Item = Mapping>
where
    P: HasBonds,
    T: HasBonds,
{
    let Indexed {
        sites: pattern_sites,
        bonds: pattern_bonds,
        adjacency: pattern_adjacency,
    } = index(pattern);
    let Indexed {
        sites: target_sites,
        bonds: target_bonds,
        adjacency: target_adjacency,
    } = index(target);
    let compat_pattern_sites = pattern_sites.clone();
    let compat_target_sites = target_sites.clone();

    embeddings(
        &pattern_adjacency,
        target_adjacency,
        move |p, t| site_match(compat_pattern_sites[p], compat_target_sites[t]),
        move |pe, te| bond_match(pattern_bonds[pe], target_bonds[te]),
    )
    .map(move |mapping| {
        let mut pairs: Vec<(SiteId, SiteId)> = mapping
            .iter()
            .enumerate()
            .map(|(p, &t)| (pattern_sites[p], target_sites[t]))
            .collect();
        pairs.sort_unstable();
        Mapping { pairs }
    })
}

/// A `(neighbour, bond index)` adjacency list over the sites `0..site_count`.
type Adjacency = Vec<Vec<(usize, usize)>>;

/// A molecule's sites and bonds in order, with the adjacency over their indices —
/// the form the matching engine consumes.
struct Indexed {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
    adjacency: Adjacency,
}

/// Indexes a molecule for the matching engine.
fn index<M: HasBonds>(mol: &M) -> Indexed {
    let sites: Vec<SiteId> = mol.sites().collect();
    let position: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let bonds: Vec<BondId> = mol.bonds().collect();
    let mut adjacency: Adjacency = vec![Vec::new(); sites.len()];
    for (edge, &bond) in bonds.iter().enumerate() {
        let (a, b) = mol.bond_endpoints(bond);
        adjacency[position[&a]].push((position[&b], edge));
        adjacency[position[&b]].push((position[&a], edge));
    }
    Indexed {
        sites,
        bonds,
        adjacency,
    }
}

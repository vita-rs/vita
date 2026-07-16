use vita_core::SiteId;

use super::{Fingerprint, combine};
use crate::HasBonds;
use crate::algorithm::topology::path::distances;

/// The atom-pair fingerprint of a molecule (Carhart).
///
/// Each feature hashes an unordered pair of atoms with the topological distance
/// between them — `(distance, lesser seed, greater seed)` — using only the atoms'
/// `site_seed` colouring; bonds enter solely through the distance. A feature's
/// multiplicity is the number of pairs realising it. Pairs in different connected
/// components, or farther apart than `max_distance`, are skipped.
///
/// Independent of input order: each unordered pair contributes once, and ordering
/// the two seeds makes the code blind to which atom is which.
///
/// # Complexity
///
/// O(V · (V + E) + V² · log V) time and O(V²) space, over the molecule's `V`
/// sites and `E` bonds — the all-pairs distances, then a scan of the `V²` pairs.
pub fn pair<M: HasBonds>(
    mol: &M,
    site_seed: impl Fn(SiteId) -> u64,
    max_distance: usize,
) -> Fingerprint {
    let matrix = distances(mol);
    let sites: Vec<SiteId> = mol.sites().collect();
    let seeds: Vec<u64> = sites.iter().map(|&site| site_seed(site)).collect();

    let mut codes: Vec<u64> = Vec::new();
    for i in 0..sites.len() {
        for j in (i + 1)..sites.len() {
            if let Some(distance) = matrix.get(sites[i], sites[j])
                && (1..=max_distance).contains(&distance)
            {
                let lesser = seeds[i].min(seeds[j]);
                let greater = seeds[i].max(seeds[j]);
                codes.push(combine([distance as u64, lesser, greater]));
            }
        }
    }
    Fingerprint::from_codes(codes)
}

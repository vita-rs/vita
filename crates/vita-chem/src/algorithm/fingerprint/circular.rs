use std::hash::{Hash, Hasher};

use vita_core::SiteId;

use super::{Fingerprint, Indexed, combine, index};
use crate::algorithm::utils::{BitSet, FxHashMap, FxHasher};
use crate::{BondId, HasBonds, StereoConfiguration, StereoDescriptor};

/// The circular fingerprint of a molecule — its ECFP.
///
/// Each feature hashes a site's environment: at radius `0` the site's own
/// `site_seed`, and at each further radius the previous feature refined by the
/// sorted `(bond_seed, neighbour feature)` pairs around it, in the manner of
/// Morgan / Weisfeiler–Leman refinement. A feature's multiplicity is the number of
/// environments realising it — every site at radius `0`, and beyond it every
/// distinct reach, so a radius that grows no further contributes once, not again
/// per site.
///
/// The chemistry is the caller's throughout. `site_seed` and `bond_seed` decide
/// which sites and bonds count as the same; `stereo` supplies the configuration
/// anchored at a site, if any, which refines that site's feature by the arrangement
/// its neighbours realise — so a centre and its mirror, alike in the graph, part
/// here. Pass `|_| None` to fingerprint the constitution alone. A configuration
/// whose neighbours the refinement cannot yet tell apart contributes nothing until a
/// radius wide enough to separate them, which is what it means for the arrangement
/// to be stereogenic at all.
///
/// The result is independent of the order the molecule presents its sites and
/// bonds. Feature codes derive from the colouring and the *sorted* neighbourhood,
/// never from indices; a configuration enters as its [`StereoDescriptor`], taken
/// against that same ranking;
/// environments that recur are pooled by the bonds they cover, keeping the
/// lowest-radius, lowest-code representative; and a site's own radius-0 code
/// bypasses that pooling so its count is preserved.
///
/// # Complexity
///
/// O(R · (E · log Δ + V · log V + (V + E) · E / w) + V · R · log(V · R)) time and
/// O(V · R · E / w) space, over the molecule's `V` sites and `E` bonds, for radius
/// `R`, maximum degree `Δ`, and word width `w = 64` — each round ranks the sites,
/// sorts every neighbourhood, and grows every site's covered-bond set, whose bits
/// both pool the environments and dominate the cost; the surviving codes are counted
/// at the end.
pub fn circular<M: HasBonds>(
    mol: &M,
    site_seed: impl Fn(SiteId) -> u64,
    bond_seed: impl Fn(BondId) -> u64,
    stereo: impl Fn(SiteId) -> Option<StereoConfiguration>,
    radius: usize,
) -> Fingerprint {
    let Indexed {
        sites,
        bonds,
        adjacency,
    } = index(mol);
    let n = sites.len();
    let edges = bonds.len();

    let position: FxHashMap<SiteId, usize> = sites
        .iter()
        .enumerate()
        .map(|(i, &site)| (site, i))
        .collect();
    let configurations: Vec<Option<StereoConfiguration>> =
        sites.iter().map(|&site| stereo(site)).collect();
    let bond_seeds: Vec<u64> = bonds.iter().map(|&bond| bond_seed(bond)).collect();
    let mut codes: Vec<u64> = sites
        .iter()
        .map(|&site| combine([site_seed(site)]))
        .collect();

    // Radius 0: each site's own code, one feature per site, never pooled.
    let radius_zero = codes.clone();

    let mut bondsets = vec![BitSet::zeros(edges); n];
    // Environments keyed by the bonds they cover, each keeping the (radius, code)
    // it is least at. The choice is decided by those values alone, never by the
    // key's bit pattern, so the survivors are the same however the molecule was
    // indexed — and the codes are counted, hence sorted, on the way out, so the
    // map's own order never reaches the result.
    let mut pooled: FxHashMap<BitSet, (usize, u64)> = FxHashMap::default();

    for r in 1..=radius {
        // The round reads the previous codes and bondsets and writes fresh ones, so
        // no site sees a neighbour already advanced this round — the guarantee that
        // the codes depend on the graph, not the iteration order.
        let mut next_codes = codes.clone();
        let mut next_bondsets = bondsets.clone();
        let ranking = ranking(&codes);
        for i in 0..n {
            let mut neighbours: Vec<(u64, u64)> = adjacency
                .neighbors(i)
                .iter()
                .map(|&(edge, j)| (bond_seeds[edge], codes[j]))
                .collect();
            neighbours.sort_unstable();

            let mut words = Vec::with_capacity(3 + 2 * neighbours.len());
            words.push(r as u64);
            words.push(codes[i]);
            for (edge_seed, neighbour_code) in neighbours {
                words.push(edge_seed);
                words.push(neighbour_code);
            }
            if let Some(configuration) = &configurations[i] {
                let descriptor = configuration.descriptor(|site| ranking[position[&site]]);
                words.push(hashed(&descriptor));
            }
            next_codes[i] = combine(words);

            let mut ball = bondsets[i].clone();
            for &(edge, j) in adjacency.neighbors(i) {
                ball.set(edge);
                ball |= &bondsets[j];
            }
            // An empty ball is an isolated site's degenerate environment, already
            // captured by its radius-0 code.
            if !ball.is_zero() {
                let candidate = (r, next_codes[i]);
                match pooled.get_mut(&ball) {
                    Some(least) if candidate < *least => *least = candidate,
                    Some(_) => {}
                    None => {
                        pooled.insert(ball.clone(), candidate);
                    }
                }
            }
            next_bondsets[i] = ball;
        }
        codes = next_codes;
        bondsets = next_bondsets;
    }

    Fingerprint::from_codes(
        radius_zero
            .into_iter()
            .chain(pooled.into_values().map(|(_, code)| code)),
    )
}

/// Ranks each site by its code, densely and by magnitude.
///
/// A configuration ranks its neighbours to fix which arrangement they realise, and
/// takes their order alone. Ranking the codes densely rather than passing them as
/// they are keeps that order faithful on a target whose `usize` is narrower than a
/// code.
fn ranking(codes: &[u64]) -> Vec<usize> {
    let mut distinct: Vec<u64> = codes.to_vec();
    distinct.sort_unstable();
    distinct.dedup();
    codes
        .iter()
        .map(|code| {
            distinct
                .binary_search(code)
                .expect("the distinct codes hold every code")
        })
        .collect()
}

/// Folds a descriptor into one word, through the same seedless hash the codes use.
fn hashed(descriptor: &StereoDescriptor) -> u64 {
    let mut hasher = FxHasher::default();
    descriptor.hash(&mut hasher);
    hasher.finish()
}

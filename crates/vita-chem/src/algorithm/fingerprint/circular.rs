use std::hash::{Hash, Hasher};

use vita_core::SiteId;

use super::{Fingerprint, Indexed, combine, index};
use crate::algorithm::utils::{BitSet, FxHashMap, FxHasher};
use crate::{BondId, HasBonds, StereoConfiguration, StereoDescriptor};

/// The circular fingerprint of a molecule — its ECFP.
///
/// Each feature hashes a site's environment: at radius `0` the site's own
/// `site_seed`, and at each further radius the previous feature refined by the
/// sorted `(bond_seed, neighbor feature)` pairs around it, in the manner of
/// Morgan / Weisfeiler–Leman refinement. A feature's multiplicity is the number of
/// environments realising it — every site at radius `0`, and beyond it every
/// distinct reach, so a radius that grows no further contributes once, not again
/// per site.
///
/// The chemistry is the caller's throughout. `site_seed` and `bond_seed` decide
/// which sites and bonds count as the same; `stereo` supplies the configuration
/// anchored at a site, if any, which refines that site's feature by the arrangement
/// its neighbors realize — so a center and its mirror, alike in the graph, part
/// here. Pass `|_| None` to fingerprint the constitution alone. A configuration
/// whose neighbors the refinement cannot yet tell apart contributes nothing until a
/// radius wide enough to separate them, which is what it means for the arrangement
/// to be stereogenic at all.
///
/// The result is independent of the order the molecule presents its sites and
/// bonds. Feature codes derive from the coloring and the *sorted* neighborhood,
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
/// sorts every neighborhood, and grows every site's covered-bond set, whose bits
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
        // no site sees a neighbor already advanced this round — the guarantee that
        // the codes depend on the graph, not the iteration order.
        let mut next_codes = codes.clone();
        let mut next_bondsets = bondsets.clone();
        let ranking = ranking(&codes);
        for i in 0..n {
            let mut neighbors: Vec<(u64, u64)> = adjacency
                .neighbors(i)
                .iter()
                .map(|&(edge, j)| (bond_seeds[edge], codes[j]))
                .collect();
            neighbors.sort_unstable();

            let mut words = Vec::with_capacity(3 + 2 * neighbors.len());
            words.push(r as u64);
            words.push(codes[i]);
            for (edge_seed, neighbor_code) in neighbors {
                words.push(edge_seed);
                words.push(neighbor_code);
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
/// A configuration ranks its neighbors to fix which arrangement they realize, and
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

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::CoordinationGeometry::*;
    use crate::fingerprint::FeatureVector;
    use crate::{CoordinationGeometry, StereoKind, StereoLocus, StereogenicGeometry};

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
        colors: Vec<u64>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
    }

    impl Mol {
        fn color(&self, site: SiteId) -> u64 {
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

    fn mol(atoms: &[(u32, u64)], bonds: &[(u32, u32, u32)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, _)| s(id)).collect(),
            colors: atoms.iter().map(|&(_, color)| color).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c)| (s(a), s(c))).collect(),
        }
    }

    fn fingerprint(m: &Mol, radius: usize) -> Fingerprint {
        circular(m, |site| m.color(site), |_| 0, |_| None, radius)
    }

    fn configured(m: &Mol, site: u32, order: [u32; 4], radius: usize) -> Fingerprint {
        let configuration = config(site, order);
        circular(
            m,
            |x| m.color(x),
            |_| 0,
            |x| (x == s(site)).then(|| configuration.clone()),
            radius,
        )
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            colors: m.colors.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().map(|&(a, c)| (c, a)).collect(),
        }
    }

    fn dumbbell() -> Mol {
        mol(&[(1, 1), (2, 1)], &[(1, 1, 2)])
    }

    fn chain() -> Mol {
        mol(&[(1, 1), (2, 1), (3, 1)], &[(1, 1, 2), (2, 2, 3)])
    }

    fn triangle() -> Mol {
        mol(
            &[(1, 1), (2, 1), (3, 1)],
            &[(1, 1, 2), (2, 2, 3), (3, 1, 3)],
        )
    }

    fn star(neighbors: [u64; 4]) -> Mol {
        mol(
            &[
                (1, 0),
                (2, neighbors[0]),
                (3, neighbors[1]),
                (4, neighbors[2]),
                (5, neighbors[3]),
            ],
            &[(1, 1, 2), (2, 1, 3), (3, 1, 4), (4, 1, 5)],
        )
    }

    fn tailed_star() -> Mol {
        mol(
            &[
                (1, 0),
                (2, 7),
                (3, 7),
                (4, 7),
                (5, 7),
                (6, 11),
                (7, 12),
                (8, 13),
                (9, 14),
            ],
            &[
                (1, 1, 2),
                (2, 1, 3),
                (3, 1, 4),
                (4, 1, 5),
                (5, 2, 6),
                (6, 3, 7),
                (7, 4, 8),
                (8, 5, 9),
            ],
        )
    }

    #[test]
    fn an_empty_molecule_has_an_empty_fingerprint() {
        assert!(fingerprint(&mol(&[], &[]), 2).is_empty());
    }

    #[test]
    fn a_radius_zero_fingerprint_counts_the_sites_by_color() {
        let print = fingerprint(&mol(&[(1, 1), (2, 2), (3, 1)], &[(1, 1, 2), (2, 2, 3)]), 0);
        assert_eq!(print.len(), 2);
        assert_eq!(print.cardinality(), 3);
    }

    #[test]
    fn sites_of_one_color_share_a_feature() {
        let print = fingerprint(&mol(&[(1, 7), (2, 7)], &[]), 0);
        assert_eq!(print.len(), 1);
        assert_eq!(print.cardinality(), 2);
    }

    #[test]
    fn sites_of_different_colors_have_distinct_features() {
        assert_eq!(fingerprint(&mol(&[(1, 7), (2, 8)], &[]), 0).len(), 2);
    }

    #[test]
    fn a_wider_radius_separates_distinct_environments() {
        let print = fingerprint(&chain(), 1);
        assert_eq!(print.len(), 3);
        assert_eq!(print.cardinality(), 6);
    }

    #[test]
    fn bond_colors_enter_the_environment() {
        let m = dumbbell();
        let plain = circular(&m, |x| m.color(x), |_| 0, |_| None, 1);
        let recolored = circular(&m, |x| m.color(x), |_| 1, |_| None, 1);
        assert_ne!(plain, recolored);
    }

    #[test]
    fn equivalent_environments_with_distinct_reaches_each_count() {
        let print = fingerprint(&triangle(), 1);
        assert_eq!(print.len(), 2);
        assert_eq!(print.cardinality(), 6);
    }

    #[test]
    fn an_isolated_site_contributes_only_its_radius_zero_feature() {
        let m = mol(&[(1, 7)], &[]);
        assert_eq!(fingerprint(&m, 3), fingerprint(&m, 0));
    }

    #[test]
    fn sites_covering_the_same_bonds_pool_into_one_environment() {
        assert_eq!(fingerprint(&dumbbell(), 1).cardinality(), 3);
    }

    #[test]
    fn a_radius_beyond_the_molecules_reach_adds_nothing() {
        assert_eq!(fingerprint(&dumbbell(), 2), fingerprint(&dumbbell(), 1));
    }

    #[test]
    fn a_configuration_is_invisible_at_radius_zero() {
        let m = star([2, 3, 4, 5]);
        assert_eq!(configured(&m, 1, [2, 3, 4, 5], 0), fingerprint(&m, 0));
    }

    #[test]
    fn a_configuration_refines_its_centers_feature() {
        let m = star([2, 3, 4, 5]);
        assert_ne!(configured(&m, 1, [2, 3, 4, 5], 1), fingerprint(&m, 1));
    }

    #[test]
    fn mirror_configurations_yield_distinct_fingerprints() {
        let m = star([2, 3, 4, 5]);
        assert_ne!(
            configured(&m, 1, [2, 3, 4, 5], 1),
            configured(&m, 1, [3, 2, 4, 5], 1),
        );
    }

    #[test]
    fn a_configuration_among_indistinguishable_neighbors_separates_nothing() {
        let m = star([7, 7, 7, 7]);
        assert_eq!(
            configured(&m, 1, [2, 3, 4, 5], 1),
            configured(&m, 1, [3, 2, 4, 5], 1),
        );
    }

    #[test]
    fn a_configuration_parts_mirror_forms_once_the_radius_resolves_its_neighbors() {
        let m = tailed_star();
        assert_eq!(
            configured(&m, 1, [2, 3, 4, 5], 1),
            configured(&m, 1, [3, 2, 4, 5], 1),
        );
        assert_ne!(
            configured(&m, 1, [2, 3, 4, 5], 2),
            configured(&m, 1, [3, 2, 4, 5], 2),
        );
    }

    #[test]
    fn the_fingerprint_is_independent_of_input_order() {
        let m = tailed_star();
        assert_eq!(
            configured(&m, 1, [2, 3, 4, 5], 2),
            configured(&reversed(&m), 1, [2, 3, 4, 5], 2),
        );
    }
}

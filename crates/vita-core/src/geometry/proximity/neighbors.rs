use super::{Arrangement, arrangement, periodic_arrangement};
use crate::units::length::{Length, LengthUnit};
use crate::{HasLattice, HasPositions, Scalar, SiteId};

/// The sites within reach of one another, as a symmetric relation over those that
/// have a neighbor.
///
/// Nearness here is a membership question: the relation answers whether it
/// [`contains`](Self::contains) a pair, which sites lie within reach
/// [`of_site`](Self::of_site), and how many pairs there are in all; the pairs
/// themselves come from [`iter`](Self::iter), in ascending order. It keeps the
/// relation alone — how far apart a pair stands is read from the arrangement that
/// gathered it.
///
/// Obtain via [`neighbors`], [`periodic_neighbors`], or [`Arrangement::neighbors`]
/// where an arrangement is already at hand.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Neighbors {
    sites: Vec<SiteId>,
    offsets: Vec<usize>,
    targets: Vec<SiteId>,
}

impl Neighbors {
    /// The number of pairs.
    pub fn len(&self) -> usize {
        self.targets.len() / 2
    }

    /// Returns `true` if there are no pairs.
    pub fn is_empty(&self) -> bool {
        self.targets.is_empty()
    }

    /// Returns `true` if `a` and `b` are within reach of one another.
    ///
    /// Symmetric, and `false` for a site paired with itself or absent altogether.
    pub fn contains(&self, a: SiteId, b: SiteId) -> bool {
        self.reach(a).binary_search(&b).is_ok()
    }

    /// Iterates the sites within reach of `site`, in ascending order.
    ///
    /// Empty for a site with no neighbor, and for one absent altogether.
    pub fn of_site(&self, site: SiteId) -> impl Iterator<Item = SiteId> + '_ {
        self.reach(site).iter().copied()
    }

    /// Iterates the pairs, each once as `(a, b)` with `a` before `b`, in ascending
    /// order.
    pub fn iter(&self) -> impl Iterator<Item = (SiteId, SiteId)> + '_ {
        self.sites.iter().enumerate().flat_map(|(index, &site)| {
            let row = &self.targets[self.offsets[index]..self.offsets[index + 1]];
            row[row.partition_point(|target| *target < site)..]
                .iter()
                .map(move |&target| (site, target))
        })
    }

    /// The sites within reach of `site`, ascending, empty if it has none.
    fn reach(&self, site: SiteId) -> &[SiteId] {
        match self.sites.binary_search(&site) {
            Ok(index) => &self.targets[self.offsets[index]..self.offsets[index + 1]],
            Err(_) => &[],
        }
    }
}

impl<V: Scalar> Arrangement<V> {
    /// The pairs within reach of one another, indexed by site.
    ///
    /// The same pairs [`pairs`](Self::pairs) yields, gathered into a relation that can
    /// be asked about a site instead of read through once. They are counted before
    /// they are placed, so nothing proportional to their number is ever held beyond
    /// the relation itself.
    pub fn neighbors(&self) -> Neighbors {
        let arranged = self.sites();
        let count = arranged.len();
        let mut order: Vec<u32> = (0..count as u32).collect();
        order.sort_unstable_by_key(|&slot| arranged[slot as usize]);
        let mut ranks = vec![0u32; count];
        for (rank, &slot) in order.iter().enumerate() {
            ranks[slot as usize] = rank as u32;
        }
        let mut sites: Vec<SiteId> = order.iter().map(|&slot| arranged[slot as usize]).collect();
        let mut offsets = vec![0usize; count + 1];
        self.couples().for_each(|(near, far, _)| {
            offsets[ranks[near] as usize + 1] += 1;
            offsets[ranks[far] as usize + 1] += 1;
        });
        for index in 1..offsets.len() {
            offsets[index] += offsets[index - 1];
        }
        let mut targets = (sites.first()).map_or_else(Vec::new, |&site| vec![site; offsets[count]]);
        self.couples().for_each(|(near, far, _)| {
            let (near, far) = (ranks[near] as usize, ranks[far] as usize);
            targets[offsets[near]] = sites[far];
            offsets[near] += 1;
            targets[offsets[far]] = sites[near];
            offsets[far] += 1;
        });
        offsets.rotate_right(1);
        offsets[0] = 0;
        let mut kept = 0;
        for index in 0..count {
            if offsets[index] < offsets[index + 1] {
                targets[offsets[index]..offsets[index + 1]].sort_unstable();
                sites[kept] = sites[index];
                offsets[kept] = offsets[index];
                kept += 1;
            }
        }
        offsets[kept] = targets.len();
        sites.truncate(kept);
        offsets.truncate(kept + 1);
        sites.shrink_to_fit();
        offsets.shrink_to_fit();
        Neighbors {
            sites,
            offsets,
            targets,
        }
    }
}

/// The pairs of sites no further apart than `cutoff`.
///
/// A negative cutoff puts nothing within reach; a cutoff of zero leaves only the
/// sites that coincide, at the cost of comparing every pair.
///
/// # Panics
///
/// Panics if any site of [`sites`](crate::HasSites::sites) is absent from
/// [`positions`](HasPositions::positions).
pub fn neighbors<S, V, U>(system: &S, cutoff: Length<V, U>) -> Neighbors
where
    S: HasPositions<V>,
    V: Scalar,
    U: LengthUnit,
{
    arrangement(system, cutoff).neighbors()
}

/// The pairs of sites no further apart than `cutoff` on the torus the lattice
/// defines, each pair counted once however many images bring them together.
///
/// A negative cutoff puts nothing within reach; a cutoff of zero leaves only the
/// sites that coincide, at the cost of comparing every pair. A site is never its own
/// neighbor, whatever the reach.
///
/// # Panics
///
/// Panics if any site of [`sites`](crate::HasSites::sites) is absent from
/// [`positions`](HasPositions::positions).
pub fn periodic_neighbors<S, V, U>(system: &S, cutoff: Length<V, U>) -> Neighbors
where
    S: HasPositions<V> + HasLattice<V>,
    V: Scalar,
    U: LengthUnit,
{
    periodic_arrangement(system, cutoff).neighbors()
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::fixture::{configuration, s};
    use crate::geometry::proximity::fixture::{Descending, chain, cube, reach, strewn, sweep};
    use crate::units::length::{Angstrom, Nanometer};

    #[test]
    fn an_empty_system_has_no_neighbors() {
        assert!(neighbors(&configuration(&[]), reach(1.0)).is_empty());
    }

    #[test]
    fn the_pairs_are_counted_once_each() {
        assert_eq!(neighbors(&chain(), reach(1.5)).len(), 3);
    }

    #[test]
    fn contains_answers_whether_a_pair_is_within_reach() {
        let found = neighbors(&chain(), reach(1.5));
        assert!(found.contains(s(1), s(2)) && !found.contains(s(1), s(3)));
    }

    #[test]
    fn contains_is_symmetric() {
        let found = neighbors(&chain(), reach(1.5));
        assert_eq!(found.contains(s(1), s(2)), found.contains(s(2), s(1)));
    }

    #[test]
    fn a_site_is_not_its_own_neighbor() {
        let found = periodic_neighbors(&cube(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), reach(6.0));
        assert!(!found.contains(s(1), s(1)));
    }

    #[test]
    fn of_site_lists_the_sites_within_reach() {
        let found = neighbors(&chain(), reach(1.5));
        assert_eq!(found.of_site(s(2)).collect::<Vec<_>>(), vec![s(1), s(3)]);
    }

    #[test]
    fn iter_lists_every_pair_once() {
        let found = neighbors(&chain(), reach(1.5));
        assert_eq!(
            found.iter().collect::<Vec<_>>(),
            vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4))]
        );
    }

    #[test]
    fn an_absent_site_reaches_nothing() {
        assert_eq!(neighbors(&chain(), reach(1.5)).of_site(s(9)).count(), 0);
    }

    #[test]
    fn a_site_out_of_reach_of_the_rest_reaches_nothing() {
        let system = configuration(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [9.0, 0.0, 0.0]]);
        assert_eq!(neighbors(&system, reach(1.5)).of_site(s(3)).count(), 0);
    }

    #[test]
    fn every_pair_is_reached_from_both_of_its_sites() {
        let found = neighbors(&configuration(&strewn()), reach(3.0));
        assert!(
            found
                .iter()
                .all(|(a, b)| found.contains(a, b) && found.of_site(b).any(|site| site == a))
        );
    }

    #[test]
    fn the_cutoff_is_taken_in_the_requested_unit() {
        let found = neighbors(&chain(), Length::<f64, Nanometer>::new(0.15));
        assert_eq!(found.len(), 3);
    }

    #[test]
    fn neighbors_are_independent_of_the_order_the_sites_arrive_in() {
        let found = neighbors(&Descending(chain()), reach(1.5));
        assert_eq!(found, neighbors(&chain(), reach(1.5)));
    }

    #[test]
    fn the_relation_holds_exactly_the_pairs_the_arrangement_yields() {
        let system = configuration(&strewn());
        let mut yielded: Vec<(SiteId, SiteId)> = (arrangement(&system, reach(3.0))
            .pairs::<Angstrom>())
        .map(|(a, b, _)| (a, b))
        .collect();
        yielded.sort_unstable();
        let found = neighbors(&system, reach(3.0));
        assert_eq!(found.iter().collect::<Vec<_>>(), yielded);
    }

    #[test]
    fn the_wrapped_relation_holds_exactly_the_pairs_a_full_sweep_would() {
        let system = cube(&strewn());
        let found = periodic_neighbors(&system, reach(3.0));
        assert_eq!(found.iter().collect::<Vec<_>>(), sweep(&system, reach(3.0)));
    }
}

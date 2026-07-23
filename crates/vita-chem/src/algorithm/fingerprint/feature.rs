use std::cmp::Ordering;
use std::ops::{Add, AddAssign};

use super::{FeatureVector, FoldedFingerprint};
use crate::algorithm::utils::SortedMap;

/// A molecule's structural fingerprint: the multiset of hashed local
/// substructure features whose overlap measures similarity.
///
/// Each feature is an opaque `u64`, a stable hash of one colored local
/// substructure, paired with the number of times it occurs. The codes are
/// reproducible across runs and machines and comparable across molecules by
/// construction, so a fingerprint persists and screens against others built the
/// same way — the same generator, parameters, and coloring, a caller contract
/// the type does not police.
///
/// The empty fingerprint is [`Default`]; [`AddAssign`] pools features — uniting
/// two molecules' or a whole dataset's — and [`fold`](Self::fold) projects onto a
/// fixed width. Similarity is measured through [`FeatureVector`]:
/// [`tanimoto`](super::FeatureVector::tanimoto),
/// [`tversky`](super::FeatureVector::tversky), or
/// [`cosine`](super::FeatureVector::cosine).
///
/// Obtain via a generator ([`circular`](super::circular), [`pair`](super::pair),
/// [`path`](super::path)) or, for a custom feature enumeration,
/// [`from_codes`](Self::from_codes).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fingerprint(SortedMap<u64, usize>);

impl Fingerprint {
    /// Counts `codes` into a fingerprint: each distinct code becomes one feature,
    /// its multiplicity the number of times the code appears.
    ///
    /// The open construction seam — a caller enumerates a molecule's substructures,
    /// hashes each to a code, and collects them here. The generators feed it.
    ///
    /// # Complexity
    ///
    /// O(K · log K) time and O(K) space, over the `K` codes.
    pub fn from_codes(codes: impl IntoIterator<Item = u64>) -> Self {
        let mut codes: Vec<u64> = codes.into_iter().collect();
        codes.sort_unstable();
        Fingerprint(SortedMap::from_pairs(
            codes.chunk_by(|a, b| a == b).map(|run| (run[0], run.len())),
        ))
    }

    /// Builds a fingerprint from `(feature, count)` pairs, summing repeated
    /// features and dropping those with zero count.
    ///
    /// Reads back what [`iter`](Self::iter) writes out: the pair carries a
    /// fingerprint through storage and returns it unchanged, no serialization
    /// format required.
    ///
    /// # Complexity
    ///
    /// O(K · log K) time and O(K) space, over the `K` pairs.
    pub fn from_counts(counts: impl IntoIterator<Item = (u64, usize)>) -> Self {
        let mut pairs: Vec<(u64, usize)> = counts.into_iter().filter(|&(_, n)| n > 0).collect();
        pairs.sort_unstable_by_key(|&(code, _)| code);
        Fingerprint(SortedMap::from_pairs(
            pairs
                .chunk_by(|a, b| a.0 == b.0)
                .map(|run| (run[0].0, run.iter().map(|&(_, n)| n).sum())),
        ))
    }

    /// The number of distinct features.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if the fingerprint has no features.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// The number of times `feature` occurs, or `0` if it is absent.
    pub fn count(&self, feature: u64) -> usize {
        self.0.get(&feature).copied().unwrap_or(0)
    }

    /// Iterates `(feature, count)` pairs in ascending feature order.
    pub fn iter(&self) -> impl Iterator<Item = (u64, usize)> + '_ {
        self.0.iter().map(|(&feature, &count)| (feature, count))
    }

    /// Folds the features into a `width`-bit vector, hashing each to bit
    /// `feature % width`.
    ///
    /// A lossy projection onto a fixed width: the counts are dropped, and features
    /// sharing a bit become indistinguishable, in return for a dense form quick to
    /// compare and compact to store. Wider folds collide less. A `width` of zero
    /// yields an empty vector.
    ///
    /// # Complexity
    ///
    /// O(A) time and O(width) space, over the `A` distinct features.
    pub fn fold(&self, width: usize) -> FoldedFingerprint {
        FoldedFingerprint::new(self.0.iter().map(|(&feature, _)| feature), width)
    }

    /// Sums `combine(count_here, count_there)` over the features present in both.
    fn matched(&self, other: &Fingerprint, combine: impl Fn(usize, usize) -> usize) -> usize {
        let mut total = 0;
        let (mut here, mut there) = (self.0.iter(), other.0.iter());
        let (mut a, mut b) = (here.next(), there.next());
        while let (Some((code_a, count_a)), Some((code_b, count_b))) = (a, b) {
            match code_a.cmp(code_b) {
                Ordering::Less => a = here.next(),
                Ordering::Greater => b = there.next(),
                Ordering::Equal => {
                    total += combine(*count_a, *count_b);
                    a = here.next();
                    b = there.next();
                }
            }
        }
        total
    }
}

impl Default for Fingerprint {
    fn default() -> Self {
        Self::from_codes([])
    }
}

impl FeatureVector for Fingerprint {
    fn cardinality(&self) -> usize {
        self.0.iter().map(|(_, &count)| count).sum()
    }

    fn intersection(&self, other: &Self) -> usize {
        self.matched(other, |here, there| here.min(there))
    }

    fn dot(&self, other: &Self) -> usize {
        self.matched(other, |here, there| here * there)
    }
}

impl AddAssign<&Fingerprint> for Fingerprint {
    fn add_assign(&mut self, other: &Fingerprint) {
        *self = Fingerprint::from_counts(self.iter().chain(other.iter()));
    }
}

impl Add<&Fingerprint> for &Fingerprint {
    type Output = Fingerprint;

    fn add(self, other: &Fingerprint) -> Fingerprint {
        Fingerprint::from_counts(self.iter().chain(other.iter()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_fingerprint_is_empty() {
        assert!(Fingerprint::default().is_empty());
    }

    #[test]
    fn an_empty_code_sequence_yields_the_default_fingerprint() {
        assert_eq!(Fingerprint::from_codes([]), Fingerprint::default());
    }

    #[test]
    fn an_empty_count_sequence_yields_the_default_fingerprint() {
        assert_eq!(Fingerprint::from_counts([]), Fingerprint::default());
    }

    #[test]
    fn an_empty_fingerprint_has_no_features() {
        assert_eq!(Fingerprint::default().len(), 0);
    }

    #[test]
    fn an_empty_fingerprint_has_zero_cardinality() {
        assert_eq!(Fingerprint::default().cardinality(), 0);
    }

    #[test]
    fn a_single_code_becomes_one_feature() {
        assert_eq!(Fingerprint::from_codes([7]).len(), 1);
    }

    #[test]
    fn a_single_code_is_counted_once() {
        assert_eq!(Fingerprint::from_codes([7]).count(7), 1);
    }

    #[test]
    fn a_fingerprint_with_a_feature_is_not_empty() {
        assert!(!Fingerprint::from_codes([7]).is_empty());
    }

    #[test]
    fn repeated_codes_accumulate_into_one_feature() {
        let print = Fingerprint::from_codes([7, 7, 7]);
        assert_eq!(print.len(), 1);
        assert_eq!(print.count(7), 3);
    }

    #[test]
    fn distinct_codes_become_distinct_features() {
        let print = Fingerprint::from_codes([3, 7]);
        assert_eq!(print.len(), 2);
        assert_eq!(print.count(3), 1);
        assert_eq!(print.count(7), 1);
    }

    #[test]
    fn from_counts_sums_repeated_features() {
        let print = Fingerprint::from_counts([(7, 2), (3, 1), (7, 3)]);
        assert_eq!(print.len(), 2);
        assert_eq!(print.count(7), 5);
    }

    #[test]
    fn iter_yields_features_in_ascending_order_with_their_counts() {
        let print = Fingerprint::from_codes([9, 1, 9, 4]);
        let pairs: Vec<(u64, usize)> = print.iter().collect();
        assert_eq!(pairs, vec![(1, 1), (4, 1), (9, 2)]);
    }

    #[test]
    fn a_fingerprint_round_trips_through_its_counts() {
        let print = Fingerprint::from_codes([9, 1, 9, 4]);
        assert_eq!(Fingerprint::from_counts(print.iter()), print);
    }

    #[test]
    fn cardinality_sums_the_counts() {
        assert_eq!(Fingerprint::from_codes([1, 1, 2]).cardinality(), 3);
    }

    #[test]
    fn intersection_takes_the_lesser_count_of_each_shared_feature() {
        let a = Fingerprint::from_codes([1, 1, 2, 5]);
        let b = Fingerprint::from_codes([1, 2, 2, 7]);
        assert_eq!(a.intersection(&b), 2);
    }

    #[test]
    fn dot_multiplies_the_counts_of_shared_features() {
        let a = Fingerprint::from_codes([1, 1, 2, 5]);
        let b = Fingerprint::from_codes([1, 2, 2, 7]);
        assert_eq!(a.dot(&b), 4);
    }

    #[test]
    fn the_sum_of_fingerprints_pools_their_features() {
        let a = Fingerprint::from_codes([1, 1, 2]);
        let b = Fingerprint::from_codes([2, 3]);
        let sum = &a + &b;
        assert_eq!(sum.len(), 3);
        assert_eq!(sum.count(1), 2);
        assert_eq!(sum.count(2), 2);
        assert_eq!(sum.count(3), 1);
    }

    #[test]
    fn add_assign_agrees_with_add() {
        let a = Fingerprint::from_codes([1, 1, 2]);
        let b = Fingerprint::from_codes([2, 3]);
        let mut pooled = a.clone();
        pooled += &b;
        assert_eq!(pooled, &a + &b);
    }

    #[test]
    fn fold_maps_each_feature_to_its_bit_modulo_the_width() {
        let folded = Fingerprint::from_codes([3, 10]).fold(8);
        assert_eq!(folded.width(), 8);
        assert!(folded.contains(3));
        assert!(folded.contains(2));
    }

    #[test]
    fn fold_drops_the_counts() {
        assert_eq!(
            Fingerprint::from_codes([5, 5]).fold(8),
            Fingerprint::from_codes([5]).fold(8),
        );
    }

    #[test]
    fn an_absent_feature_has_a_count_of_zero() {
        assert_eq!(Fingerprint::from_codes([1]).count(2), 0);
    }

    #[test]
    fn zero_count_features_are_dropped() {
        let print = Fingerprint::from_counts([(1, 0), (2, 1)]);
        assert_eq!(print.len(), 1);
        assert_eq!(print.count(1), 0);
    }

    #[test]
    fn disjoint_fingerprints_share_nothing() {
        let a = Fingerprint::from_codes([1, 2]);
        let b = Fingerprint::from_codes([3, 4]);
        assert_eq!(a.intersection(&b), 0);
        assert_eq!(a.dot(&b), 0);
    }

    #[test]
    fn features_sharing_a_bit_fold_indistinguishably() {
        assert_eq!(
            Fingerprint::from_codes([1, 9]).fold(8),
            Fingerprint::from_codes([9]).fold(8),
        );
    }

    #[test]
    fn a_zero_width_fold_is_empty() {
        let folded = Fingerprint::from_codes([1]).fold(0);
        assert_eq!(folded.width(), 0);
        assert_eq!(folded.cardinality(), 0);
    }

    #[test]
    fn from_codes_and_from_counts_agree_on_the_same_multiset() {
        assert_eq!(
            Fingerprint::from_codes([7, 7, 3]),
            Fingerprint::from_counts([(3, 1), (7, 2)]),
        );
    }

    #[test]
    fn fingerprints_differing_only_in_counts_are_not_equal() {
        assert_ne!(
            Fingerprint::from_codes([1]),
            Fingerprint::from_codes([1, 1]),
        );
    }

    #[test]
    fn the_fingerprint_is_independent_of_input_order() {
        assert_eq!(
            Fingerprint::from_codes([3, 1, 2, 1]),
            Fingerprint::from_codes([1, 1, 2, 3]),
        );
        assert_eq!(
            Fingerprint::from_counts([(3, 1), (1, 2), (2, 1)]),
            Fingerprint::from_counts([(1, 2), (2, 1), (3, 1)]),
        );
    }
}

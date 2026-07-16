use std::cmp::Ordering;
use std::ops::{Add, AddAssign};

use super::{FeatureVector, FoldedFingerprint};
use crate::algorithm::utils::SortedMap;

/// A molecule's structural fingerprint: the multiset of hashed local
/// substructure features whose overlap measures similarity.
///
/// Each feature is an opaque `u64`, a stable hash of one coloured local
/// substructure, paired with the number of times it occurs. The codes are
/// reproducible across runs and machines and comparable across molecules by
/// construction, so a fingerprint persists and screens against others built the
/// same way — the same generator, parameters, and colouring, a caller contract
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
    /// fingerprint through storage and returns it unchanged, no serialisation
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

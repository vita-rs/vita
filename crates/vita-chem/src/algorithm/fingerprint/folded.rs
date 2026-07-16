use super::FeatureVector;
use crate::algorithm::utils::BitSet;

/// A fingerprint folded into a fixed-width bit vector.
///
/// [`Fingerprint::fold`](super::Fingerprint::fold) hashes every feature into one
/// of `width` bits, trading the sparse fingerprint's exact counts for a dense,
/// bounded form — quick to compare word-by-word and compact to store, at the cost
/// of the collisions folding introduces. It measures similarity through the same
/// [`FeatureVector`] metrics as the sparse form, where its
/// [`cardinality`](FeatureVector::cardinality) is the number of bits a feature
/// reached.
///
/// Folds of different widths lie in different spaces: measuring one against
/// another is undefined, and panics.
///
/// Obtain via [`Fingerprint::fold`](super::Fingerprint::fold).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FoldedFingerprint(BitSet);

impl FoldedFingerprint {
    /// Folds `codes` into `width` bits, setting bit `code % width` for each.
    pub(super) fn new(codes: impl IntoIterator<Item = u64>, width: usize) -> Self {
        let mut bits = BitSet::zeros(width);
        if width > 0 {
            let modulus = width as u64;
            for code in codes {
                bits.set((code % modulus) as usize);
            }
        }
        FoldedFingerprint(bits)
    }

    /// The width folded into — the number of bits the features were hashed onto.
    pub fn width(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if a feature folded onto bit `i`.
    pub fn contains(&self, i: usize) -> bool {
        self.0.test(i)
    }

    /// Iterates the bits features folded onto, ascending.
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.0.len()).filter(|&i| self.0.test(i))
    }
}

impl FeatureVector for FoldedFingerprint {
    fn cardinality(&self) -> usize {
        self.0.count_ones()
    }

    /// # Panics
    ///
    /// Panics if `other` was folded to a different width. The two then index
    /// different spaces, and no meet between them is defined.
    fn intersection(&self, other: &Self) -> usize {
        assert_eq!(
            self.0.len(),
            other.0.len(),
            "folded fingerprints of different widths do not compare",
        );
        self.0.count_and(&other.0)
    }

    /// # Panics
    ///
    /// Panics if `other` was folded to a different width, as
    /// [`intersection`](Self::intersection) does.
    fn dot(&self, other: &Self) -> usize {
        self.intersection(other)
    }
}

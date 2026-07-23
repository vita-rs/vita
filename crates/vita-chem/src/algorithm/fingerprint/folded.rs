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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_fold_of_no_codes_has_zero_cardinality() {
        assert_eq!(FoldedFingerprint::new([], 8).cardinality(), 0);
    }

    #[test]
    fn a_fold_of_no_codes_yields_no_bits() {
        assert!(FoldedFingerprint::new([], 8).iter().next().is_none());
    }

    #[test]
    fn width_reports_the_folded_width() {
        assert_eq!(FoldedFingerprint::new([], 8).width(), 8);
    }

    #[test]
    fn each_code_reaches_its_bit_modulo_the_width() {
        let folded = FoldedFingerprint::new([3, 10], 8);
        assert!(folded.contains(3));
        assert!(folded.contains(2));
    }

    #[test]
    fn iter_yields_the_reached_bits_in_ascending_order() {
        let folded = FoldedFingerprint::new([10, 3], 8);
        let bits: Vec<usize> = folded.iter().collect();
        assert_eq!(bits, vec![2, 3]);
    }

    #[test]
    fn cardinality_counts_the_reached_bits() {
        assert_eq!(FoldedFingerprint::new([1, 3, 10], 8).cardinality(), 3);
    }

    #[test]
    fn colliding_codes_reach_one_bit() {
        assert_eq!(FoldedFingerprint::new([2, 10], 8).cardinality(), 1);
    }

    #[test]
    fn intersection_counts_the_shared_bits() {
        let a = FoldedFingerprint::new([1, 2], 8);
        let b = FoldedFingerprint::new([2, 3], 8);
        assert_eq!(a.intersection(&b), 1);
    }

    #[test]
    fn dot_agrees_with_intersection() {
        let a = FoldedFingerprint::new([1, 2], 8);
        let b = FoldedFingerprint::new([2, 3], 8);
        assert_eq!(a.dot(&b), a.intersection(&b));
    }

    #[test]
    fn an_unreached_bit_is_not_contained() {
        assert!(!FoldedFingerprint::new([3], 8).contains(0));
    }

    #[test]
    fn folds_sharing_no_bit_have_zero_intersection() {
        let a = FoldedFingerprint::new([1], 8);
        let b = FoldedFingerprint::new([2], 8);
        assert_eq!(a.intersection(&b), 0);
    }

    #[test]
    fn a_zero_width_fold_ignores_its_codes() {
        let folded = FoldedFingerprint::new([1, 9], 0);
        assert_eq!(folded.width(), 0);
        assert_eq!(folded.cardinality(), 0);
    }

    #[test]
    #[should_panic(expected = "different widths")]
    fn intersection_across_widths_panics() {
        let _ = FoldedFingerprint::new([1], 8).intersection(&FoldedFingerprint::new([1], 16));
    }

    #[test]
    #[should_panic(expected = "different widths")]
    fn dot_across_widths_panics() {
        let _ = FoldedFingerprint::new([1], 8).dot(&FoldedFingerprint::new([1], 16));
    }

    #[test]
    fn folds_reaching_the_same_bits_are_equal() {
        assert_eq!(
            FoldedFingerprint::new([2, 10], 8),
            FoldedFingerprint::new([10], 8),
        );
    }

    #[test]
    fn folds_reaching_different_bits_are_not_equal() {
        assert_ne!(
            FoldedFingerprint::new([1], 8),
            FoldedFingerprint::new([2], 8),
        );
    }

    #[test]
    fn folds_of_different_widths_are_not_equal() {
        assert_ne!(
            FoldedFingerprint::new([], 8),
            FoldedFingerprint::new([], 16),
        );
    }

    #[test]
    fn the_fold_is_independent_of_input_order() {
        assert_eq!(
            FoldedFingerprint::new([3, 10, 1], 8),
            FoldedFingerprint::new([1, 3, 10], 8),
        );
    }
}

use vita_core::Scalar;

/// A feature vector and the similarities measured on it.
///
/// A fingerprint compares through three quantities of the feature counts it
/// carries — the [`cardinality`](Self::cardinality) of one vector, the
/// [`intersection`](Self::intersection) of two, and their [`dot`](Self::dot)
/// product — from which the metrics follow: the overlap metrics
/// [`tanimoto`](Self::tanimoto) and [`tversky`](Self::tversky) from the first two,
/// the vector metric [`cosine`](Self::cosine) from the third; a vector's squared
/// norm is its self-[`dot`](Self::dot).
///
/// Implementing the three invariants brings the metrics with it. Both the sparse
/// [`Fingerprint`](super::Fingerprint), through a merge of its sorted features,
/// and the dense [`FoldedFingerprint`](super::FoldedFingerprint), through a
/// population count, do — so the metrics measure either representation, and any
/// further feature vector joins them by implementing the same three.
pub trait FeatureVector {
    /// The number of features, counted with multiplicity — `|a|`, the sum of the
    /// counts.
    fn cardinality(&self) -> usize;

    /// The size of the meet with `other` — `|a ∧ b|`, summing the smaller of the
    /// two counts over every feature.
    fn intersection(&self, other: &Self) -> usize;

    /// The inner product with `other` — `⟨a, b⟩`, summing the product of the two
    /// counts over every feature.
    fn dot(&self, other: &Self) -> usize;

    /// The Tanimoto similarity to `other`.
    ///
    /// [`tversky`](Self::tversky) with unit weights — `|a ∧ b| / |a ∨ b|`, the
    /// multiset Jaccard index. Runs from `0` for disjoint vectors to `1` for equal
    /// ones; two empty vectors are maximally similar, `1`.
    ///
    /// # Complexity
    ///
    /// O(A + B) time and O(1) space, over the two vectors' `A` and `B` distinct
    /// features — O(width / w) for folded vectors of word width `w = 64`.
    fn tanimoto<S: Scalar>(&self, other: &Self) -> S {
        self.tversky(other, S::ONE, S::ONE)
    }

    /// The Tversky similarity to `other`, weighting each side's surplus.
    ///
    /// `|a ∧ b| / (|a ∧ b| + alpha · |a ∖ b| + beta · |b ∖ a|)`. Symmetric weights
    /// (`alpha == beta`) measure mutual resemblance — unit weights are
    /// [`tanimoto`](Self::tanimoto), a half each the Dice coefficient; asymmetric
    /// weights measure containment — `a.tversky(b, 1, 0)` is the fraction of `a`
    /// found in `b`, a substructure screen. Two empty vectors are maximally
    /// similar, `1`.
    ///
    /// # Complexity
    ///
    /// O(A + B) time and O(1) space, over the two vectors' `A` and `B` distinct
    /// features — O(width / w) for folded vectors of word width `w = 64`.
    fn tversky<S: Scalar>(&self, other: &Self, alpha: S, beta: S) -> S {
        if self.cardinality() == 0 && other.cardinality() == 0 {
            return S::ONE;
        }
        let shared = self.intersection(other);
        let intersection = S::from_f64(shared as f64);
        let self_surplus = S::from_f64((self.cardinality() - shared) as f64);
        let other_surplus = S::from_f64((other.cardinality() - shared) as f64);
        let denominator = intersection + alpha * self_surplus + beta * other_surplus;
        if denominator == S::ZERO {
            S::ZERO
        } else {
            intersection / denominator
        }
    }

    /// The cosine similarity to `other`.
    ///
    /// `⟨a, b⟩ / (‖a‖ · ‖b‖)`, the cosine of the angle the vectors subtend. Unlike
    /// the overlap-based [`tanimoto`](Self::tanimoto), it weighs counts as vector
    /// components, so the two diverge whenever features repeat. Runs from `0` to
    /// `1`; two empty vectors are maximally similar, `1`.
    ///
    /// # Complexity
    ///
    /// O(A + B) time and O(1) space, over the two vectors' `A` and `B` distinct
    /// features — O(width / w) for folded vectors of word width `w = 64`.
    fn cosine<S: Scalar>(&self, other: &Self) -> S {
        if self.cardinality() == 0 && other.cardinality() == 0 {
            return S::ONE;
        }
        let numerator = S::from_f64(self.dot(other) as f64);
        let norm =
            S::from_f64(self.dot(self) as f64).sqrt() * S::from_f64(other.dot(other) as f64).sqrt();
        if norm == S::ZERO {
            S::ZERO
        } else {
            numerator / norm
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Counts(Vec<usize>);

    impl FeatureVector for Counts {
        fn cardinality(&self) -> usize {
            self.0.iter().sum()
        }

        fn intersection(&self, other: &Self) -> usize {
            self.0.iter().zip(&other.0).map(|(&a, &b)| a.min(b)).sum()
        }

        fn dot(&self, other: &Self) -> usize {
            self.0.iter().zip(&other.0).map(|(&a, &b)| a * b).sum()
        }
    }

    fn counts(values: &[usize]) -> Counts {
        Counts(values.to_vec())
    }

    #[test]
    fn tanimoto_of_two_empty_vectors_is_one() {
        assert_eq!(counts(&[]).tanimoto::<f64>(&counts(&[])), 1.0);
    }

    #[test]
    fn tversky_of_two_empty_vectors_is_one() {
        assert_eq!(counts(&[]).tversky(&counts(&[]), 2.0, 3.0), 1.0);
    }

    #[test]
    fn cosine_of_two_empty_vectors_is_one() {
        assert_eq!(counts(&[]).cosine::<f64>(&counts(&[])), 1.0);
    }

    #[test]
    fn tanimoto_is_the_shared_count_over_the_union() {
        let a = counts(&[2, 1, 0]);
        let b = counts(&[1, 1, 1]);
        assert_eq!(a.tanimoto::<f64>(&b), 0.5);
    }

    #[test]
    fn tversky_weights_each_sides_surplus() {
        let a = counts(&[3, 1]);
        let b = counts(&[1, 1]);
        assert_eq!(a.tversky(&b, 1.0, 0.0), 0.5);
        assert_eq!(a.tversky(&b, 0.0, 1.0), 1.0);
    }

    #[test]
    fn half_weights_give_the_dice_coefficient() {
        let a = counts(&[2, 1]);
        let b = counts(&[1, 1]);
        assert_eq!(a.tversky(&b, 0.5, 0.5), 0.8);
    }

    #[test]
    fn cosine_is_the_dot_over_the_norms() {
        let a = counts(&[3, 4]);
        let b = counts(&[4, 3]);
        assert_eq!(a.cosine::<f64>(&b), 0.96);
    }

    #[test]
    fn tanimoto_of_disjoint_vectors_is_zero() {
        let a = counts(&[1, 0]);
        let b = counts(&[0, 1]);
        assert_eq!(a.tanimoto::<f64>(&b), 0.0);
    }

    #[test]
    fn cosine_of_disjoint_vectors_is_zero() {
        let a = counts(&[1, 0]);
        let b = counts(&[0, 1]);
        assert_eq!(a.cosine::<f64>(&b), 0.0);
    }

    #[test]
    fn tanimoto_of_an_empty_and_a_nonempty_vector_is_zero() {
        assert_eq!(counts(&[]).tanimoto::<f64>(&counts(&[1])), 0.0);
    }

    #[test]
    fn cosine_of_an_empty_and_a_nonempty_vector_is_zero() {
        assert_eq!(counts(&[]).cosine::<f64>(&counts(&[1])), 0.0);
    }

    #[test]
    fn tversky_with_zero_weights_on_disjoint_vectors_is_zero() {
        let a = counts(&[1, 0]);
        let b = counts(&[0, 1]);
        assert_eq!(a.tversky(&b, 0.0, 0.0), 0.0);
    }

    #[test]
    fn cosine_and_tanimoto_diverge_on_repeated_features() {
        let a = counts(&[2]);
        let b = counts(&[1]);
        assert_eq!(a.cosine::<f64>(&b), 1.0);
        assert_eq!(a.tanimoto::<f64>(&b), 0.5);
    }

    #[test]
    fn tanimoto_of_a_vector_with_itself_is_one() {
        let a = counts(&[2, 1]);
        assert_eq!(a.tanimoto::<f64>(&a), 1.0);
    }

    #[test]
    fn cosine_of_a_vector_with_itself_is_one() {
        let a = counts(&[3, 4]);
        assert_eq!(a.cosine::<f64>(&a), 1.0);
    }

    #[test]
    fn tanimoto_is_symmetric() {
        let a = counts(&[3, 0, 1]);
        let b = counts(&[1, 2, 1]);
        assert_eq!(a.tanimoto::<f64>(&b), b.tanimoto::<f64>(&a));
    }

    #[test]
    fn cosine_is_symmetric() {
        let a = counts(&[3, 0, 1]);
        let b = counts(&[1, 2, 1]);
        assert_eq!(a.cosine::<f64>(&b), b.cosine::<f64>(&a));
    }

    #[test]
    fn tversky_with_equal_weights_is_symmetric() {
        let a = counts(&[3, 0]);
        let b = counts(&[1, 1]);
        assert_eq!(a.tversky(&b, 2.0, 2.0), b.tversky(&a, 2.0, 2.0));
    }
}

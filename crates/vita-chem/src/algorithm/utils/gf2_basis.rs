use super::BitSet;

/// An incrementally built basis of a GF(2) vector space.
///
/// Collects a maximal linearly independent subset of the [`BitSet`] vectors
/// offered to [`insert`](Self::insert), held in reduced row echelon form: every
/// stored vector owns a distinct pivot — its lowest set bit — that is clear in
/// all the others. [`reduce`](Self::reduce) clears each pivot from a vector,
/// leaving the one representative of its coset modulo the span, which is zero
/// exactly when the vector already lies in the span.
///
/// Vectors span the same `dimension` coordinates, fixed at construction.
///
/// Obtain via [`new`](Self::new).
pub struct Gf2Basis {
    rows: Vec<Option<BitSet>>,
    rank: usize,
}

impl Gf2Basis {
    /// Creates an empty basis for vectors of `dimension` coordinates.
    pub fn new(dimension: usize) -> Self {
        Gf2Basis {
            rows: vec![None; dimension],
            rank: 0,
        }
    }

    /// The number of basis vectors — the rank of the spanned subspace.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.rank
    }

    /// Returns `true` if the basis holds no vectors.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.rank == 0
    }

    /// Reduces `vector` modulo the span, returning the canonical representative
    /// of its coset.
    ///
    /// The result is zero exactly when `vector` lies in the span.
    pub fn reduce(&self, vector: &BitSet) -> BitSet {
        let mut residue = vector.clone();
        self.reduce_in_place(&mut residue);
        residue
    }

    /// Inserts `vector` if it is independent of the span, returning `true` if it
    /// was added.
    ///
    /// A vector already in the span is rejected and leaves the basis unchanged.
    pub fn insert(&mut self, mut vector: BitSet) -> bool {
        self.reduce_in_place(&mut vector);
        let Some(pivot) = vector.lowest_set() else {
            return false;
        };
        for row in self.rows.iter_mut().flatten() {
            if row.test(pivot) {
                *row ^= &vector;
            }
        }
        self.rows[pivot] = Some(vector);
        self.rank += 1;
        true
    }

    /// Clears every basis pivot from `vector` in place.
    fn reduce_in_place(&self, vector: &mut BitSet) {
        for (pivot, row) in self.rows.iter().enumerate() {
            let Some(row) = row else { continue };
            if vector.test(pivot) {
                *vector ^= row;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vector(dimension: usize, bits: &[usize]) -> BitSet {
        let mut v = BitSet::zeros(dimension);
        for &bit in bits {
            v.set(bit);
        }
        v
    }

    #[test]
    fn empty_basis_has_rank_zero() {
        assert_eq!(Gf2Basis::new(4).len(), 0);
    }

    #[test]
    fn empty_basis_is_empty() {
        assert!(Gf2Basis::new(4).is_empty());
    }

    #[test]
    fn reducing_against_an_empty_basis_returns_the_vector() {
        let v = vector(4, &[1, 3]);
        assert_eq!(Gf2Basis::new(4).reduce(&v), v);
    }

    #[test]
    fn the_zero_vector_is_never_independent() {
        let mut b = Gf2Basis::new(4);
        assert!(!b.insert(vector(4, &[])));
        assert_eq!(b.len(), 0);
    }

    #[test]
    fn inserting_an_independent_vector_returns_true() {
        let mut b = Gf2Basis::new(4);
        assert!(b.insert(vector(4, &[0])));
    }

    #[test]
    fn inserting_an_independent_vector_increases_the_rank() {
        let mut b = Gf2Basis::new(4);
        b.insert(vector(4, &[0]));
        assert_eq!(b.len(), 1);
    }

    #[test]
    fn an_inserted_vector_lies_in_the_span() {
        let mut b = Gf2Basis::new(4);
        let v = vector(4, &[0, 2]);
        b.insert(v.clone());
        assert!(b.reduce(&v).is_zero());
    }

    #[test]
    fn inserting_a_dependent_vector_returns_false() {
        let mut b = Gf2Basis::new(4);
        b.insert(vector(4, &[0]));
        assert!(!b.insert(vector(4, &[0])));
    }

    #[test]
    fn inserting_a_dependent_vector_leaves_the_rank_unchanged() {
        let mut b = Gf2Basis::new(4);
        b.insert(vector(4, &[0]));
        b.insert(vector(4, &[0]));
        assert_eq!(b.len(), 1);
    }

    #[test]
    fn a_vector_outside_the_span_reduces_to_nonzero() {
        let mut b = Gf2Basis::new(4);
        b.insert(vector(4, &[0]));
        assert!(!b.reduce(&vector(4, &[1])).is_zero());
    }

    #[test]
    fn a_linear_combination_of_basis_vectors_reduces_to_zero() {
        let mut b = Gf2Basis::new(4);
        b.insert(vector(4, &[0]));
        b.insert(vector(4, &[1]));
        assert!(b.reduce(&vector(4, &[0, 1])).is_zero());
    }

    #[test]
    fn reduce_yields_a_canonical_coset_representative() {
        let mut b = Gf2Basis::new(4);
        b.insert(vector(4, &[0, 2]));
        let combined = b.reduce(&vector(4, &[0, 1, 2]));
        assert_eq!(combined, b.reduce(&vector(4, &[1])));
        assert_eq!(combined, vector(4, &[1]));
    }

    #[test]
    fn rank_counts_only_the_independent_insertions() {
        let mut b = Gf2Basis::new(4);
        b.insert(vector(4, &[0]));
        b.insert(vector(4, &[1]));
        b.insert(vector(4, &[0, 1]));
        b.insert(vector(4, &[2]));
        assert_eq!(b.len(), 3);
    }

    #[test]
    fn a_full_rank_basis_contains_every_vector() {
        let mut b = Gf2Basis::new(2);
        b.insert(vector(2, &[0]));
        b.insert(vector(2, &[1]));
        for bits in [vec![0], vec![1], vec![0, 1]] {
            assert!(b.reduce(&vector(2, &bits)).is_zero());
        }
    }

    #[test]
    fn the_span_is_independent_of_insertion_order() {
        let vectors = [vector(4, &[0, 1]), vector(4, &[1, 2]), vector(4, &[0, 2])];
        let mut forward = Gf2Basis::new(4);
        for v in &vectors {
            forward.insert(v.clone());
        }
        let mut reversed = Gf2Basis::new(4);
        for v in vectors.iter().rev() {
            reversed.insert(v.clone());
        }
        assert_eq!(forward.len(), reversed.len());
        let probe = vector(4, &[0, 1]);
        assert_eq!(
            forward.reduce(&probe).is_zero(),
            reversed.reduce(&probe).is_zero(),
        );
    }
}

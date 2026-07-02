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
    pub fn len(&self) -> usize {
        self.rank
    }

    /// Returns `true` if the basis holds no vectors.
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

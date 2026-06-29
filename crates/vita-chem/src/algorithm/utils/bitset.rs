/// A fixed-length bit vector with GF(2) arithmetic.
///
/// Stores `n` bits packed into 64-bit words. All bit positions start at zero.
/// Supports in-place XOR (`^=`) as GF(2) addition (symmetric difference),
/// which is the fundamental row operation for Gaussian elimination over GF(2)
/// vector spaces.
///
/// Obtain via [`zeros`](Self::zeros).
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct BitSet {
    words: Vec<u64>,
    len: usize,
}

impl BitSet {
    /// Returns a bit set of `len` zeroed bits.
    pub fn zeros(len: usize) -> Self {
        BitSet {
            words: vec![0u64; len.div_ceil(64)],
            len,
        }
    }

    /// Number of bit positions.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the bit set has no positions.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Sets bit `i` to one.
    pub fn set(&mut self, i: usize) {
        self.words[i >> 6] |= 1u64 << (i & 63);
    }

    /// Flips bit `i`.
    pub fn toggle(&mut self, i: usize) {
        self.words[i >> 6] ^= 1u64 << (i & 63);
    }

    /// Returns `true` if bit `i` is one.
    pub fn test(&self, i: usize) -> bool {
        self.words[i >> 6] >> (i & 63) & 1 == 1
    }

    /// Returns `true` if every bit is zero.
    pub fn is_zero(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }

    /// Number of bits set to one.
    pub fn count_ones(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    /// Index of the lowest set bit, or `None` if all bits are zero.
    pub fn lowest_set(&self) -> Option<usize> {
        self.words
            .iter()
            .enumerate()
            .find_map(|(k, &w)| (w != 0).then(|| k * 64 + w.trailing_zeros() as usize))
    }
}

impl std::ops::BitXorAssign<&BitSet> for BitSet {
    fn bitxor_assign(&mut self, other: &BitSet) {
        for (a, b) in self.words.iter_mut().zip(&other.words) {
            *a ^= b;
        }
    }
}

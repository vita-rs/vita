/// A fixed-length bit vector with set and GF(2) algebra.
///
/// Stores `n` bits packed into 64-bit words. All bit positions start at zero.
/// In-place XOR (`^=`) is GF(2) addition — the symmetric difference, and the row
/// operation for Gaussian elimination over GF(2) vector spaces; in-place OR (`|=`)
/// is set union, and [`count_and`](Self::count_and) sizes the intersection without
/// materialising it.
///
/// Obtain via [`zeros`](Self::zeros).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
    #[allow(dead_code)]
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

    /// Number of positions set in both this bit set and `other`.
    pub fn count_and(&self, other: &BitSet) -> usize {
        self.words
            .iter()
            .zip(&other.words)
            .map(|(a, b)| (a & b).count_ones() as usize)
            .sum()
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

impl std::ops::BitOrAssign<&BitSet> for BitSet {
    fn bitor_assign(&mut self, other: &BitSet) {
        for (a, b) in self.words.iter_mut().zip(&other.words) {
            *a |= b;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zero(len: usize) -> BitSet {
        BitSet::zeros(len)
    }

    fn single(len: usize, i: usize) -> BitSet {
        let mut b = BitSet::zeros(len);
        b.set(i);
        b
    }

    #[test]
    fn empty_bitset_has_zero_length() {
        assert_eq!(zero(0).len(), 0);
    }

    #[test]
    fn empty_bitset_is_empty() {
        assert!(zero(0).is_empty());
    }

    #[test]
    fn nonempty_bitset_is_not_empty() {
        assert!(!zero(1).is_empty());
    }

    #[test]
    fn zeros_is_zero() {
        assert!(zero(8).is_zero());
    }

    #[test]
    fn zeros_spanning_multiple_words_is_zero() {
        assert!(zero(70).is_zero());
    }

    #[test]
    fn zeros_has_no_set_bits() {
        assert_eq!(zero(8).count_ones(), 0);
    }

    #[test]
    fn zeros_lowest_set_is_none() {
        assert_eq!(zero(8).lowest_set(), None);
    }

    #[test]
    fn len_matches_capacity() {
        assert_eq!(zero(0).len(), 0);
        assert_eq!(zero(1).len(), 1);
        assert_eq!(zero(64).len(), 64);
        assert_eq!(zero(100).len(), 100);
    }

    #[test]
    fn set_bit_tests_true() {
        assert!(single(8, 3).test(3));
    }

    #[test]
    fn set_bit_increments_count() {
        assert_eq!(single(8, 3).count_ones(), 1);
    }

    #[test]
    fn set_bit_is_lowest_set() {
        assert_eq!(single(8, 3).lowest_set(), Some(3));
    }

    #[test]
    fn set_is_idempotent() {
        let mut b = zero(8);
        b.set(4);
        b.set(4);
        assert_eq!(b.count_ones(), 1);
        assert!(b.test(4));
    }

    #[test]
    fn toggle_on_zero_sets_bit() {
        let mut b = zero(8);
        b.toggle(5);
        assert!(b.test(5));
    }

    #[test]
    fn unset_bit_tests_false() {
        assert!(!zero(8).test(3));
    }

    #[test]
    fn toggle_on_set_bit_clears_it() {
        let mut b = single(8, 5);
        b.toggle(5);
        assert!(!b.test(5));
    }

    #[test]
    fn toggle_on_set_bit_restores_zero() {
        let mut b = single(8, 5);
        b.toggle(5);
        assert!(b.is_zero());
    }

    #[test]
    fn first_bit_is_settable() {
        assert!(single(8, 0).test(0));
    }

    #[test]
    fn last_bit_in_first_word_is_settable() {
        assert!(single(128, 63).test(63));
    }

    #[test]
    fn first_bit_in_second_word_is_settable() {
        assert!(single(128, 64).test(64));
    }

    #[test]
    fn lowest_set_crosses_word_boundary_correctly() {
        let mut b = zero(130);
        b.set(100);
        b.set(120);
        assert_eq!(b.lowest_set(), Some(100));
        b.set(3);
        assert_eq!(b.lowest_set(), Some(3));
    }

    #[test]
    fn multiple_set_bits_count_is_correct() {
        let mut b = zero(8);
        b.set(1);
        b.set(3);
        b.set(7);
        assert_eq!(b.count_ones(), 3);
    }

    #[test]
    fn count_ones_spans_multiple_words() {
        let mut b = zero(130);
        b.set(3);
        b.set(65);
        b.set(129);
        assert_eq!(b.count_ones(), 3);
    }

    #[test]
    fn count_and_is_the_intersection_size() {
        let mut a = zero(8);
        a.set(1);
        a.set(2);
        let mut b = zero(8);
        b.set(2);
        b.set(3);
        assert_eq!(a.count_and(&b), 1);
    }

    #[test]
    fn count_and_of_disjoint_sets_is_zero() {
        let mut a = zero(8);
        a.set(1);
        let mut b = zero(8);
        b.set(2);
        assert_eq!(a.count_and(&b), 0);
    }

    #[test]
    fn count_and_with_self_is_count_ones() {
        let mut a = zero(8);
        a.set(1);
        a.set(3);
        let copy = a.clone();
        assert_eq!(a.count_and(&copy), a.count_ones());
    }

    #[test]
    fn count_and_spans_multiple_words() {
        let mut a = zero(130);
        a.set(3);
        a.set(65);
        a.set(129);
        let mut b = zero(130);
        b.set(65);
        b.set(129);
        assert_eq!(a.count_and(&b), 2);
    }

    #[test]
    fn lowest_set_returns_minimum_index() {
        let mut b = zero(8);
        b.set(7);
        b.set(2);
        b.set(5);
        assert_eq!(b.lowest_set(), Some(2));
    }

    #[test]
    fn xor_assign_is_symmetric_difference() {
        let mut a = zero(8);
        a.set(1);
        a.set(2);
        let mut b = zero(8);
        b.set(2);
        b.set(3);
        a ^= &b;
        assert!(a.test(1));
        assert!(!a.test(2));
        assert!(a.test(3));
    }

    #[test]
    fn xor_assign_with_self_yields_zero() {
        let mut a = zero(80);
        a.set(5);
        a.set(70);
        let copy = a.clone();
        a ^= &copy;
        assert!(a.is_zero());
    }

    #[test]
    fn xor_assign_is_commutative() {
        let mut a = zero(8);
        a.set(1);
        a.set(3);
        let mut b = zero(8);
        b.set(3);
        b.set(5);

        let mut a_xor_b = a.clone();
        a_xor_b ^= &b;

        let mut b_xor_a = b;
        b_xor_a ^= &a;

        assert_eq!(a_xor_b, b_xor_a);
    }

    #[test]
    fn or_assign_is_union() {
        let mut a = zero(8);
        a.set(1);
        a.set(2);
        let mut b = zero(8);
        b.set(2);
        b.set(3);
        a |= &b;
        assert!(a.test(1));
        assert!(a.test(2));
        assert!(a.test(3));
    }

    #[test]
    fn or_assign_with_self_leaves_it_unchanged() {
        let mut a = zero(80);
        a.set(5);
        a.set(70);
        let copy = a.clone();
        a |= &copy;
        assert_eq!(a, copy);
    }

    #[test]
    fn or_assign_is_commutative() {
        let mut a = zero(8);
        a.set(1);
        a.set(3);
        let mut b = zero(8);
        b.set(3);
        b.set(5);

        let mut a_or_b = a.clone();
        a_or_b |= &b;

        let mut b_or_a = b;
        b_or_a |= &a;

        assert_eq!(a_or_b, b_or_a);
    }
}

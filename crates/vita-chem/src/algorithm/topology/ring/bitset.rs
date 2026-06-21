/// A fixed-length bit vector over a molecule's bonds, used for GF(2)
/// arithmetic in the cycle space.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(super) struct Bits {
    words: Vec<u64>,
}

impl Bits {
    /// Returns a bit vector of `len` zeroed bits.
    pub(super) fn zeros(len: usize) -> Self {
        Bits {
            words: vec![0; len.div_ceil(64)],
        }
    }

    /// Sets bit `i` to one.
    pub(super) fn set(&mut self, i: usize) {
        self.words[i >> 6] |= 1 << (i & 63);
    }

    /// Flips bit `i`.
    pub(super) fn toggle(&mut self, i: usize) {
        self.words[i >> 6] ^= 1 << (i & 63);
    }

    /// Returns `true` if bit `i` is set.
    pub(super) fn test(&self, i: usize) -> bool {
        self.words[i >> 6] >> (i & 63) & 1 == 1
    }

    /// Adds `other` in place (GF(2) addition; symmetric difference).
    pub(super) fn xor(&mut self, other: &Bits) {
        for (a, b) in self.words.iter_mut().zip(&other.words) {
            *a ^= b;
        }
    }

    /// Returns the index of the lowest set bit, or `None` if all bits are zero.
    pub(super) fn lowest_set(&self) -> Option<usize> {
        self.words
            .iter()
            .enumerate()
            .find(|&(_, &w)| w != 0)
            .map(|(k, &w)| k * 64 + w.trailing_zeros() as usize)
    }

    /// Returns `true` if no bit is set.
    pub(super) fn is_zero(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }

    /// Returns the number of set bits.
    pub(super) fn count_ones(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zeros_is_zero() {
        let b = Bits::zeros(70);
        assert!(b.is_zero());
        assert_eq!(b.count_ones(), 0);
        assert_eq!(b.lowest_set(), None);
    }

    #[test]
    fn set_and_test() {
        let mut b = Bits::zeros(130);
        b.set(0);
        b.set(65);
        b.set(129);
        assert!(b.test(0));
        assert!(b.test(65));
        assert!(b.test(129));
        assert!(!b.test(1));
        assert!(!b.test(64));
        assert_eq!(b.count_ones(), 3);
    }

    #[test]
    fn toggle_flips() {
        let mut b = Bits::zeros(8);
        b.toggle(3);
        assert!(b.test(3));
        b.toggle(3);
        assert!(!b.test(3));
    }

    #[test]
    fn lowest_set_crosses_words() {
        let mut b = Bits::zeros(130);
        b.set(100);
        b.set(120);
        assert_eq!(b.lowest_set(), Some(100));
        b.set(3);
        assert_eq!(b.lowest_set(), Some(3));
    }

    #[test]
    fn xor_is_symmetric_difference() {
        let mut a = Bits::zeros(8);
        a.set(1);
        a.set(2);
        let mut b = Bits::zeros(8);
        b.set(2);
        b.set(3);
        a.xor(&b);
        assert!(a.test(1));
        assert!(!a.test(2));
        assert!(a.test(3));
    }

    #[test]
    fn xor_with_self_is_zero() {
        let mut a = Bits::zeros(80);
        a.set(5);
        a.set(70);
        let copy = a.clone();
        a.xor(&copy);
        assert!(a.is_zero());
    }

    #[test]
    fn ordering_is_lexicographic_by_lowest_bit() {
        let mut a = Bits::zeros(8);
        a.set(1);
        let mut b = Bits::zeros(8);
        b.set(2);
        assert!(a < b);
    }

    #[test]
    fn equal_and_hashable() {
        use std::collections::HashSet;
        let mut a = Bits::zeros(8);
        a.set(4);
        let mut b = Bits::zeros(8);
        b.set(4);
        assert_eq!(a, b);
        let mut seen = HashSet::new();
        seen.insert(a);
        assert!(seen.contains(&b));
    }
}

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};

/// The FxHash multiplicative constant: an odd 64-bit integer derived from the
/// golden ratio, which scatters each absorbed word across the whole hash.
const MULTIPLIER: u64 = 0x51_7c_c1_b7_27_22_0a_95;

/// A fast, deterministic hasher for the small integer keys of molecular graphs.
///
/// Implements the FxHash algorithm: each word is folded in by a rotate, an XOR,
/// and a multiply by [`MULTIPLIER`]. The state carries no random seed, so equal
/// keys hash equally across runs.
///
/// Use through [`FxHashSet`] and [`FxHashMap`] rather than directly.
#[derive(Default)]
pub struct FxHasher {
    hash: u64,
}

impl FxHasher {
    /// Folds one word into the running hash.
    #[inline]
    fn absorb(&mut self, word: u64) {
        self.hash = (self.hash.rotate_left(5) ^ word).wrapping_mul(MULTIPLIER);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.absorb(byte as u64);
        }
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.absorb(i as u64);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.absorb(i as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

/// A [`HashSet`] keyed by [`FxHasher`].
pub type FxHashSet<T> = HashSet<T, BuildHasherDefault<FxHasher>>;

/// A [`HashMap`] keyed by [`FxHasher`].
pub type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::hash::{Hash, Hasher};

    fn hash(value: impl Hash) -> u64 {
        let mut hasher = FxHasher::default();
        value.hash(&mut hasher);
        hasher.finish()
    }

    #[test]
    fn default_set_is_empty() {
        let set: FxHashSet<u32> = FxHashSet::default();
        assert!(set.is_empty());
    }

    #[test]
    fn set_reports_an_inserted_key_as_present() {
        let mut set: FxHashSet<u32> = FxHashSet::default();
        set.insert(7);
        assert!(set.contains(&7));
    }

    #[test]
    fn map_returns_the_value_for_an_inserted_key() {
        let mut map: FxHashMap<u32, &str> = FxHashMap::default();
        map.insert(1, "a");
        assert_eq!(map.get(&1), Some(&"a"));
    }

    #[test]
    fn set_reports_an_absent_key_as_missing() {
        let mut set: FxHashSet<u32> = FxHashSet::default();
        set.insert(7);
        assert!(!set.contains(&8));
    }

    #[test]
    fn map_returns_none_for_an_absent_key() {
        let mut map: FxHashMap<u32, &str> = FxHashMap::default();
        map.insert(1, "a");
        assert_eq!(map.get(&2), None);
    }

    #[test]
    fn distinct_values_hash_differently() {
        assert_ne!(hash(5u32), hash(6u32));
    }

    #[test]
    fn distinct_byte_strings_hash_differently() {
        assert_ne!(hash("ab"), hash("cd"));
    }

    #[test]
    fn reordering_a_compound_value_changes_the_hash() {
        assert_ne!(hash((1u32, 2u32)), hash((2u32, 1u32)));
    }

    #[test]
    fn set_holds_every_distinct_key() {
        let mut set: FxHashSet<u32> = FxHashSet::default();
        for key in [1, 2, 3, 100, 5000] {
            set.insert(key);
        }
        assert_eq!(set.len(), 5);
        assert!([1, 2, 3, 100, 5000].iter().all(|key| set.contains(key)));
    }

    #[test]
    fn hashing_is_deterministic() {
        assert_eq!(hash(0x1234_5678_u32), hash(0x1234_5678_u32));
    }
}

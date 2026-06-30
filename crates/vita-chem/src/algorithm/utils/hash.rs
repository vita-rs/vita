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

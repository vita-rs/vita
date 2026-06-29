/// A map backed by a sorted array of key–value pairs.
///
/// Stores entries sorted by key in one contiguous `Vec` and answers lookups by
/// binary search. For the small key sets of molecular graphs this matches a
/// hash map's lookup speed while staying allocation-light, fully deterministic,
/// and free of any hashing.
///
/// Keys are expected to be unique; built from duplicate keys, [`get`](Self::get)
/// returns the value of an unspecified one.
///
/// Obtain via [`from_pairs`](Self::from_pairs).
pub struct SortedMap<K, V> {
    entries: Vec<(K, V)>,
}

impl<K: Ord, V> SortedMap<K, V> {
    /// Builds a map from key–value pairs, sorting them by key.
    pub fn from_pairs(pairs: impl IntoIterator<Item = (K, V)>) -> Self {
        let mut entries: Vec<(K, V)> = pairs.into_iter().collect();
        entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        SortedMap { entries }
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if the map has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Returns a reference to the value for `key`.
    ///
    /// Returns `None` if `key` is absent.
    pub fn get(&self, key: &K) -> Option<&V> {
        self.entries
            .binary_search_by(|(k, _)| k.cmp(key))
            .ok()
            .map(|i| &self.entries[i].1)
    }

    /// Returns `true` if `key` is present.
    pub fn contains_key(&self, key: &K) -> bool {
        self.entries.binary_search_by(|(k, _)| k.cmp(key)).is_ok()
    }

    /// Iterates the entries as `(key, value)` references, ordered by key.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> + '_ {
        self.entries.iter().map(|(k, v)| (k, v))
    }
}

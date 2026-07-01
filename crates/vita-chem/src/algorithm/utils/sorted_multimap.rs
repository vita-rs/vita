/// A multimap backed by sorted arrays, associating each key with several values.
///
/// Concatenates the values of every key into one contiguous `Vec`, addressed by
/// a sorted array of the distinct keys, and answers lookups by binary search.
/// Keys and the values under each key are both kept in ascending order, so the
/// contents depend only on the set of pairs supplied, not on their order. For
/// the small key sets of molecular graphs this matches a hash map's lookup speed
/// while staying allocation-light, fully deterministic, and free of any hashing.
///
/// Obtain via [`from_pairs`](Self::from_pairs).
pub struct SortedMultimap<K, V> {
    keys: Vec<K>,
    offsets: Vec<usize>,
    values: Vec<V>,
}

impl<K: Ord, V: Ord> SortedMultimap<K, V> {
    /// Builds a multimap from key–value pairs, grouping each key's values in
    /// ascending order.
    pub fn from_pairs(pairs: impl IntoIterator<Item = (K, V)>) -> Self {
        let mut pairs: Vec<(K, V)> = pairs.into_iter().collect();
        pairs.sort_unstable();

        let mut keys = Vec::new();
        let mut offsets = Vec::new();
        let mut values = Vec::with_capacity(pairs.len());
        for (key, value) in pairs {
            if keys.last() != Some(&key) {
                offsets.push(values.len());
                keys.push(key);
            }
            values.push(value);
        }
        offsets.push(values.len());

        SortedMultimap {
            keys,
            offsets,
            values,
        }
    }

    /// Number of distinct keys.
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Returns `true` if the multimap has no keys.
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Returns the values associated with `key`, in ascending order.
    ///
    /// Returns an empty slice if `key` is absent.
    pub fn get(&self, key: &K) -> &[V] {
        match self.keys.binary_search_by(|k| k.cmp(key)) {
            Ok(i) => &self.values[self.offsets[i]..self.offsets[i + 1]],
            Err(_) => &[],
        }
    }

    /// Returns `true` if `key` is present.
    pub fn contains_key(&self, key: &K) -> bool {
        self.keys.binary_search_by(|k| k.cmp(key)).is_ok()
    }

    /// Iterates each key with its values as `(key, values)`, ordered by key.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &[V])> + '_ {
        (0..self.keys.len()).map(move |i| {
            (
                &self.keys[i],
                &self.values[self.offsets[i]..self.offsets[i + 1]],
            )
        })
    }
}

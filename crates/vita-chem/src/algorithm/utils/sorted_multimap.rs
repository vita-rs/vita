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

#[cfg(test)]
mod tests {
    use super::*;

    fn empty() -> SortedMultimap<i32, char> {
        SortedMultimap::from_pairs(Vec::new())
    }

    fn single() -> SortedMultimap<i32, char> {
        SortedMultimap::from_pairs([(5, 'e')])
    }

    fn grouped() -> SortedMultimap<i32, char> {
        SortedMultimap::from_pairs([(3, 'c'), (1, 'd'), (2, 'b'), (1, 'a')])
    }

    fn gapped() -> SortedMultimap<i32, char> {
        SortedMultimap::from_pairs([(1, 'a'), (4, 'e')])
    }

    #[test]
    fn empty_multimap_has_no_keys() {
        assert_eq!(empty().len(), 0);
    }

    #[test]
    fn empty_multimap_is_empty() {
        assert!(empty().is_empty());
    }

    #[test]
    fn empty_multimap_get_yields_no_values() {
        assert!(empty().get(&0).is_empty());
    }

    #[test]
    fn single_key_multimap_has_one_key() {
        assert_eq!(single().len(), 1);
    }

    #[test]
    fn single_key_multimap_is_not_empty() {
        assert!(!single().is_empty());
    }

    #[test]
    fn get_returns_all_of_a_keys_values_in_ascending_order() {
        assert_eq!(grouped().get(&1), &['a', 'd']);
    }

    #[test]
    fn contains_key_is_true_for_a_present_key() {
        assert!(grouped().contains_key(&2));
    }

    #[test]
    fn get_yields_no_values_for_a_key_in_a_gap() {
        assert!(gapped().get(&2).is_empty());
    }

    #[test]
    fn contains_key_is_false_for_an_absent_key() {
        assert!(!gapped().contains_key(&2));
    }

    #[test]
    fn get_finds_the_minimum_key() {
        assert_eq!(gapped().get(&1), &['a']);
    }

    #[test]
    fn get_finds_the_maximum_key() {
        assert_eq!(gapped().get(&4), &['e']);
    }

    #[test]
    fn get_yields_no_values_below_the_minimum_key() {
        assert!(gapped().get(&0).is_empty());
    }

    #[test]
    fn get_yields_no_values_above_the_maximum_key() {
        assert!(gapped().get(&5).is_empty());
    }

    #[test]
    fn iter_yields_each_key_with_its_values_in_key_order() {
        let mm = grouped();
        let entries: Vec<(&i32, &[char])> = mm.iter().collect();
        assert_eq!(
            entries,
            vec![(&1, &['a', 'd'][..]), (&2, &['b'][..]), (&3, &['c'][..])],
        );
    }

    #[test]
    fn len_counts_distinct_keys_not_values() {
        assert_eq!(grouped().len(), 3);
    }

    #[test]
    fn lookup_is_independent_of_input_order() {
        let a = SortedMultimap::from_pairs([(3, 'c'), (1, 'd'), (2, 'b'), (1, 'a')]);
        let b = SortedMultimap::from_pairs([(1, 'a'), (2, 'b'), (1, 'd'), (3, 'c')]);
        assert_eq!(a.iter().collect::<Vec<_>>(), b.iter().collect::<Vec<_>>());
        for key in [1, 2, 3] {
            assert_eq!(a.get(&key), b.get(&key), "values for {key} differ");
        }
    }
}

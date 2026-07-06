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
#[derive(Clone, Debug, PartialEq, Eq)]
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
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if the map has no entries.
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn contains_key(&self, key: &K) -> bool {
        self.entries.binary_search_by(|(k, _)| k.cmp(key)).is_ok()
    }

    /// Iterates the entries as `(key, value)` references, ordered by key.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> + '_ {
        self.entries.iter().map(|(k, v)| (k, v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty() -> SortedMap<i32, char> {
        SortedMap::from_pairs(Vec::new())
    }

    fn single() -> SortedMap<i32, char> {
        SortedMap::from_pairs([(5, 'e')])
    }

    fn unsorted() -> SortedMap<i32, char> {
        SortedMap::from_pairs([(3, 'c'), (1, 'a'), (2, 'b')])
    }

    fn gapped() -> SortedMap<i32, char> {
        SortedMap::from_pairs([(1, 'a'), (5, 'e'), (3, 'c')])
    }

    #[test]
    fn empty_map_has_no_entries() {
        assert_eq!(empty().len(), 0);
    }

    #[test]
    fn empty_map_is_empty() {
        assert!(empty().is_empty());
    }

    #[test]
    fn empty_map_get_is_none() {
        assert_eq!(empty().get(&0), None);
    }

    #[test]
    fn single_entry_map_has_length_one() {
        assert_eq!(single().len(), 1);
    }

    #[test]
    fn single_entry_map_is_not_empty() {
        assert!(!single().is_empty());
    }

    #[test]
    fn get_returns_value_for_present_key() {
        assert_eq!(single().get(&5), Some(&'e'));
    }

    #[test]
    fn contains_key_is_true_for_present_key() {
        assert!(single().contains_key(&5));
    }

    #[test]
    fn get_returns_none_for_absent_key() {
        assert_eq!(single().get(&9), None);
    }

    #[test]
    fn contains_key_is_false_for_absent_key() {
        assert!(!single().contains_key(&9));
    }

    #[test]
    fn get_is_none_for_key_in_a_gap() {
        assert_eq!(gapped().get(&2), None);
    }

    #[test]
    fn get_finds_minimum_key() {
        assert_eq!(gapped().get(&1), Some(&'a'));
    }

    #[test]
    fn get_finds_maximum_key() {
        assert_eq!(gapped().get(&5), Some(&'e'));
    }

    #[test]
    fn get_is_none_below_minimum_key() {
        assert_eq!(gapped().get(&0), None);
    }

    #[test]
    fn get_is_none_above_maximum_key() {
        assert_eq!(gapped().get(&6), None);
    }

    #[test]
    fn iter_yields_entries_in_ascending_key_order() {
        let map = unsorted();
        let entries: Vec<(&i32, &char)> = map.iter().collect();
        assert_eq!(entries, vec![(&1, &'a'), (&2, &'b'), (&3, &'c')]);
    }

    #[test]
    fn get_resolves_each_key_to_its_value() {
        let map = unsorted();
        assert_eq!(map.get(&1), Some(&'a'));
        assert_eq!(map.get(&2), Some(&'b'));
        assert_eq!(map.get(&3), Some(&'c'));
    }

    #[test]
    fn len_counts_all_entries() {
        assert_eq!(unsorted().len(), 3);
        assert_eq!(gapped().len(), 3);
    }

    #[test]
    fn lookup_is_independent_of_insertion_order() {
        let a = SortedMap::from_pairs([(3, 'c'), (1, 'a'), (2, 'b')]);
        let b = SortedMap::from_pairs([(2, 'b'), (3, 'c'), (1, 'a')]);
        assert_eq!(
            a.iter().collect::<Vec<_>>(),
            b.iter().collect::<Vec<_>>(),
            "iteration order differs",
        );
        for k in [1, 2, 3] {
            assert_eq!(a.get(&k), b.get(&k), "lookup of {k} differs");
        }
    }
}

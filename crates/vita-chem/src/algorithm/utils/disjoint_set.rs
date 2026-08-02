/// A partition of the integers `0..n` into disjoint sets.
///
/// Starts with `n` singletons and merges them with [`union`](Self::union);
/// [`find`](Self::find) reports which set an element belongs to by returning a
/// representative shared by exactly that set's elements. Path halving and union
/// by size keep every operation within O(α(n)) amortised time, where α is the
/// inverse Ackermann function — effectively constant.
///
/// Elements are the integers `0..n`; every method taking an element requires it
/// to be less than [`len`](Self::len).
///
/// Obtain via [`new`](Self::new).
#[derive(Debug)]
pub struct DisjointSet {
    parent: Vec<usize>,
    size: Vec<usize>,
    sets: usize,
}

impl DisjointSet {
    /// Creates a partition of `0..n` into `n` singleton sets.
    pub fn new(n: usize) -> Self {
        DisjointSet {
            parent: (0..n).collect(),
            size: vec![1; n],
            sets: n,
        }
    }

    /// The number of elements.
    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    /// Returns `true` if there are no elements.
    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    /// The number of disjoint sets.
    #[cfg(test)]
    pub fn set_count(&self) -> usize {
        self.sets
    }

    /// The representative of the set containing `x`, shared by exactly that
    /// set's elements.
    pub fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    /// Merges the sets containing `a` and `b`, returning `true` if they were
    /// previously separate.
    pub fn union(&mut self, a: usize, b: usize) -> bool {
        let a = self.find(a);
        let b = self.find(b);
        if a == b {
            return false;
        }
        let (small, large) = if self.size[a] < self.size[b] {
            (a, b)
        } else {
            (b, a)
        };
        self.parent[small] = large;
        self.size[large] += self.size[small];
        self.sets -= 1;
        true
    }

    /// Returns `true` if `a` and `b` are in the same set.
    #[cfg(test)]
    pub fn connected(&mut self, a: usize, b: usize) -> bool {
        self.find(a) == self.find(b)
    }

    /// The sets, each as its elements in ascending order, ordered by least
    /// element.
    pub fn groups(&mut self) -> Vec<Vec<usize>> {
        let n = self.parent.len();
        let mut group_of = vec![usize::MAX; n];
        let mut groups: Vec<Vec<usize>> = Vec::with_capacity(self.sets);
        for i in 0..n {
            let root = self.find(i);
            let g = if group_of[root] == usize::MAX {
                group_of[root] = groups.len();
                groups.push(Vec::new());
                group_of[root]
            } else {
                group_of[root]
            };
            groups[g].push(i);
        }
        groups
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one_pair() -> DisjointSet {
        let mut ds = DisjointSet::new(2);
        ds.union(0, 1);
        ds
    }

    fn two_pairs() -> DisjointSet {
        let mut ds = DisjointSet::new(4);
        ds.union(0, 1);
        ds.union(2, 3);
        ds
    }

    #[test]
    fn empty_universe_has_length_zero() {
        assert_eq!(DisjointSet::new(0).len(), 0);
    }

    #[test]
    fn empty_universe_is_empty() {
        assert!(DisjointSet::new(0).is_empty());
    }

    #[test]
    fn empty_universe_has_no_sets() {
        assert_eq!(DisjointSet::new(0).set_count(), 0);
    }

    #[test]
    fn empty_universe_has_no_groups() {
        assert!(DisjointSet::new(0).groups().is_empty());
    }

    #[test]
    fn lone_element_is_its_own_representative() {
        assert_eq!(DisjointSet::new(1).find(0), 0);
    }

    #[test]
    fn universe_length_is_the_element_count() {
        assert_eq!(DisjointSet::new(3).len(), 3);
    }

    #[test]
    fn nonempty_universe_is_not_empty() {
        assert!(!DisjointSet::new(3).is_empty());
    }

    #[test]
    fn fresh_universe_has_one_set_per_element() {
        assert_eq!(DisjointSet::new(3).set_count(), 3);
    }

    #[test]
    fn union_of_separate_sets_returns_true() {
        let mut ds = DisjointSet::new(2);
        assert!(ds.union(0, 1));
    }

    #[test]
    fn union_connects_its_two_elements() {
        let mut ds = one_pair();
        assert!(ds.connected(0, 1));
    }

    #[test]
    fn joined_elements_share_a_representative() {
        let mut ds = one_pair();
        assert_eq!(ds.find(0), ds.find(1));
    }

    #[test]
    fn union_reduces_the_set_count_by_one() {
        let mut ds = DisjointSet::new(3);
        ds.union(0, 1);
        assert_eq!(ds.set_count(), 2);
    }

    #[test]
    fn distinct_elements_are_not_connected() {
        let mut ds = DisjointSet::new(2);
        assert!(!ds.connected(0, 1));
    }

    #[test]
    fn union_of_already_joined_sets_returns_false() {
        let mut ds = one_pair();
        assert!(!ds.union(0, 1));
    }

    #[test]
    fn redundant_union_leaves_the_set_count_unchanged() {
        let mut ds = one_pair();
        ds.union(0, 1);
        assert_eq!(ds.set_count(), 1);
    }

    #[test]
    fn an_element_is_connected_to_itself() {
        let mut ds = DisjointSet::new(1);
        assert!(ds.connected(0, 0));
    }

    #[test]
    fn self_union_returns_false() {
        let mut ds = DisjointSet::new(2);
        assert!(!ds.union(0, 0));
    }

    #[test]
    fn merging_a_chain_yields_a_single_set() {
        let mut ds = DisjointSet::new(4);
        ds.union(0, 1);
        ds.union(1, 2);
        ds.union(2, 3);
        assert_eq!(ds.set_count(), 1);
    }

    #[test]
    fn connectivity_is_transitive() {
        let mut ds = DisjointSet::new(3);
        ds.union(0, 1);
        ds.union(1, 2);
        assert!(ds.connected(0, 2));
    }

    #[test]
    fn elements_of_different_sets_are_not_connected() {
        let mut ds = two_pairs();
        assert!(!ds.connected(0, 2));
    }

    #[test]
    fn disjoint_unions_yield_separate_sets() {
        assert_eq!(two_pairs().set_count(), 2);
    }

    #[test]
    fn groups_are_the_connected_sets() {
        let mut ds = DisjointSet::new(6);
        ds.union(0, 2);
        ds.union(2, 4);
        ds.union(1, 3);
        assert_eq!(ds.groups(), vec![vec![0, 2, 4], vec![1, 3], vec![5]]);
    }

    #[test]
    fn groups_are_independent_of_union_order() {
        let mut forward = DisjointSet::new(5);
        forward.union(0, 1);
        forward.union(1, 2);
        forward.union(3, 4);

        let mut reversed = DisjointSet::new(5);
        reversed.union(4, 3);
        reversed.union(2, 1);
        reversed.union(1, 0);

        assert_eq!(forward.groups(), reversed.groups());
    }
}

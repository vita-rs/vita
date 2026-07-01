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
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    /// Returns `true` if there are no elements.
    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    /// The number of disjoint sets.
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

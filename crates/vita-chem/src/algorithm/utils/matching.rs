use std::collections::VecDeque;

/// Absent vertex: an unmatched mate or a missing tree parent.
const NIL: usize = usize::MAX;

/// A maximum-cardinality matching of an undirected graph.
///
/// Pairs vertices across edges so that no vertex lies in two pairs and no other
/// matching pairs up more. A perfect matching — one leaving no vertex unpaired —
/// is therefore found whenever the graph admits one.
///
/// Obtain via [`maximum_matching`].
#[derive(Debug)]
pub struct MaximumMatching {
    mate: Vec<Option<usize>>,
}

impl MaximumMatching {
    /// The vertex matched with `vertex`, or `None` if it is left unmatched.
    pub fn mate(&self, vertex: usize) -> Option<usize> {
        self.mate[vertex]
    }

    /// Number of matched pairs.
    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.mate.iter().flatten().count() / 2
    }

    /// Returns `true` if no vertex is matched.
    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.mate.iter().all(Option::is_none)
    }

    /// Returns `true` if every vertex is matched.
    pub fn is_perfect(&self) -> bool {
        self.mate.iter().all(Option::is_some)
    }
}

/// The maximum-cardinality matching of a graph given as an adjacency list over
/// its vertices `0..adjacency.len()`.
///
/// Runs Edmonds' blossom algorithm: it grows an alternating tree from each
/// unmatched vertex and, on meeting another, flips the path between them to add
/// a pair. The odd cycles ("blossoms") that foil a naive search are contracted
/// to a single vertex and expanded once matched, so the matching is maximum even
/// on the non-bipartite graphs — odd rings and fused systems — that chemistry
/// poses.
///
/// # Complexity
///
/// O(V³) time and O(V) space, over the graph's `V` vertices.
pub fn maximum_matching(adjacency: &[Vec<usize>]) -> MaximumMatching {
    let n = adjacency.len();
    let mut blossom = Blossom {
        adjacency,
        mate: vec![NIL; n],
        parent: vec![NIL; n],
        base: vec![0; n],
        even: vec![false; n],
        in_blossom: vec![false; n],
        queue: VecDeque::new(),
    };
    for root in 0..n {
        if blossom.mate[root] == NIL {
            blossom.grow(root);
        }
    }
    MaximumMatching {
        mate: blossom
            .mate
            .into_iter()
            .map(|v| (v != NIL).then_some(v))
            .collect(),
    }
}

/// The working state of one blossom run over a borrowed graph.
struct Blossom<'a> {
    adjacency: &'a [Vec<usize>],
    /// The vertex each is matched with, or [`NIL`].
    mate: Vec<usize>,
    /// The alternating-tree parent of each vertex, or [`NIL`].
    parent: Vec<usize>,
    /// The base each vertex's contracted blossom hangs from.
    base: Vec<usize>,
    /// Whether each vertex sits at even depth in the tree (an outer vertex).
    even: Vec<bool>,
    /// Scratch flags marking the vertices of the blossom being contracted.
    in_blossom: Vec<bool>,
    /// Outer vertices still to expand.
    queue: VecDeque<usize>,
}

impl Blossom<'_> {
    /// Grows an alternating tree from the unmatched `root`. On reaching another
    /// unmatched vertex it flips the augmenting path between them, adding one
    /// pair to the matching; otherwise the tree exhausts and `root` stays free.
    fn grow(&mut self, root: usize) {
        self.even.iter_mut().for_each(|e| *e = false);
        self.parent.iter_mut().for_each(|p| *p = NIL);
        self.base.iter_mut().enumerate().for_each(|(i, b)| *b = i);
        self.even[root] = true;
        self.queue.clear();
        self.queue.push_back(root);

        while let Some(v) = self.queue.pop_front() {
            for i in 0..self.adjacency[v].len() {
                let to = self.adjacency[v][i];
                if self.base[v] == self.base[to] || self.mate[v] == to {
                    continue;
                }
                if to == root || (self.mate[to] != NIL && self.parent[self.mate[to]] != NIL) {
                    let base = self.lca(v, to);
                    self.in_blossom.iter_mut().for_each(|f| *f = false);
                    self.mark_path(v, base, to);
                    self.mark_path(to, base, v);
                    for u in 0..self.adjacency.len() {
                        if self.in_blossom[self.base[u]] {
                            self.base[u] = base;
                            if !self.even[u] {
                                self.even[u] = true;
                                self.queue.push_back(u);
                            }
                        }
                    }
                } else if self.parent[to] == NIL {
                    self.parent[to] = v;
                    if self.mate[to] == NIL {
                        self.augment(to);
                        return;
                    }
                    self.even[self.mate[to]] = true;
                    self.queue.push_back(self.mate[to]);
                }
            }
        }
    }

    /// Returns the base of the lowest common ancestor of `a` and `b` in the
    /// alternating tree — the base the contracted blossom hangs from.
    fn lca(&self, mut a: usize, mut b: usize) -> usize {
        let mut seen = vec![false; self.adjacency.len()];
        loop {
            a = self.base[a];
            seen[a] = true;
            if self.mate[a] == NIL {
                break;
            }
            a = self.parent[self.mate[a]];
        }
        loop {
            b = self.base[b];
            if seen[b] {
                return b;
            }
            b = self.parent[self.mate[b]];
        }
    }

    /// Walks the alternating path from `v` up to the blossom `base`, flagging the
    /// base of every vertex on it and threading parent pointers back toward
    /// `child` so the contracted cycle stays traversable.
    fn mark_path(&mut self, mut v: usize, base: usize, mut child: usize) {
        while self.base[v] != base {
            self.in_blossom[self.base[v]] = true;
            self.in_blossom[self.base[self.mate[v]]] = true;
            self.parent[v] = child;
            child = self.mate[v];
            v = self.parent[self.mate[v]];
        }
    }

    /// Flips matched and unmatched edges along the augmenting path from the newly
    /// reached unmatched vertex `v` back to the root.
    fn augment(&mut self, mut v: usize) {
        while v != NIL {
            let parent = self.parent[v];
            let next = self.mate[parent];
            self.mate[v] = parent;
            self.mate[parent] = v;
            v = next;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty() -> Vec<Vec<usize>> {
        vec![]
    }

    fn isolated() -> Vec<Vec<usize>> {
        vec![vec![]]
    }

    fn edge() -> Vec<Vec<usize>> {
        vec![vec![1], vec![0]]
    }

    fn even_cycle() -> Vec<Vec<usize>> {
        vec![vec![1, 3], vec![0, 2], vec![1, 3], vec![0, 2]]
    }

    fn odd_cycle() -> Vec<Vec<usize>> {
        vec![vec![1, 2], vec![0, 2], vec![0, 1]]
    }

    fn blossom() -> Vec<Vec<usize>> {
        vec![vec![1, 2, 3], vec![0, 2], vec![0, 1], vec![0]]
    }

    fn disconnected_edges() -> Vec<Vec<usize>> {
        vec![vec![1], vec![0], vec![3], vec![2]]
    }

    fn joined_triangles() -> Vec<Vec<usize>> {
        vec![
            vec![1, 2],
            vec![0, 2],
            vec![0, 1, 3],
            vec![2, 4, 5],
            vec![3, 5],
            vec![3, 4],
        ]
    }

    #[test]
    fn empty_graph_matches_nothing() {
        let m = maximum_matching(&empty());
        assert!(m.is_empty());
        assert_eq!(m.len(), 0);
    }

    #[test]
    fn isolated_vertex_is_unmatched() {
        assert_eq!(maximum_matching(&isolated()).mate(0), None);
    }

    #[test]
    fn single_edge_matches_its_two_vertices() {
        let m = maximum_matching(&edge());
        assert_eq!(m.mate(0), Some(1));
        assert_eq!(m.mate(1), Some(0));
        assert_eq!(m.len(), 1);
    }

    #[test]
    fn even_cycle_is_perfectly_matched() {
        let m = maximum_matching(&even_cycle());
        assert!(m.is_perfect());
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn odd_cycle_leaves_one_vertex_unmatched() {
        let m = maximum_matching(&odd_cycle());
        assert!(!m.is_perfect());
        assert_eq!(m.len(), 1);
    }

    #[test]
    fn blossom_is_contracted_to_reach_a_perfect_matching() {
        let m = maximum_matching(&blossom());
        assert!(m.is_perfect());
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn disconnected_components_are_each_matched() {
        let m = maximum_matching(&disconnected_edges());
        assert!(m.is_perfect());
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn joined_triangles_are_perfectly_matched() {
        let m = maximum_matching(&joined_triangles());
        assert!(m.is_perfect());
        assert_eq!(m.len(), 3);
    }

    #[test]
    fn every_matched_pair_is_a_symmetric_edge() {
        let graph = joined_triangles();
        let m = maximum_matching(&graph);
        for (v, _) in graph.iter().enumerate() {
            if let Some(u) = m.mate(v) {
                assert_eq!(m.mate(u), Some(v), "pairing of {v} and {u} is asymmetric");
                assert!(
                    graph[v].contains(&u),
                    "{v} and {u} are matched but not adjacent"
                );
            }
        }
    }
}

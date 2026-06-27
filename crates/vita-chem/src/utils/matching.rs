use std::collections::VecDeque;

/// No vertex: an unmatched mate, or a missing tree parent.
const NIL: usize = usize::MAX;

/// A maximum-cardinality matching of an undirected graph.
///
/// Takes the graph as an adjacency list over the vertices `0..adjacency.len()`
/// and returns, for each vertex, the vertex it is matched with, or `None` where
/// it is left unmatched. The result is symmetric: vertex `a` is matched with `b`
/// exactly when `b` is matched with `a`. The matching is maximum — no other
/// matching pairs up more vertices — so a perfect one is found whenever the
/// graph admits it.
///
/// # Complexity
///
/// O(V³) time, O(V) space.
pub(crate) fn maximum_matching(adjacency: &[Vec<usize>]) -> Vec<Option<usize>> {
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
    blossom
        .mate
        .into_iter()
        .map(|v| (v != NIL).then_some(v))
        .collect()
}

/// The working state of one run of the blossom algorithm.
struct Blossom<'a> {
    adjacency: &'a [Vec<usize>],
    mate: Vec<usize>,
    parent: Vec<usize>,
    base: Vec<usize>,
    even: Vec<bool>,
    in_blossom: Vec<bool>,
    queue: VecDeque<usize>,
}

impl Blossom<'_> {
    /// Grows an alternating tree from the unmatched `root`. On reaching another
    /// unmatched vertex it flips the augmenting path between them, adding one
    /// edge to the matching; otherwise the tree exhausts and `root` stays free.
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

    fn matched(matching: &[Option<usize>]) -> usize {
        matching.iter().filter(|m| m.is_some()).count()
    }

    fn is_symmetric(matching: &[Option<usize>]) -> bool {
        matching.iter().enumerate().all(|(v, &m)| match m {
            Some(u) => matching[u] == Some(v),
            None => true,
        })
    }

    #[test]
    fn empty_graph_has_no_matching() {
        assert_eq!(maximum_matching(&[]), Vec::<Option<usize>>::new());
    }

    #[test]
    fn isolated_vertex_is_unmatched() {
        assert_eq!(maximum_matching(&[vec![]]), vec![None]);
    }

    #[test]
    fn single_edge_is_matched() {
        assert_eq!(
            maximum_matching(&[vec![1], vec![0]]),
            vec![Some(1), Some(0)]
        );
    }

    #[test]
    fn path_matches_all_but_one() {
        let m = maximum_matching(&[vec![1], vec![0, 2], vec![1]]);
        assert_eq!(matched(&m), 2);
        assert!(is_symmetric(&m));
    }

    #[test]
    fn even_cycle_is_perfectly_matched() {
        let m = maximum_matching(&[vec![1, 3], vec![0, 2], vec![1, 3], vec![0, 2]]);
        assert_eq!(matched(&m), 4);
        assert!(is_symmetric(&m));
    }

    #[test]
    fn odd_cycle_leaves_one_unmatched() {
        let m = maximum_matching(&[vec![1, 2], vec![0, 2], vec![0, 1]]);
        assert_eq!(matched(&m), 2);
        assert!(is_symmetric(&m));
    }

    #[test]
    fn blossom_is_contracted_to_match_a_pendant() {
        let m = maximum_matching(&[vec![1, 2, 3], vec![0, 2], vec![0, 1], vec![0]]);
        assert_eq!(matched(&m), 4);
        assert!(is_symmetric(&m));
    }

    #[test]
    fn joined_triangles_are_perfectly_matched() {
        let m = maximum_matching(&[
            vec![1, 2],
            vec![0, 2],
            vec![0, 1, 3],
            vec![2, 4, 5],
            vec![3, 5],
            vec![3, 4],
        ]);
        assert_eq!(matched(&m), 6);
        assert!(is_symmetric(&m));
    }
}

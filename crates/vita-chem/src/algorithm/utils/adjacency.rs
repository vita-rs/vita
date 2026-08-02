/// An indexed undirected adjacency list.
///
/// Represents an undirected graph over integer nodes `0..n`. Each undirected
/// edge supplied to [`build`](Self::build) is stored in both directions:
/// node `a` records `(edge, b)` and node `b` records `(edge, a)`.
///
/// Obtain via [`build`](Self::build).
#[derive(Clone, Debug)]
pub struct AdjacencyList {
    adj: Vec<Vec<(usize, usize)>>,
}

impl AdjacencyList {
    /// Build an adjacency list over `n` nodes from an edge iterator.
    ///
    /// Each yielded triple `(edge, a, b)` inserts the undirected edge between
    /// nodes `a` and `b` with identifier `edge`. Both `a` and `b` must be
    /// less than `n`.
    pub fn build(n: usize, edges: impl IntoIterator<Item = (usize, usize, usize)>) -> Self {
        let mut adj = vec![vec![]; n];
        for (e, a, b) in edges {
            adj[a].push((e, b));
            adj[b].push((e, a));
        }
        AdjacencyList { adj }
    }

    /// Number of nodes.
    pub fn len(&self) -> usize {
        self.adj.len()
    }

    /// Returns `true` if the graph has no nodes.
    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.adj.is_empty()
    }

    /// Edges incident to node `i` as `(edge, neighbor)` pairs.
    ///
    /// The order reflects the order in which edges were supplied to
    /// [`build`](Self::build).
    pub fn neighbors(&self, i: usize) -> &[(usize, usize)] {
        &self.adj[i]
    }
}

impl std::ops::Index<usize> for AdjacencyList {
    type Output = [(usize, usize)];

    fn index(&self, i: usize) -> &Self::Output {
        &self.adj[i]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty() -> AdjacencyList {
        AdjacencyList::build(0, [])
    }

    fn isolated(n: usize) -> AdjacencyList {
        AdjacencyList::build(n, [])
    }

    fn single_edge() -> AdjacencyList {
        AdjacencyList::build(2, [(0, 0, 1)])
    }

    fn chain() -> AdjacencyList {
        AdjacencyList::build(3, [(0, 0, 1), (1, 1, 2)])
    }

    fn triangle() -> AdjacencyList {
        AdjacencyList::build(3, [(0, 0, 1), (1, 1, 2), (2, 0, 2)])
    }

    fn star() -> AdjacencyList {
        AdjacencyList::build(4, [(0, 0, 1), (1, 0, 2), (2, 0, 3)])
    }

    #[test]
    fn empty_graph_has_no_nodes() {
        assert_eq!(empty().len(), 0);
    }

    #[test]
    fn empty_graph_is_empty() {
        assert!(empty().is_empty());
    }

    #[test]
    fn single_node_is_not_empty() {
        assert!(!isolated(1).is_empty());
    }

    #[test]
    fn single_node_has_no_neighbors() {
        assert!(isolated(1).neighbors(0).is_empty());
    }

    #[test]
    fn edgeless_graph_has_no_neighbors_at_any_node() {
        let g = isolated(3);
        assert!(g.neighbors(0).is_empty());
        assert!(g.neighbors(1).is_empty());
        assert!(g.neighbors(2).is_empty());
    }

    #[test]
    fn node_count_matches_n() {
        assert_eq!(isolated(0).len(), 0);
        assert_eq!(isolated(1).len(), 1);
        assert_eq!(isolated(5).len(), 5);
        assert_eq!(triangle().len(), 3);
    }

    #[test]
    fn edge_appears_at_first_endpoint() {
        assert!(single_edge().neighbors(0).contains(&(0, 1)));
    }

    #[test]
    fn edge_appears_at_second_endpoint() {
        assert!(single_edge().neighbors(1).contains(&(0, 0)));
    }

    #[test]
    fn edge_id_is_preserved_at_first_endpoint() {
        let g = AdjacencyList::build(2, [(42, 0, 1)]);
        assert!(g.neighbors(0).contains(&(42, 1)));
    }

    #[test]
    fn edge_id_is_preserved_at_second_endpoint() {
        let g = AdjacencyList::build(2, [(42, 0, 1)]);
        assert!(g.neighbors(1).contains(&(42, 0)));
    }

    #[test]
    fn disconnected_node_has_no_neighbors() {
        let g = AdjacencyList::build(3, [(0, 0, 1)]);
        assert!(g.neighbors(2).is_empty());
    }

    #[test]
    fn len_counts_nodes_not_edges() {
        let g = AdjacencyList::build(5, [(0, 0, 1), (1, 2, 3)]);
        assert_eq!(g.len(), 5);
    }

    #[test]
    fn single_edge_first_endpoint_has_one_neighbor() {
        assert_eq!(single_edge().neighbors(0).len(), 1);
    }

    #[test]
    fn single_edge_second_endpoint_has_one_neighbor() {
        assert_eq!(single_edge().neighbors(1).len(), 1);
    }

    #[test]
    fn self_loop_appears_twice_at_same_node() {
        let g = AdjacencyList::build(1, [(0, 0, 0)]);
        assert_eq!(g.neighbors(0).len(), 2);
        assert_eq!(g.neighbors(0)[0], (0, 0));
        assert_eq!(g.neighbors(0)[1], (0, 0));
    }

    #[test]
    fn triangle_each_node_has_two_neighbors() {
        let g = triangle();
        assert_eq!(g.neighbors(0).len(), 2);
        assert_eq!(g.neighbors(1).len(), 2);
        assert_eq!(g.neighbors(2).len(), 2);
    }

    #[test]
    fn triangle_node_zero_neighbor_set_is_correct() {
        let mut n: Vec<_> = triangle().neighbors(0).to_vec();
        n.sort();
        assert_eq!(n, [(0, 1), (2, 2)]);
    }

    #[test]
    fn triangle_node_one_neighbor_set_is_correct() {
        let mut n: Vec<_> = triangle().neighbors(1).to_vec();
        n.sort();
        assert_eq!(n, [(0, 0), (1, 2)]);
    }

    #[test]
    fn chain_endpoints_have_one_neighbor() {
        assert_eq!(chain().neighbors(0).len(), 1);
        assert_eq!(chain().neighbors(2).len(), 1);
    }

    #[test]
    fn chain_internal_node_has_two_neighbors() {
        assert_eq!(chain().neighbors(1).len(), 2);
    }

    #[test]
    fn star_center_degree_equals_leaf_count() {
        assert_eq!(star().neighbors(0).len(), 3);
    }

    #[test]
    fn star_each_leaf_has_one_neighbor() {
        let g = star();
        assert_eq!(g.neighbors(1).len(), 1);
        assert_eq!(g.neighbors(2).len(), 1);
        assert_eq!(g.neighbors(3).len(), 1);
    }

    #[test]
    fn index_operator_yields_same_slice_as_neighbors() {
        let g = triangle();
        for i in 0..3 {
            assert_eq!(&g[i], g.neighbors(i));
        }
    }

    #[test]
    fn adjacency_is_independent_of_edge_order() {
        let forward = AdjacencyList::build(3, [(0, 0, 1), (1, 1, 2), (2, 0, 2)]);
        let reversed = AdjacencyList::build(3, [(2, 2, 0), (1, 2, 1), (0, 1, 0)]);
        for i in 0..3 {
            let mut fa: Vec<_> = forward.neighbors(i).to_vec();
            let mut rb: Vec<_> = reversed.neighbors(i).to_vec();
            fa.sort();
            rb.sort();
            assert_eq!(fa, rb, "node {i} neighbor sets differ");
        }
    }
}

/// An indexed undirected adjacency list.
///
/// Represents an undirected graph over integer nodes `0..n`. Each undirected
/// edge supplied to [`build`](Self::build) is stored in both directions:
/// node `a` records `(edge, b)` and node `b` records `(edge, a)`.
///
/// Obtain via [`build`](Self::build).
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

use std::collections::VecDeque;

/// No target vertex: an unmapped pattern vertex, or a free target vertex.
const UNMAPPED: usize = usize::MAX;

/// Lazily enumerates the subgraph monomorphisms of `pattern` into `target`.
///
/// Both graphs are adjacency lists over the vertices `0..len`, each entry a
/// `(neighbour, edge)` pair whose `edge` indexes the caller's edge attributes. A
/// monomorphism is an injective map of every pattern vertex to a distinct target
/// vertex under which every pattern edge meets a target edge — the target may
/// carry further edges among the image, so the match is a subgraph, not an
/// induced one. `vertex_compat(p, t)` and `edge_compat(pe, te)` gate which
/// vertices and edges may correspond.
///
/// The returned iterator yields each complete mapping (`mapping[pattern_vertex] =
/// target_vertex`) in turn, so the first match, a count, or all of them each cost
/// only the search they need. Pattern vertices are matched in a connected order —
/// each after the first adjacent to one already placed — so the candidates at
/// every step are the neighbours of an image already chosen, the restriction that
/// keeps the search fast on the sparse, connected patterns chemistry poses.
///
/// # Complexity
///
/// Exponential in the worst case, as subgraph isomorphism is NP-complete;
/// near-linear in practice for connected patterns.
pub fn embeddings<VC, EC>(
    pattern: &[Vec<(usize, usize)>],
    target: Vec<Vec<(usize, usize)>>,
    vertex_compat: VC,
    edge_compat: EC,
) -> Embeddings<VC, EC>
where
    VC: Fn(usize, usize) -> bool,
    EC: Fn(usize, usize) -> bool,
{
    let (order, parents) = match_order(pattern);
    Embeddings {
        forward: vec![UNMAPPED; order.len()],
        reverse: vec![UNMAPPED; target.len()],
        target,
        vertex_compat,
        edge_compat,
        order,
        parents,
        stack: Vec::new(),
        started: false,
    }
}

/// The depth-first search over subgraph monomorphisms, one yielded per step.
///
/// Obtain via [`embeddings`].
pub struct Embeddings<VC, EC> {
    target: Vec<Vec<(usize, usize)>>,
    vertex_compat: VC,
    edge_compat: EC,
    order: Vec<usize>,
    parents: Vec<Vec<(usize, usize)>>,
    forward: Vec<usize>,
    reverse: Vec<usize>,
    stack: Vec<Frame>,
    started: bool,
}

/// One level of the search: the target vertices left to try for a pattern vertex.
struct Frame {
    candidates: Vec<usize>,
    cursor: usize,
}

impl<VC: Fn(usize, usize) -> bool, EC: Fn(usize, usize) -> bool> Iterator for Embeddings<VC, EC> {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Vec<usize>> {
        let n = self.order.len();
        if !self.started {
            self.started = true;
            if n == 0 {
                return Some(Vec::new());
            }
            if n > self.target.len() {
                return None;
            }
            let candidates = self.candidates(0);
            self.stack.push(Frame {
                candidates,
                cursor: 0,
            });
        } else if self.stack.is_empty() {
            return None;
        } else {
            self.unplace(self.order[self.stack.len() - 1]);
        }

        while !self.stack.is_empty() {
            let depth = self.stack.len() - 1;
            let p = self.order[depth];
            let mut placed = false;
            while self.stack[depth].cursor < self.stack[depth].candidates.len() {
                let t = self.stack[depth].candidates[self.stack[depth].cursor];
                self.stack[depth].cursor += 1;
                if self.reverse[t] == UNMAPPED && self.feasible(p, t, depth) {
                    self.forward[p] = t;
                    self.reverse[t] = p;
                    placed = true;
                    break;
                }
            }
            if !placed {
                self.stack.pop();
                if !self.stack.is_empty() {
                    self.unplace(self.order[self.stack.len() - 1]);
                }
                continue;
            }
            if depth + 1 == n {
                return Some(self.forward.clone());
            }
            let candidates = self.candidates(depth + 1);
            self.stack.push(Frame {
                candidates,
                cursor: 0,
            });
        }
        None
    }
}

impl<VC: Fn(usize, usize) -> bool, EC: Fn(usize, usize) -> bool> Embeddings<VC, EC> {
    /// The target vertices to try for the pattern vertex due at `depth`: the
    /// neighbours of an already-mapped neighbour's image, or — for the first
    /// vertex of a component — every vertex.
    fn candidates(&self, depth: usize) -> Vec<usize> {
        match self.parents[depth].first() {
            Some(&(parent, _)) => self.target[self.forward[parent]]
                .iter()
                .map(|&(t, _)| t)
                .collect(),
            None => (0..self.target.len()).collect(),
        }
    }

    /// Whether the pattern vertex due at `depth` may map to `t`: their colours
    /// agree, and every edge back to an already-mapped neighbour meets a
    /// compatible target edge.
    fn feasible(&self, p: usize, t: usize, depth: usize) -> bool {
        (self.vertex_compat)(p, t)
            && self.parents[depth].iter().all(|&(parent, edge)| {
                match self.edge(t, self.forward[parent]) {
                    Some(target_edge) => (self.edge_compat)(edge, target_edge),
                    None => false,
                }
            })
    }

    /// The edge joining target vertices `t` and `other`, if one exists.
    fn edge(&self, t: usize, other: usize) -> Option<usize> {
        self.target[t]
            .iter()
            .find(|&&(neighbour, _)| neighbour == other)
            .map(|&(_, edge)| edge)
    }

    /// Frees the target vertex pattern vertex `p` was mapped to.
    fn unplace(&mut self, p: usize) {
        let t = self.forward[p];
        self.forward[p] = UNMAPPED;
        self.reverse[t] = UNMAPPED;
    }
}

/// A connected order to match the pattern vertices in, with each position's
/// already-ordered neighbours.
///
/// Breadth-first from the lowest-indexed vertex of each component, so every
/// vertex after a component's first carries a neighbour earlier in the order.
/// Returns the order, and for each of its positions the `(neighbour, edge)` pairs
/// of the vertex there that precede it — the edges a candidate must honour and,
/// in the first, the image whose neighbourhood the candidates are drawn from.
fn match_order(pattern: &[Vec<(usize, usize)>]) -> (Vec<usize>, Vec<Vec<(usize, usize)>>) {
    let n = pattern.len();
    let mut order = Vec::with_capacity(n);
    let mut position = vec![UNMAPPED; n];
    for root in 0..n {
        if position[root] != UNMAPPED {
            continue;
        }
        position[root] = order.len();
        order.push(root);
        let mut queue = VecDeque::from([root]);
        while let Some(v) = queue.pop_front() {
            for &(neighbour, _) in &pattern[v] {
                if position[neighbour] == UNMAPPED {
                    position[neighbour] = order.len();
                    order.push(neighbour);
                    queue.push_back(neighbour);
                }
            }
        }
    }

    let parents = order
        .iter()
        .enumerate()
        .map(|(depth, &v)| {
            pattern[v]
                .iter()
                .copied()
                .filter(|&(neighbour, _)| position[neighbour] < depth)
                .collect()
        })
        .collect();
    (order, parents)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all(pattern: &[Vec<(usize, usize)>], target: &[Vec<(usize, usize)>]) -> Vec<Vec<usize>> {
        embeddings(pattern, target.to_vec(), |_, _| true, |_, _| true).collect()
    }

    fn single() -> Vec<Vec<(usize, usize)>> {
        vec![vec![]]
    }

    fn edge() -> Vec<Vec<(usize, usize)>> {
        vec![vec![(1, 0)], vec![(0, 0)]]
    }

    fn path() -> Vec<Vec<(usize, usize)>> {
        vec![vec![(1, 0)], vec![(0, 0), (2, 1)], vec![(1, 1)]]
    }

    fn pair() -> Vec<Vec<(usize, usize)>> {
        vec![vec![], vec![]]
    }

    fn triangle() -> Vec<Vec<(usize, usize)>> {
        vec![
            vec![(1, 0), (2, 1)],
            vec![(0, 0), (2, 2)],
            vec![(0, 1), (1, 2)],
        ]
    }

    #[test]
    fn empty_pattern_matches_once() {
        assert_eq!(all(&[], &triangle()), vec![Vec::<usize>::new()]);
    }

    #[test]
    fn pattern_larger_than_target_does_not_match() {
        assert!(all(&triangle(), &edge()).is_empty());
    }

    #[test]
    fn single_vertex_matches_every_vertex() {
        assert_eq!(all(&single(), &triangle()), vec![vec![0], vec![1], vec![2]]);
    }

    #[test]
    fn edge_matches_in_both_directions() {
        assert_eq!(all(&edge(), &triangle()).len(), 6);
    }

    #[test]
    fn path_matches_every_walk() {
        assert_eq!(all(&path(), &triangle()).len(), 6);
    }

    #[test]
    fn disconnected_pattern_matches_distinct_vertices() {
        let pairs = all(&pair(), &triangle());
        assert_eq!(pairs.len(), 6);
        assert!(pairs.iter().all(|m| m[0] != m[1]));
    }

    #[test]
    fn vertex_colours_filter_matches() {
        let pcolour = [0usize];
        let tcolour = [0usize, 1, 0];
        let found: Vec<_> = embeddings(
            &single(),
            triangle(),
            |p, t| pcolour[p] == tcolour[t],
            |_, _| true,
        )
        .collect();
        assert_eq!(found, vec![vec![0], vec![2]]);
    }

    #[test]
    fn edge_colours_filter_matches() {
        let pcolour = [0usize];
        let tcolour = [0usize, 1, 1];
        let found: Vec<_> = embeddings(
            &edge(),
            triangle(),
            |_, _| true,
            |pe, te| pcolour[pe] == tcolour[te],
        )
        .collect();
        assert_eq!(found, vec![vec![0, 1], vec![1, 0]]);
    }

    #[test]
    fn the_search_is_lazy() {
        let mut search = embeddings(&single(), triangle(), |_, _| true, |_, _| true);
        assert_eq!(search.next(), Some(vec![0]));
        assert_eq!(search.next(), Some(vec![1]));
    }
}

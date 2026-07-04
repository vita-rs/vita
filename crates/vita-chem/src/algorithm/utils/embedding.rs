use std::collections::VecDeque;

use super::AdjacencyList;

/// No target vertex: an unplaced pattern vertex, or a free target vertex.
const UNPLACED: usize = usize::MAX;

/// Lazily enumerates the subgraph monomorphisms of `pattern` into `target`.
///
/// Both graphs are [`AdjacencyList`]s over the vertices `0..len`, their edges
/// carrying the indices the caller reads its own attributes by. A monomorphism
/// maps every pattern vertex to a distinct target vertex so that every pattern
/// edge meets a target edge; the target may hold further edges among the image,
/// so the match is a subgraph, not an induced one. `vertex_compat(p, t)` gates
/// which vertices may correspond, and `edge_compat(pe, te)` which edges.
///
/// The iterator yields each mapping (`mapping[pattern_vertex] = target_vertex`)
/// in turn, so the first match, a count, or all of them each cost only the
/// search they need. Pattern vertices are placed in a connected order — each
/// after the first adjacent to one already placed — so every step draws its
/// candidates from a placed neighbour's images, the restriction that keeps the
/// search fast on the sparse, connected patterns chemistry poses.
///
/// # Complexity
///
/// O(T^P) time in the worst case and O(P · T) auxiliary space, over a pattern of
/// `P` vertices and a target of `T` vertices; near-linear in practice for
/// connected patterns, whose candidates each come from one placed image.
pub fn embeddings<VC, EC>(
    pattern: &AdjacencyList,
    target: AdjacencyList,
    vertex_compat: VC,
    edge_compat: EC,
) -> impl Iterator<Item = Vec<usize>> + use<VC, EC>
where
    VC: Fn(usize, usize) -> bool,
    EC: Fn(usize, usize) -> bool,
{
    let (order, parents) = match_order(pattern);
    Embeddings {
        forward: vec![UNPLACED; order.len()],
        reverse: vec![UNPLACED; target.len()],
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
struct Embeddings<VC, EC> {
    target: AdjacencyList,
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
        let count = self.order.len();
        if !self.started {
            self.started = true;
            if count == 0 {
                return Some(Vec::new());
            }
            if count > self.target.len() {
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
            let vertex = self.order[depth];
            let mut placed = false;
            while self.stack[depth].cursor < self.stack[depth].candidates.len() {
                let target_vertex = self.stack[depth].candidates[self.stack[depth].cursor];
                self.stack[depth].cursor += 1;
                if self.reverse[target_vertex] == UNPLACED
                    && self.feasible(vertex, target_vertex, depth)
                {
                    self.forward[vertex] = target_vertex;
                    self.reverse[target_vertex] = vertex;
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
            if depth + 1 == count {
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
    /// neighbours of an already-placed neighbour's image, or — for the first
    /// vertex of a component — every target vertex.
    fn candidates(&self, depth: usize) -> Vec<usize> {
        match self.parents[depth].first() {
            Some(&(parent, _)) => self
                .target
                .neighbors(self.forward[parent])
                .iter()
                .map(|&(_, neighbour)| neighbour)
                .collect(),
            None => (0..self.target.len()).collect(),
        }
    }

    /// Whether the pattern vertex due at `depth` may map to `target_vertex`:
    /// their colours agree, and every edge back to an already-placed neighbour
    /// meets a compatible target edge.
    fn feasible(&self, vertex: usize, target_vertex: usize, depth: usize) -> bool {
        (self.vertex_compat)(vertex, target_vertex)
            && self.parents[depth].iter().all(|&(parent, edge)| {
                match self.target_edge(target_vertex, self.forward[parent]) {
                    Some(target_edge) => (self.edge_compat)(edge, target_edge),
                    None => false,
                }
            })
    }

    /// The edge joining target vertices `from` and `to`, if one exists.
    fn target_edge(&self, from: usize, to: usize) -> Option<usize> {
        self.target
            .neighbors(from)
            .iter()
            .find(|&&(_, neighbour)| neighbour == to)
            .map(|&(edge, _)| edge)
    }

    /// Frees the target vertex the pattern vertex `vertex` was placed on.
    fn unplace(&mut self, vertex: usize) {
        let target_vertex = self.forward[vertex];
        self.forward[vertex] = UNPLACED;
        self.reverse[target_vertex] = UNPLACED;
    }
}

/// A connected order to place the pattern vertices in, with each position's
/// already-ordered neighbours.
///
/// Breadth-first from the lowest-indexed vertex of each component, so every
/// vertex after a component's first has a neighbour earlier in the order.
/// Returns the order and, for each of its positions, the `(neighbour, edge)`
/// pairs of the vertex there that precede it — the edges a candidate must
/// honour, the first also naming the placed image its candidates are drawn from.
fn match_order(pattern: &AdjacencyList) -> (Vec<usize>, Vec<Vec<(usize, usize)>>) {
    let count = pattern.len();
    let mut order = Vec::with_capacity(count);
    let mut position = vec![UNPLACED; count];
    for root in 0..count {
        if position[root] != UNPLACED {
            continue;
        }
        position[root] = order.len();
        order.push(root);
        let mut queue = VecDeque::from([root]);
        while let Some(vertex) = queue.pop_front() {
            for &(_, neighbour) in pattern.neighbors(vertex) {
                if position[neighbour] == UNPLACED {
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
        .map(|(depth, &vertex)| {
            pattern
                .neighbors(vertex)
                .iter()
                .filter(|&&(_, neighbour)| position[neighbour] < depth)
                .map(|&(edge, neighbour)| (neighbour, edge))
                .collect()
        })
        .collect();
    (order, parents)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all(pattern: &AdjacencyList, target: AdjacencyList) -> Vec<Vec<usize>> {
        embeddings(pattern, target, |_, _| true, |_, _| true).collect()
    }

    fn adjacent(graph: &AdjacencyList, from: usize, to: usize) -> bool {
        graph.neighbors(from).iter().any(|&(_, nb)| nb == to)
    }

    fn empty() -> AdjacencyList {
        AdjacencyList::build(0, [])
    }

    fn vertex() -> AdjacencyList {
        AdjacencyList::build(1, [])
    }

    fn pair() -> AdjacencyList {
        AdjacencyList::build(2, [])
    }

    fn edge() -> AdjacencyList {
        AdjacencyList::build(2, [(0, 0, 1)])
    }

    fn path() -> AdjacencyList {
        AdjacencyList::build(3, [(0, 0, 1), (1, 1, 2)])
    }

    fn triangle() -> AdjacencyList {
        AdjacencyList::build(3, [(0, 0, 1), (1, 1, 2), (2, 0, 2)])
    }

    #[test]
    fn empty_pattern_embeds_once() {
        assert_eq!(all(&empty(), triangle()), vec![Vec::<usize>::new()]);
    }

    #[test]
    fn single_vertex_embeds_at_every_target_vertex() {
        assert_eq!(all(&vertex(), triangle()), vec![vec![0], vec![1], vec![2]]);
    }

    #[test]
    fn an_edge_embeds_at_every_adjacent_pair() {
        assert_eq!(all(&edge(), triangle()).len(), 6);
    }

    #[test]
    fn pattern_larger_than_target_does_not_embed() {
        assert!(all(&triangle(), edge()).is_empty());
    }

    #[test]
    fn an_edge_does_not_embed_without_a_target_edge() {
        assert!(all(&edge(), pair()).is_empty());
    }

    #[test]
    fn incompatible_vertices_are_rejected() {
        let pattern_colour = [0usize];
        let target_colour = [0usize, 1, 0];
        let found: Vec<Vec<usize>> = embeddings(
            &vertex(),
            triangle(),
            |p, t| pattern_colour[p] == target_colour[t],
            |_, _| true,
        )
        .collect();
        assert_eq!(found, vec![vec![0], vec![2]]);
    }

    #[test]
    fn incompatible_edges_are_rejected() {
        let pattern_colour = [0usize];
        let target_colour = [0usize, 1, 1];
        let found: Vec<Vec<usize>> = embeddings(
            &edge(),
            triangle(),
            |_, _| true,
            |pe, te| pattern_colour[pe] == target_colour[te],
        )
        .collect();
        assert_eq!(found, vec![vec![0, 1], vec![1, 0]]);
    }

    #[test]
    fn a_path_embeds_along_every_adjacent_walk() {
        assert_eq!(all(&path(), triangle()).len(), 6);
    }

    #[test]
    fn a_cyclic_pattern_embeds_honouring_every_edge() {
        assert_eq!(all(&triangle(), triangle()).len(), 6);
    }

    #[test]
    fn a_disjoint_pattern_embeds_at_distinct_vertices() {
        let found = all(&pair(), triangle());
        assert_eq!(found.len(), 6);
        assert!(found.iter().all(|m| m[0] != m[1]));
    }

    #[test]
    fn every_embedding_preserves_edges() {
        let found = all(&path(), triangle());
        let target = triangle();
        for m in &found {
            assert!(adjacent(&target, m[0], m[1]));
            assert!(adjacent(&target, m[1], m[2]));
        }
    }

    #[test]
    fn the_search_is_lazy() {
        let mut search = embeddings(&vertex(), triangle(), |_, _| true, |_, _| true);
        assert_eq!(search.next(), Some(vec![0]));
        assert_eq!(search.next(), Some(vec![1]));
    }
}

/// A canonical labelling of a coloured graph, with its automorphism orbits.
pub struct Labelling {
    ranks: Vec<usize>,
    orbits: Vec<usize>,
}

impl Labelling {
    /// The canonical rank of each vertex: a permutation of `0..n`.
    pub fn ranks(&self) -> &[usize] {
        &self.ranks
    }

    /// The orbit of each vertex, named by its least member: two vertices share a
    /// value exactly when an automorphism carries one onto the other.
    pub fn orbits(&self) -> &[usize] {
        &self.orbits
    }
}

/// Canonically labels a coloured graph, with its automorphism orbits.
///
/// Takes the graph as an adjacency list over the vertices `0..adjacency.len()`,
/// each entry a `(neighbour, edge colour)` pair, together with an initial `seed`
/// colour per vertex. Returns every vertex's rank in a total order fixed by the
/// graph and the colouring, and the orbit it shares with the vertices an
/// automorphism can carry it onto — neither depending on the order the vertices
/// were given in.
///
/// Colour refinement settles the vertices into classes by their coloured
/// neighbourhoods; where symmetry leaves a class unsplit the search
/// individualises each member in turn, refines again, and keeps the labelling of
/// least certificate. Two leaves of equal certificate exhibit an automorphism,
/// which both joins orbits and prunes the symmetric branches it makes redundant.
/// Taking the least over every branch — rather than committing to a greedy
/// tie-break — is what frees the ranks from the input order.
///
/// # Complexity
///
/// O(V · (V + E) · log V) per refinement — one for a rigid graph, one per search
/// node where symmetry forces a branch. Near-linear in practice, exponential in
/// the worst case.
pub fn labelling(adjacency: &[Vec<(usize, usize)>], seed: &[usize]) -> Labelling {
    let n = adjacency.len();
    let mut search = Search {
        adjacency,
        seed,
        n,
        first: None,
        best: None,
        parent: (0..n).collect(),
        generators: Vec::new(),
    };
    search.descend(seed.to_vec(), &[]);
    let ranks = search.best.expect("the search reaches a leaf").1;
    let orbits = (0..n)
        .map(|vertex| find(&mut search.parent, vertex))
        .collect();
    Labelling { ranks, orbits }
}

/// The individualisation–refinement search for the least certificate and the
/// automorphisms found along the way.
struct Search<'a> {
    adjacency: &'a [Vec<(usize, usize)>],
    seed: &'a [usize],
    n: usize,
    first: Option<(Vec<usize>, Vec<usize>)>,
    best: Option<(Vec<usize>, Vec<usize>)>,
    parent: Vec<usize>,
    generators: Vec<Vec<usize>>,
}

impl Search<'_> {
    /// Refines, then records a discrete colouring as a leaf or individualises a
    /// target cell vertex by vertex, skipping any an automorphism maps onto one
    /// already taken.
    fn descend(&mut self, mut colours: Vec<usize>, path: &[usize]) {
        let count = refine(self.adjacency, &mut colours);
        if count == self.n {
            self.leaf(colours);
            return;
        }

        let mut quotient: Vec<usize> = (0..self.n).collect();
        for generator in &self.generators {
            if path.iter().all(|&fixed| generator[fixed] == fixed) {
                for (vertex, &image) in generator.iter().enumerate() {
                    union(&mut quotient, vertex, image);
                }
            }
        }

        let mut taken = vec![false; self.n];
        for vertex in self.target(&colours, count) {
            let orbit = find(&mut quotient, vertex);
            if taken[orbit] {
                continue;
            }
            taken[orbit] = true;
            let mut next: Vec<usize> = colours.iter().map(|&colour| colour * 2).collect();
            next[vertex] += 1;
            let mut child = path.to_vec();
            child.push(vertex);
            self.descend(next, &child);
        }
    }

    /// Records a discrete colouring: the first leaf sets the reference, a leaf of
    /// equal certificate yields an automorphism, and the least certificate keeps
    /// the canonical labelling.
    fn leaf(&mut self, colours: Vec<usize>) {
        let certificate = self.certificate(&colours);
        let generator = match &self.first {
            None => {
                self.best = Some((certificate.clone(), colours.clone()));
                self.first = Some((certificate, colours));
                return;
            }
            Some((first_certificate, first)) => {
                (certificate == *first_certificate).then(|| automorphism(first, &colours))
            }
        };
        if let Some(generator) = generator {
            for (vertex, &image) in generator.iter().enumerate() {
                union(&mut self.parent, vertex, image);
            }
            self.generators.push(generator);
        }
        if certificate < self.best.as_ref().expect("the first leaf set the best").0 {
            self.best = Some((certificate, colours));
        }
    }

    /// The vertices of the smallest non-singleton colour class. Ties fall to the
    /// lowest colour, itself canonical, so the branch taken — and the least
    /// certificate it leads to — does not depend on input order.
    fn target(&self, colours: &[usize], count: usize) -> Vec<usize> {
        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); count];
        for (vertex, &colour) in colours.iter().enumerate() {
            cells[colour].push(vertex);
        }
        cells
            .into_iter()
            .filter(|cell| cell.len() > 1)
            .min_by_key(Vec::len)
            .expect("a non-discrete colouring has a non-singleton class")
    }

    /// The certificate of a discrete colouring: the seeded, labelled graph
    /// written out in rank order. The least certificate over all leaves picks the
    /// canonical labelling.
    fn certificate(&self, colours: &[usize]) -> Vec<usize> {
        let mut order = vec![0; self.n];
        for (vertex, &rank) in colours.iter().enumerate() {
            order[rank] = vertex;
        }
        let mut certificate = Vec::new();
        for &vertex in &order {
            certificate.push(self.seed[vertex]);
            let mut incident: Vec<(usize, usize)> = self.adjacency[vertex]
                .iter()
                .map(|&(neighbour, edge)| (edge, colours[neighbour]))
                .collect();
            incident.sort_unstable();
            certificate.push(incident.len());
            for (edge, neighbour) in incident {
                certificate.push(edge);
                certificate.push(neighbour);
            }
        }
        certificate
    }
}

/// Refines a vertex colouring to the coarsest equitable partition.
///
/// Replaces each vertex's colour, round on round, by its own colour paired with
/// the sorted multiset of its incident `(edge colour, neighbour colour)`, until
/// no class splits further; returns the number of resulting classes. The refined
/// colours are dense and canonical — they run `0..count`, ordered by the
/// signature that produced them — so the partition depends only on the initial
/// colouring and the graph, never on vertex order. A finer initial colouring is
/// honoured: classes only split, never merge.
fn refine(adjacency: &[Vec<(usize, usize)>], colours: &mut [usize]) -> usize {
    let n = adjacency.len();
    if n == 0 {
        return 0;
    }

    let mut count = 0;
    loop {
        let signatures: Vec<(usize, Vec<(usize, usize)>)> = (0..n)
            .map(|vertex| {
                let mut incident: Vec<(usize, usize)> = adjacency[vertex]
                    .iter()
                    .map(|&(neighbour, edge)| (edge, colours[neighbour]))
                    .collect();
                incident.sort_unstable();
                (colours[vertex], incident)
            })
            .collect();

        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by(|&a, &b| signatures[a].cmp(&signatures[b]));
        let mut rank = 0;
        for i in 0..n {
            if i > 0 && signatures[order[i]] != signatures[order[i - 1]] {
                rank += 1;
            }
            colours[order[i]] = rank;
        }

        let split = rank + 1;
        if split == count {
            return count;
        }
        count = split;
    }
}

/// The vertex permutation carrying the labelling `leaf` onto `first`: an
/// automorphism, as both label the graph into the same certificate.
fn automorphism(first: &[usize], leaf: &[usize]) -> Vec<usize> {
    let mut inverse = vec![0; first.len()];
    for (vertex, &rank) in first.iter().enumerate() {
        inverse[rank] = vertex;
    }
    leaf.iter().map(|&rank| inverse[rank]).collect()
}

/// Returns the representative of `vertex`'s set, halving the path to it.
fn find(parent: &mut [usize], mut vertex: usize) -> usize {
    while parent[vertex] != vertex {
        parent[vertex] = parent[parent[vertex]];
        vertex = parent[vertex];
    }
    vertex
}

/// Joins the sets of `a` and `b` under the lesser representative.
fn union(parent: &mut [usize], a: usize, b: usize) {
    let (a, b) = (find(parent, a), find(parent, b));
    if a != b {
        parent[a.max(b)] = a.min(b);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ranks(adjacency: &[Vec<(usize, usize)>], seed: &[usize]) -> Vec<usize> {
        labelling(adjacency, seed).ranks
    }

    fn orbits(adjacency: &[Vec<(usize, usize)>], seed: &[usize]) -> Vec<usize> {
        labelling(adjacency, seed).orbits
    }

    #[test]
    fn empty_graph_has_no_labelling() {
        let labelling = labelling(&[], &[]);
        assert!(labelling.ranks.is_empty());
        assert!(labelling.orbits.is_empty());
    }

    #[test]
    fn single_vertex_ranks_zero_in_its_own_orbit() {
        let labelling = labelling(&[vec![]], &[0]);
        assert_eq!(labelling.ranks, [0]);
        assert_eq!(labelling.orbits, [0]);
    }

    #[test]
    fn ranks_are_a_permutation() {
        let path = [vec![(1, 0)], vec![(0, 0), (2, 0)], vec![(1, 0)]];
        let mut sorted = ranks(&path, &[0, 0, 0]);
        sorted.sort_unstable();
        assert_eq!(sorted, [0, 1, 2]);
    }

    #[test]
    fn path_ends_share_an_orbit_apart_from_the_centre() {
        let path = [vec![(1, 0)], vec![(0, 0), (2, 0)], vec![(1, 0)]];
        let orbits = orbits(&path, &[0, 0, 0]);
        assert_eq!(orbits[0], orbits[2]);
        assert_ne!(orbits[0], orbits[1]);
    }

    #[test]
    fn cycle_is_one_orbit() {
        let triangle = [
            vec![(1, 0), (2, 0)],
            vec![(0, 0), (2, 0)],
            vec![(0, 0), (1, 0)],
        ];
        let orbits = orbits(&triangle, &[0, 0, 0]);
        assert!(orbits.iter().all(|&orbit| orbit == orbits[0]));
    }

    #[test]
    fn rigid_graph_has_singleton_orbits() {
        let path = [vec![(1, 0)], vec![(0, 0), (2, 0)], vec![(1, 0)]];
        let orbits = orbits(&path, &[1, 0, 2]);
        assert_eq!(orbits, [0, 1, 2]);
    }

    #[test]
    fn seed_breaks_an_otherwise_symmetric_ranking() {
        let path = [vec![(1, 0)], vec![(0, 0), (2, 0)], vec![(1, 0)]];
        assert_eq!(ranks(&path, &[1, 0, 2]), [1, 0, 2]);
    }

    #[test]
    fn edge_colours_distinguish_neighbours() {
        let path = [vec![(1, 1)], vec![(0, 1), (2, 0)], vec![(1, 0)]];
        assert_eq!(orbits(&path, &[0, 0, 0]), [0, 1, 2]);
    }

    #[test]
    fn labelling_is_independent_of_vertex_order() {
        let edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)];
        let build = |perm: &[usize; 4]| {
            let mut adjacency = vec![Vec::new(); 4];
            for &(a, b) in &edges {
                adjacency[perm[a]].push((perm[b], 0));
                adjacency[perm[b]].push((perm[a], 0));
            }
            labelling(&adjacency, &[0; 4])
        };
        let plain = build(&[0, 1, 2, 3]);
        let shuffled = build(&[2, 0, 3, 1]);
        for vertex in 0..4 {
            assert_eq!(plain.ranks[vertex], shuffled.ranks[[2, 0, 3, 1][vertex]]);
            assert_eq!(
                plain.orbits[vertex] == plain.orbits[0],
                shuffled.orbits[[2, 0, 3, 1][vertex]] == shuffled.orbits[2],
            );
        }
    }

    #[test]
    fn complete_graph_is_one_orbit() {
        let n = 6;
        let complete: Vec<Vec<(usize, usize)>> = (0..n)
            .map(|v| (0..n).filter(|&u| u != v).map(|u| (u, 0)).collect())
            .collect();
        let orbits = orbits(&complete, &[0; 6]);
        assert!(orbits.iter().all(|&orbit| orbit == orbits[0]));
    }

    #[test]
    fn cube_is_one_orbit() {
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        let mut cube = vec![Vec::new(); 8];
        for &(a, b) in &edges {
            cube[a].push((b, 0));
            cube[b].push((a, 0));
        }
        let orbits = orbits(&cube, &[0; 8]);
        assert!(orbits.iter().all(|&orbit| orbit == orbits[0]));
    }
}

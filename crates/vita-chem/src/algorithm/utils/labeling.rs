/// A canonical labeling of a colored graph, paired with its automorphism orbits.
///
/// [`ranks`](Self::ranks) place the vertices in a total order fixed by the graph
/// and its coloring; [`orbits`](Self::orbits) name each vertex's symmetry class.
/// Neither depends on the order the vertices were given in.
///
/// Obtain via [`labeling`].
pub struct Labeling {
    ranks: Vec<usize>,
    orbits: Vec<usize>,
}

impl Labeling {
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

/// Canonically labels a colored graph, with its automorphism orbits.
///
/// The graph is an adjacency list over the vertices `0..adjacency.len()`, each
/// entry a `(neighbor, edge color)` pair, together with an initial `seed` color
/// per vertex. Returns every vertex's rank in a total order fixed by the graph
/// and the coloring, and the orbit it shares with the vertices an automorphism
/// can carry it onto — neither depending on the order the vertices were given in.
///
/// Color refinement settles the vertices into classes by their colored
/// neighborhoods; where symmetry leaves a class unsplit, the search
/// individualizes each member in turn, refines again, and keeps the labeling of
/// least certificate. Two leaves of equal certificate exhibit an automorphism,
/// which both joins orbits and prunes the symmetric branches it makes redundant.
/// Taking the least over every branch — rather than a greedy tie-break — is what
/// frees the labeling from the input order.
///
/// # Complexity
///
/// O(V · (V + E) · log V) per refinement — one for a rigid graph, one per search
/// node where symmetry forces a branch. Near-linear in practice, exponential in
/// the worst case.
pub fn labeling(adjacency: &[Vec<(usize, usize)>], seed: &[usize]) -> Labeling {
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
    Labeling { ranks, orbits }
}

/// The individualization–refinement search for the least certificate and the
/// automorphisms met along the way.
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
    /// Refines, then records a discrete coloring as a leaf, or individualizes a
    /// target-cell vertex by vertex, skipping any an automorphism maps onto one
    /// already taken.
    fn descend(&mut self, mut colors: Vec<usize>, path: &[usize]) {
        let count = refine(self.adjacency, &mut colors);
        if count == self.n {
            self.leaf(colors);
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
        for vertex in self.target(&colors, count) {
            let orbit = find(&mut quotient, vertex);
            if taken[orbit] {
                continue;
            }
            taken[orbit] = true;
            let mut refined: Vec<usize> = colors.iter().map(|&color| color * 2).collect();
            refined[vertex] += 1;
            let mut child = path.to_vec();
            child.push(vertex);
            self.descend(refined, &child);
        }
    }

    /// Records a discrete coloring: the first leaf sets the reference, a leaf of
    /// equal certificate yields an automorphism, and the least certificate keeps
    /// the canonical labeling.
    fn leaf(&mut self, ranks: Vec<usize>) {
        let certificate = self.certificate(&ranks);
        let generator = match &self.first {
            None => {
                self.best = Some((certificate.clone(), ranks.clone()));
                self.first = Some((certificate, ranks));
                return;
            }
            Some((first_certificate, first)) => {
                (certificate == *first_certificate).then(|| automorphism(first, &ranks))
            }
        };
        if let Some(generator) = generator {
            for (vertex, &image) in generator.iter().enumerate() {
                union(&mut self.parent, vertex, image);
            }
            self.generators.push(generator);
        }
        if certificate < self.best.as_ref().expect("the first leaf set the best").0 {
            self.best = Some((certificate, ranks));
        }
    }

    /// The vertices of the smallest non-singleton color class. Ties fall to the
    /// lowest color, itself canonical, so the branch taken — and the least
    /// certificate it leads to — does not depend on input order.
    fn target(&self, colors: &[usize], count: usize) -> Vec<usize> {
        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); count];
        for (vertex, &color) in colors.iter().enumerate() {
            cells[color].push(vertex);
        }
        cells
            .into_iter()
            .filter(|cell| cell.len() > 1)
            .min_by_key(Vec::len)
            .expect("a non-discrete coloring has a non-singleton class")
    }

    /// The certificate of a discrete coloring: the seeded, labeled graph written
    /// out in rank order. The least certificate over all leaves picks the
    /// canonical labeling.
    fn certificate(&self, ranks: &[usize]) -> Vec<usize> {
        let mut order = vec![0; self.n];
        for (vertex, &rank) in ranks.iter().enumerate() {
            order[rank] = vertex;
        }
        let mut certificate = Vec::new();
        for &vertex in &order {
            certificate.push(self.seed[vertex]);
            let mut incident: Vec<(usize, usize)> = self.adjacency[vertex]
                .iter()
                .map(|&(neighbor, edge)| (edge, ranks[neighbor]))
                .collect();
            incident.sort_unstable();
            certificate.push(incident.len());
            for (edge, neighbor) in incident {
                certificate.push(edge);
                certificate.push(neighbor);
            }
        }
        certificate
    }
}

/// Refines a vertex coloring to the coarsest equitable partition.
///
/// Replaces each vertex's color, round on round, by its own color paired with
/// the sorted multiset of its incident `(edge color, neighbor color)`, until no
/// class splits further; returns the number of resulting classes. The refined
/// colors are dense and canonical — they run `0..count`, ordered by the signature
/// that produced them — so the partition depends only on the initial coloring and
/// the graph, never on vertex order. A finer initial coloring is honored: classes
/// only split, never merge.
fn refine(adjacency: &[Vec<(usize, usize)>], colors: &mut [usize]) -> usize {
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
                    .map(|&(neighbor, edge)| (edge, colors[neighbor]))
                    .collect();
                incident.sort_unstable();
                (colors[vertex], incident)
            })
            .collect();

        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by(|&a, &b| signatures[a].cmp(&signatures[b]));
        let mut rank = 0;
        for i in 0..n {
            if i > 0 && signatures[order[i]] != signatures[order[i - 1]] {
                rank += 1;
            }
            colors[order[i]] = rank;
        }

        let split = rank + 1;
        if split == count {
            return count;
        }
        count = split;
    }
}

/// The vertex permutation carrying the labeling `leaf` onto `first`: an
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

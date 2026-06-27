/// Refines a vertex colouring to the coarsest equitable partition.
///
/// Takes the graph as an adjacency list over the vertices `0..adjacency.len()`,
/// each entry a `(neighbour, edge colour)` pair, and refines the initial
/// `colours` in place by 1-dimensional Weisfeiler–Leman: a vertex's colour is
/// replaced, round on round, by its own colour paired with the sorted multiset
/// of its incident `(edge colour, neighbour colour)`, until no class splits
/// further. Returns the number of resulting classes.
///
/// The refined colours are dense and canonical — they run `0..count`, ordered by
/// the signature that produced them — so both the partition and the colour
/// values depend only on the initial colouring and the graph, never on vertex
/// order. A finer initial colouring is honoured: classes only split, never merge.
///
/// # Complexity
///
/// O(V · (V + E) · log V) time, O(V + E) space.
pub(crate) fn refine(adjacency: &[Vec<(usize, usize)>], colours: &mut [usize]) -> usize {
    let n = adjacency.len();
    if n == 0 {
        return 0;
    }

    let mut count = 0;
    loop {
        let signatures: Vec<(usize, Vec<(usize, usize)>)> = (0..n)
            .map(|v| {
                let mut incident: Vec<(usize, usize)> =
                    adjacency[v].iter().map(|&(u, e)| (e, colours[u])).collect();
                incident.sort_unstable();
                (colours[v], incident)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_graph_has_no_classes() {
        assert_eq!(refine(&[], &mut []), 0);
    }

    #[test]
    fn single_vertex_is_one_class() {
        let mut colours = [0];
        assert_eq!(refine(&[vec![]], &mut colours), 1);
        assert_eq!(colours, [0]);
    }

    #[test]
    fn path_splits_ends_from_middle() {
        let adjacency = vec![vec![(1, 0)], vec![(0, 0), (2, 0)], vec![(1, 0)]];
        let mut colours = [0, 0, 0];
        assert_eq!(refine(&adjacency, &mut colours), 2);
        assert_eq!(colours[0], colours[2]);
        assert_ne!(colours[0], colours[1]);
    }

    #[test]
    fn cycle_is_one_class() {
        let adjacency = vec![
            vec![(1, 0), (2, 0)],
            vec![(0, 0), (2, 0)],
            vec![(0, 0), (1, 0)],
        ];
        let mut colours = [0, 0, 0];
        assert_eq!(refine(&adjacency, &mut colours), 1);
    }

    #[test]
    fn initial_colours_are_honoured() {
        let adjacency = vec![vec![], vec![]];
        let mut apart = [0, 1];
        assert_eq!(refine(&adjacency, &mut apart), 2);
        let mut alike = [0, 0];
        assert_eq!(refine(&adjacency, &mut alike), 1);
    }

    #[test]
    fn edge_colours_distinguish_neighbours() {
        let adjacency = vec![vec![(1, 1)], vec![(0, 1), (2, 0)], vec![(1, 0)]];
        let mut colours = [0, 0, 0];
        assert_eq!(refine(&adjacency, &mut colours), 3);
    }

    #[test]
    fn colours_are_dense() {
        let adjacency = vec![vec![(1, 0)], vec![(0, 0), (2, 0)], vec![(1, 0)]];
        let mut colours = [0, 0, 0];
        let count = refine(&adjacency, &mut colours);
        assert!(colours.iter().all(|&c| c < count));
        assert_eq!(colours.iter().copied().max().unwrap(), count - 1);
    }

    #[test]
    fn refinement_is_independent_of_vertex_order() {
        let edges = [(0, 1), (1, 2), (2, 3), (1, 4)];
        let build = |perm: &[usize; 5]| {
            let mut adjacency = vec![Vec::new(); 5];
            for &(a, b) in &edges {
                adjacency[perm[a]].push((perm[b], 0));
                adjacency[perm[b]].push((perm[a], 0));
            }
            let mut colours = vec![0; 5];
            refine(&adjacency, &mut colours);
            colours
        };

        let identity = [0, 1, 2, 3, 4];
        let shuffled = [3, 0, 4, 2, 1];
        let plain = build(&identity);
        let relabelled = build(&shuffled);
        for v in 0..5 {
            assert_eq!(plain[v], relabelled[shuffled[v]]);
        }
    }
}

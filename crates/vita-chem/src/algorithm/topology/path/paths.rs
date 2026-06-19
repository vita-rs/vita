use std::collections::{HashMap, VecDeque};

use vita_core::SiteId;

use crate::HasBonds;

/// All shortest topological paths between two sites, each inclusive of both endpoints.
///
/// Uses BFS to determine the minimum distance, then collects every path of
/// that length. Returns an empty vector when `start` and `end` lie in
/// different connected components. Returns `vec![[start]]` when `start == end`.
/// The order of paths in the returned vector is unspecified.
///
/// # Complexity
///
/// O(V + E + P) time and O(V + P) auxiliary space, where P is the total
/// number of sites across all returned paths.
pub fn paths<M: HasBonds>(mol: &M, start: SiteId, end: SiteId) -> Vec<Vec<SiteId>> {
    if start == end {
        return vec![vec![start]];
    }

    let mut dist: HashMap<SiteId, usize> = HashMap::new();
    let mut queue: VecDeque<SiteId> = VecDeque::new();

    dist.insert(start, 0);
    queue.push_back(start);

    while let Some(site) = queue.pop_front() {
        let d = dist[&site];
        for nb in mol.neighbors(site) {
            if !dist.contains_key(&nb) {
                dist.insert(nb, d + 1);
                queue.push_back(nb);
            }
        }
    }

    if !dist.contains_key(&end) {
        return vec![];
    }

    let mut result: Vec<Vec<SiteId>> = Vec::new();
    let mut stack: Vec<(SiteId, Vec<SiteId>)> = vec![(end, vec![end])];

    while let Some((site, partial)) = stack.pop() {
        if site == start {
            let mut path = partial;
            path.reverse();
            result.push(path);
            continue;
        }
        let d = dist[&site];
        for nb in mol.neighbors(site) {
            if dist.get(&nb) == Some(&(d - 1)) {
                let mut next = partial.clone();
                next.push(nb);
                stack.push((nb, next));
            }
        }
    }

    result
}

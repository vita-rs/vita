use std::collections::{HashMap, VecDeque};

use vita_core::SiteId;

use crate::HasBonds;

/// Shortest topological path between two sites, inclusive of both endpoints.
///
/// Uses BFS, guaranteeing the minimum number of bonds. The returned sequence
/// starts at `start` and ends at `end`. Returns `None` when `start` and `end`
/// lie in different connected components. Returns `Some` containing only
/// `start` when `start == end`.
///
/// # Complexity
///
/// O(V + E) time, O(V) auxiliary space.
pub fn path<M: HasBonds>(mol: &M, start: SiteId, end: SiteId) -> Option<Vec<SiteId>> {
    if start == end {
        return Some(vec![start]);
    }

    let mut parent: HashMap<SiteId, SiteId> = HashMap::new();
    let mut queue: VecDeque<SiteId> = VecDeque::new();

    parent.insert(start, start);
    queue.push_back(start);

    while let Some(site) = queue.pop_front() {
        for nb in mol.neighbors(site) {
            if parent.contains_key(&nb) {
                continue;
            }
            parent.insert(nb, site);
            if nb == end {
                let mut result = vec![end];
                let mut cur = end;
                while cur != start {
                    cur = parent[&cur];
                    result.push(cur);
                }
                result.reverse();
                return Some(result);
            }
            queue.push_back(nb);
        }
    }

    None
}

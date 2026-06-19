use std::collections::{HashSet, VecDeque};

use vita_core::SiteId;

use crate::HasBonds;

/// Breadth-first traversal order from `start`.
///
/// Yields every site reachable from `start` in breadth-first (level) order,
/// visiting sites closer to `start` before sites farther away. Sites in
/// disconnected components are never yielded. `start` is always the first
/// site yielded.
///
/// # Complexity
///
/// O(V + E) time, O(V) auxiliary space.
pub fn bfs<M: HasBonds>(mol: &M, start: SiteId) -> impl Iterator<Item = SiteId> + '_ {
    BfsIter::new(mol, start)
}

struct BfsIter<'a, M> {
    mol: &'a M,
    queue: VecDeque<SiteId>,
    visited: HashSet<SiteId>,
}

impl<'a, M: HasBonds> BfsIter<'a, M> {
    fn new(mol: &'a M, start: SiteId) -> Self {
        let mut visited = HashSet::new();
        visited.insert(start);
        let mut queue = VecDeque::new();
        queue.push_back(start);
        Self {
            mol,
            queue,
            visited,
        }
    }
}

impl<'a, M: HasBonds> Iterator for BfsIter<'a, M> {
    type Item = SiteId;

    fn next(&mut self) -> Option<SiteId> {
        let site = self.queue.pop_front()?;
        let mol = self.mol;
        for nb in mol.neighbors(site) {
            if self.visited.insert(nb) {
                self.queue.push_back(nb);
            }
        }
        Some(site)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.queue.len(), None)
    }
}

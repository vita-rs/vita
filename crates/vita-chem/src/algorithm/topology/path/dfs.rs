use std::collections::HashSet;

use vita_core::SiteId;

use crate::HasBonds;

/// Depth-first traversal order from `start`.
///
/// Yields every site reachable from `start` in depth-first order, following
/// each path to its terminus before backtracking. Sites in disconnected
/// components are never yielded. `start` is always the first site yielded.
///
/// # Complexity
///
/// O(V + E) time, O(V) auxiliary space.
pub fn dfs<M: HasBonds>(mol: &M, start: SiteId) -> impl Iterator<Item = SiteId> + '_ {
    DfsIter::new(mol, start)
}

struct DfsIter<'a, M> {
    mol: &'a M,
    stack: Vec<SiteId>,
    visited: HashSet<SiteId>,
}

impl<'a, M: HasBonds> DfsIter<'a, M> {
    fn new(mol: &'a M, start: SiteId) -> Self {
        let mut visited = HashSet::new();
        visited.insert(start);
        Self {
            mol,
            stack: vec![start],
            visited,
        }
    }
}

impl<'a, M: HasBonds> Iterator for DfsIter<'a, M> {
    type Item = SiteId;

    fn next(&mut self) -> Option<SiteId> {
        let site = self.stack.pop()?;
        let mol = self.mol;
        for nb in mol.neighbors(site) {
            if self.visited.insert(nb) {
                self.stack.push(nb);
            }
        }
        Some(site)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.stack.len(), None)
    }
}

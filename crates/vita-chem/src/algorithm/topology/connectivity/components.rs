use std::collections::HashMap;

use vita_core::{HasSites, SiteId};

use crate::HasBonds;

/// The connected components of a molecule.
///
/// Each component is a maximal set of sites that are mutually reachable
/// through bonds. Sites with no bonds form singleton components. An empty
/// molecule has no components.
///
/// Obtain via [`components`].
pub struct Components {
    groups: Vec<Vec<SiteId>>,
    index: HashMap<SiteId, usize>,
}

impl Components {
    /// Number of connected components.
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Returns `true` if the molecule contains no sites.
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// Returns `true` if the molecule is a single connected component.
    pub fn is_connected(&self) -> bool {
        self.groups.len() == 1
    }

    /// Iterates all components, each as a slice of site identifiers.
    pub fn iter(&self) -> impl Iterator<Item = &[SiteId]> + '_ {
        self.groups.iter().map(|g| g.as_slice())
    }

    /// Returns the component containing `site`.
    ///
    /// Returns `None` if `site` is not present in the molecule.
    pub fn get(&self, site: SiteId) -> Option<&[SiteId]> {
        let &g = self.index.get(&site)?;
        Some(&self.groups[g])
    }

    /// Returns `true` if `a` and `b` belong to the same connected component.
    ///
    /// Returns `false` if either site is absent from the molecule.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.index.get(&a), self.index.get(&b)) {
            (Some(ia), Some(ib)) => ia == ib,
            _ => false,
        }
    }
}

/// Connected components of a molecule.
///
/// Returns every maximal set of mutually reachable sites. Sites with no bonds
/// form singleton components. The order of components follows `mol.sites()`;
/// the order of sites within each component is DFS discovery order.
///
/// # Complexity
///
/// O(V + E) time and space.
pub fn components<M: HasBonds + HasSites>(mol: &M) -> Components {
    let mut index: HashMap<SiteId, usize> = HashMap::new();
    let mut groups: Vec<Vec<SiteId>> = Vec::new();

    for start in mol.sites() {
        if index.contains_key(&start) {
            continue;
        }

        let g = groups.len();
        let mut group = Vec::new();
        let mut stack = vec![start];
        index.insert(start, g);

        while let Some(site) = stack.pop() {
            group.push(site);
            for nb in mol.neighbors(site) {
                if let std::collections::hash_map::Entry::Vacant(e) = index.entry(nb) {
                    e.insert(g);
                    stack.push(nb);
                }
            }
        }

        groups.push(group);
    }

    Components { groups, index }
}

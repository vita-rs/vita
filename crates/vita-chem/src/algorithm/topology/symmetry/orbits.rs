use std::collections::HashMap;

use vita_core::{HasSites, SiteId};

use crate::HasBonds;

/// The symmetry classes of a molecule.
///
/// Each class is a maximal set of interchangeable sites: a symmetry of the
/// molecular graph — a relabelling that leaves it unchanged — maps any member
/// of a class onto any other. An empty molecule has no classes.
///
/// Obtain via [`orbits`].
pub struct Orbits {
    groups: Vec<Vec<SiteId>>,
    index: HashMap<SiteId, usize>,
}

impl Orbits {
    /// Number of symmetry classes.
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Returns `true` if the molecule contains no sites.
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// Iterates all classes, each as a slice of site identifiers.
    pub fn iter(&self) -> impl Iterator<Item = &[SiteId]> + '_ {
        self.groups.iter().map(|g| g.as_slice())
    }

    /// Returns the class containing `site`.
    ///
    /// Returns `None` if `site` is not present in the molecule.
    pub fn get(&self, site: SiteId) -> Option<&[SiteId]> {
        let &g = self.index.get(&site)?;
        Some(&self.groups[g])
    }

    /// Returns `true` if `a` and `b` belong to the same symmetry class.
    ///
    /// Returns `false` if either site is absent from the molecule.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.index.get(&a), self.index.get(&b)) {
            (Some(ia), Some(ib)) => ia == ib,
            _ => false,
        }
    }
}

/// Symmetry classes of a molecule.
///
/// Refines the sites by 1-dimensional Weisfeiler–Leman colouring until no class
/// splits further: two sites stay together only while their neighbours' classes
/// match. Classes are ordered by their sites, which within each class are
/// ascending.
///
/// # Complexity
///
/// O(V · (V + E)) time.
pub fn orbits<M: HasBonds + HasSites>(mol: &M) -> Orbits {
    let sites: Vec<SiteId> = mol.sites().collect();
    let n = sites.len();
    if n == 0 {
        return Orbits {
            groups: Vec::new(),
            index: HashMap::new(),
        };
    }

    let pos: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for bond in mol.bonds() {
        let (a, b) = mol.bond_endpoints(bond);
        adj[pos[&a]].push(pos[&b]);
        adj[pos[&b]].push(pos[&a]);
    }

    let mut class = vec![0usize; n];
    let mut count = 1;
    loop {
        let mut ids: HashMap<(usize, Vec<usize>), usize> = HashMap::new();
        let mut next = vec![0usize; n];
        for v in 0..n {
            let mut neighbours: Vec<usize> = adj[v].iter().map(|&u| class[u]).collect();
            neighbours.sort_unstable();
            let id = ids.len();
            next[v] = *ids.entry((class[v], neighbours)).or_insert(id);
        }
        if ids.len() == count {
            break;
        }
        count = ids.len();
        class = next;
    }

    let mut groups: Vec<Vec<SiteId>> = vec![Vec::new(); count];
    for (v, &site) in sites.iter().enumerate() {
        groups[class[v]].push(site);
    }
    for group in &mut groups {
        group.sort_unstable();
    }
    groups.sort_unstable();

    let mut index: HashMap<SiteId, usize> = HashMap::new();
    for (i, group) in groups.iter().enumerate() {
        for &site in group {
            index.insert(site, i);
        }
    }

    Orbits { groups, index }
}

use std::collections::HashMap;

use vita_core::{HasSites, SiteId};

use crate::{BondId, HasBonds};

/// A biconnected component of a molecule.
///
/// A biconnected component is a maximal set of sites and bonds that remains
/// connected after removing any single site. In a molecular graph, each block
/// is either a ring or a single bridge bond.
///
/// Obtain via [`Blocks`].
pub struct Block {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
    is_ring: bool,
}

impl Block {
    /// Sites in this block.
    pub fn sites(&self) -> &[SiteId] {
        &self.sites
    }

    /// Bonds in this block.
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }

    /// Returns `true` if the block contains a cycle, `false` if it is a
    /// single bridge bond.
    pub fn is_ring(&self) -> bool {
        self.is_ring
    }
}

/// The biconnected components (blocks) of a molecule.
///
/// Decomposes the molecular graph into maximal 2-connected subgraphs. Each
/// block is either a ring or a single bridge bond. Isolated sites (with no
/// bonds) belong to no block.
///
/// Obtain via [`blocks`].
pub struct Blocks {
    blocks: Vec<Block>,
    site_index: HashMap<SiteId, Vec<usize>>,
    bond_index: HashMap<BondId, usize>,
}

impl Blocks {
    /// Number of biconnected components.
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    /// Returns `true` if the molecule has no bonds.
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// Returns `true` if the molecule is a single biconnected component.
    pub fn is_biconnected(&self) -> bool {
        self.blocks.len() == 1
    }

    /// Iterates all biconnected components.
    pub fn iter(&self) -> impl Iterator<Item = &Block> + '_ {
        self.blocks.iter()
    }

    /// Iterates articulation points (cut sites) of the molecule.
    ///
    /// A site is an articulation point if removing it would increase the number
    /// of connected components. Equivalently, it appears in two or more blocks.
    pub fn cuts(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.site_index
            .iter()
            .filter(|(_, v)| v.len() >= 2)
            .map(|(&s, _)| s)
    }

    /// Iterates bridge bonds of the molecule.
    ///
    /// A bond is a bridge if removing it would increase the number of connected
    /// components. Every bridge forms its own single-bond block
    /// (`is_ring = false`).
    pub fn bridges(&self) -> impl Iterator<Item = BondId> + '_ {
        self.blocks
            .iter()
            .filter(|b| !b.is_ring())
            .flat_map(|b| b.bonds().iter().copied())
    }

    /// Returns `true` if `site` is an articulation point.
    ///
    /// Returns `false` if `site` is absent from the molecule or not a cut.
    pub fn contains_cut(&self, site: SiteId) -> bool {
        self.site_index.get(&site).is_some_and(|v| v.len() >= 2)
    }

    /// Returns `true` if `bond` is a bridge.
    ///
    /// Returns `false` if `bond` is absent from the molecule or not a bridge.
    pub fn contains_bridge(&self, bond: BondId) -> bool {
        self.bond_index
            .get(&bond)
            .is_some_and(|&i| !self.blocks[i].is_ring())
    }

    /// Iterates all blocks containing `site`.
    ///
    /// Returns an empty iterator if `site` is absent from the molecule or is
    /// isolated (has no bonds).
    pub fn of_site(&self, site: SiteId) -> impl Iterator<Item = &Block> + '_ {
        self.site_index
            .get(&site)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
            .iter()
            .map(|&i| &self.blocks[i])
    }

    /// Returns the block containing `bond`.
    ///
    /// Returns `None` if `bond` is not present in the molecule.
    pub fn of_bond(&self, bond: BondId) -> Option<&Block> {
        let &i = self.bond_index.get(&bond)?;
        Some(&self.blocks[i])
    }
}

/// Biconnected components of a molecule.
///
/// Decomposes the molecular graph into maximal 2-connected subgraphs using
/// Tarjan's algorithm. Each block is either a ring or a single bridge bond.
/// Isolated sites belong to no block. The order of blocks follows DFS
/// discovery order.
///
/// # Complexity
///
/// O(V + E) time and space.
pub fn blocks<M: HasBonds + HasSites>(mol: &M) -> Blocks {
    let sites: Vec<SiteId> = mol.sites().collect();
    let n = sites.len();

    if n == 0 {
        return Blocks {
            blocks: vec![],
            site_index: HashMap::new(),
            bond_index: HashMap::new(),
        };
    }

    let site_pos: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

    let mut adj: Vec<Vec<(BondId, usize)>> = vec![vec![]; n];
    for bond in mol.bonds() {
        let (a, b) = mol.bond_endpoints(bond);
        let ai = site_pos[&a];
        let bi = site_pos[&b];
        adj[ai].push((bond, bi));
        adj[bi].push((bond, ai));
    }

    let mut disc = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut timer = 0usize;
    let mut edge_stack: Vec<(BondId, usize, usize)> = Vec::new();
    let mut result_blocks: Vec<Block> = Vec::new();

    let mut dfs_stack: Vec<(usize, Option<BondId>, usize)> = Vec::new();

    for start in 0..n {
        if disc[start] != usize::MAX {
            continue;
        }

        disc[start] = timer;
        low[start] = timer;
        timer += 1;
        dfs_stack.push((start, None, 0));

        while !dfs_stack.is_empty() {
            let (u, parent_bond, adj_pos) = *dfs_stack.last().unwrap();

            if adj_pos < adj[u].len() {
                let (bond, v) = adj[u][adj_pos];
                dfs_stack.last_mut().unwrap().2 += 1;

                if Some(bond) == parent_bond {
                    continue;
                }

                if disc[v] == usize::MAX {
                    edge_stack.push((bond, u, v));
                    disc[v] = timer;
                    low[v] = timer;
                    timer += 1;
                    dfs_stack.push((v, Some(bond), 0));
                } else if disc[v] < disc[u] {
                    edge_stack.push((bond, u, v));
                    if disc[v] < low[u] {
                        low[u] = disc[v];
                    }
                }
            } else {
                dfs_stack.pop();

                if let Some(&(pu, _, _)) = dfs_stack.last() {
                    if low[u] < low[pu] {
                        low[pu] = low[u];
                    }
                    let pb = parent_bond.unwrap();
                    if low[u] >= disc[pu] {
                        let mut block_bonds: Vec<BondId> = Vec::new();
                        let mut block_sites: Vec<usize> = Vec::new();
                        loop {
                            let (b, a, c) = edge_stack.pop().unwrap();
                            block_bonds.push(b);
                            if !block_sites.contains(&a) {
                                block_sites.push(a);
                            }
                            if !block_sites.contains(&c) {
                                block_sites.push(c);
                            }
                            if b == pb {
                                break;
                            }
                        }
                        let is_ring = block_bonds.len() > 1;
                        result_blocks.push(Block {
                            sites: block_sites.into_iter().map(|i| sites[i]).collect(),
                            bonds: block_bonds,
                            is_ring,
                        });
                    }
                }
            }
        }
    }

    let mut site_index: HashMap<SiteId, Vec<usize>> = HashMap::new();
    let mut bond_index: HashMap<BondId, usize> = HashMap::new();
    for (i, block) in result_blocks.iter().enumerate() {
        for &s in block.sites() {
            site_index.entry(s).or_default().push(i);
        }
        for &b in block.bonds() {
            bond_index.insert(b, i);
        }
    }

    Blocks {
        blocks: result_blocks,
        site_index,
        bond_index,
    }
}

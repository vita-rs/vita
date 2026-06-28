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
    /// Sites in this block, in ascending order.
    pub fn sites(&self) -> &[SiteId] {
        &self.sites
    }

    /// Bonds in this block, in ascending order.
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

    /// Iterates all biconnected components, ordered by their sites.
    pub fn iter(&self) -> impl Iterator<Item = &Block> + '_ {
        self.blocks.iter()
    }

    /// Iterates the articulation points (cut sites) of the molecule, in ascending
    /// order.
    ///
    /// A site is an articulation point if removing it would increase the number
    /// of connected components. Equivalently, it appears in two or more blocks.
    pub fn cuts(&self) -> impl Iterator<Item = SiteId> + '_ {
        let mut cuts: Vec<SiteId> = self
            .site_index
            .iter()
            .filter(|(_, v)| v.len() >= 2)
            .map(|(&s, _)| s)
            .collect();
        cuts.sort_unstable();
        cuts.into_iter()
    }

    /// Iterates the bridge bonds of the molecule, in ascending order.
    ///
    /// A bond is a bridge if removing it would increase the number of connected
    /// components. Every bridge forms its own single-bond block
    /// (`is_ring = false`).
    pub fn bridges(&self) -> impl Iterator<Item = BondId> + '_ {
        let mut bridges: Vec<BondId> = self
            .blocks
            .iter()
            .filter(|b| !b.is_ring())
            .flat_map(|b| b.bonds().iter().copied())
            .collect();
        bridges.sort_unstable();
        bridges.into_iter()
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

    /// Iterates all blocks containing `site`, ordered by their sites.
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
/// Isolated sites belong to no block. Blocks are ordered by their sites,
/// ascending within each.
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
                        let mut block_sites: Vec<SiteId> =
                            block_sites.into_iter().map(|i| sites[i]).collect();
                        block_sites.sort_unstable();
                        block_bonds.sort_unstable();
                        result_blocks.push(Block {
                            sites: block_sites,
                            bonds: block_bonds,
                            is_ring,
                        });
                    }
                }
            }
        }
    }

    result_blocks.sort_by(|a, b| a.sites.cmp(&b.sites));

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use vita_core::HasSites;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasBonds for Mol {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.endpoints[i]
        }
    }

    fn empty() -> Mol {
        Mol {
            sites: vec![],
            bonds: vec![],
            endpoints: vec![],
        }
    }

    fn single() -> Mol {
        Mol {
            sites: vec![s(1)],
            bonds: vec![],
            endpoints: vec![],
        }
    }

    fn chain() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
        }
    }

    fn triangle() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3))],
        }
    }

    fn lollipop() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3)), (s(1), s(4))],
        }
    }

    fn dumbbell() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(4), s(6)),
            ],
        }
    }

    #[test]
    fn empty_molecule_has_no_blocks() {
        let blks = blocks(&empty());
        assert_eq!(blks.len(), 0);
        assert!(blks.is_empty());
        assert!(!blks.is_biconnected());
    }

    #[test]
    fn single_site_has_no_blocks() {
        assert_eq!(blocks(&single()).len(), 0);
    }

    #[test]
    fn chain_has_two_blocks() {
        assert_eq!(blocks(&chain()).len(), 2);
    }

    #[test]
    fn triangle_has_one_block() {
        assert_eq!(blocks(&triangle()).len(), 1);
    }

    #[test]
    fn lollipop_has_two_blocks() {
        assert_eq!(blocks(&lollipop()).len(), 2);
    }

    #[test]
    fn triangle_is_biconnected() {
        assert!(blocks(&triangle()).is_biconnected());
    }

    #[test]
    fn chain_is_not_biconnected() {
        assert!(!blocks(&chain()).is_biconnected());
    }

    #[test]
    fn triangle_blocks_are_rings() {
        assert!(blocks(&triangle()).iter().all(|b| b.is_ring()));
    }

    #[test]
    fn chain_blocks_are_not_rings() {
        assert!(blocks(&chain()).iter().all(|b| !b.is_ring()));
    }

    #[test]
    fn triangle_has_no_cuts() {
        assert_eq!(blocks(&triangle()).cuts().count(), 0);
    }

    #[test]
    fn triangle_has_no_bridges() {
        assert_eq!(blocks(&triangle()).bridges().count(), 0);
    }

    #[test]
    fn chain_cuts() {
        let cuts: Vec<SiteId> = blocks(&chain()).cuts().collect();
        assert_eq!(cuts, vec![s(2)]);
    }

    #[test]
    fn lollipop_cuts() {
        let cuts: Vec<SiteId> = blocks(&lollipop()).cuts().collect();
        assert_eq!(cuts, vec![s(1)]);
    }

    #[test]
    fn dumbbell_cuts() {
        let cuts: Vec<SiteId> = blocks(&dumbbell()).cuts().collect();
        assert_eq!(cuts, vec![s(3), s(4)]);
    }

    #[test]
    fn chain_bridges() {
        let bridges: Vec<BondId> = blocks(&chain()).bridges().collect();
        assert_eq!(bridges, vec![b(1), b(2)]);
    }

    #[test]
    fn lollipop_bridges() {
        let bridges: Vec<BondId> = blocks(&lollipop()).bridges().collect();
        assert_eq!(bridges, vec![b(4)]);
    }

    #[test]
    fn dumbbell_bridges() {
        let bridges: Vec<BondId> = blocks(&dumbbell()).bridges().collect();
        assert_eq!(bridges, vec![b(4)]);
    }

    #[test]
    fn blocks_are_independent_of_input_order() {
        let shuffled = Mol {
            sites: vec![s(6), s(4), s(2), s(5), s(1), s(3)],
            bonds: vec![b(7), b(4), b(1), b(6), b(2), b(5), b(3)],
            endpoints: vec![
                (s(4), s(6)),
                (s(3), s(4)),
                (s(1), s(2)),
                (s(5), s(6)),
                (s(2), s(3)),
                (s(4), s(5)),
                (s(1), s(3)),
            ],
        };
        let shape = |m: &Mol| -> (Vec<Vec<SiteId>>, Vec<SiteId>, Vec<BondId>) {
            let bl = blocks(m);
            (
                bl.iter().map(|blk| blk.sites().to_vec()).collect(),
                bl.cuts().collect(),
                bl.bridges().collect(),
            )
        };
        assert_eq!(shape(&dumbbell()), shape(&shuffled));
    }

    #[test]
    fn lollipop_contains_cut() {
        let blks = blocks(&lollipop());
        assert!(blks.contains_cut(s(1)));
        assert!(!blks.contains_cut(s(2)));
        assert!(!blks.contains_cut(s(99)));
    }

    #[test]
    fn lollipop_contains_bridge() {
        let blks = blocks(&lollipop());
        assert!(blks.contains_bridge(b(4)));
        assert!(!blks.contains_bridge(b(1)));
        assert!(!blks.contains_bridge(b(99)));
    }

    #[test]
    fn site_belongs_to_one_block_in_triangle() {
        let blks = blocks(&triangle());
        for site in [s(1), s(2), s(3)] {
            assert_eq!(blks.of_site(site).count(), 1);
        }
    }

    #[test]
    fn articulation_site_belongs_to_multiple_blocks() {
        assert_eq!(blocks(&lollipop()).of_site(s(1)).count(), 2);
    }

    #[test]
    fn isolated_site_has_no_blocks() {
        assert_eq!(blocks(&single()).of_site(s(1)).count(), 0);
    }

    #[test]
    fn unknown_site_has_no_blocks() {
        assert_eq!(blocks(&chain()).of_site(s(99)).count(), 0);
    }

    #[test]
    fn bond_belongs_to_ring_block() {
        let blks = blocks(&triangle());
        assert!(blks.of_bond(b(1)).is_some_and(|blk| blk.is_ring()));
    }

    #[test]
    fn bridge_bond_belongs_to_non_ring_block() {
        let blks = blocks(&lollipop());
        assert!(blks.of_bond(b(4)).is_some_and(|blk| !blk.is_ring()));
    }

    #[test]
    fn unknown_bond_returns_none() {
        assert!(blocks(&chain()).of_bond(b(99)).is_none());
    }

    #[test]
    fn blocks_partition_all_bonds() {
        let mol = dumbbell();
        let blks = blocks(&mol);
        let via_blocks: HashSet<BondId> = blks
            .iter()
            .flat_map(|b| b.bonds().iter().copied())
            .collect();
        let via_mol: HashSet<BondId> = mol.bonds().collect();
        assert_eq!(via_blocks, via_mol);
    }
}

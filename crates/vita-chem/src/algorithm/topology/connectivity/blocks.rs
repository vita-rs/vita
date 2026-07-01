use vita_core::SiteId;

use crate::algorithm::utils::{AdjacencyList, FxHashMap, SortedMap, SortedMultimap};
use crate::{BondId, HasBonds};

/// A biconnected component of a molecule.
///
/// A maximal set of bonds any two of which lie on a common cycle, together with
/// the sites they touch. Each block is either a single bridge bond or a
/// 2-connected subgraph holding at least one cycle.
///
/// Obtain via [`Blocks::iter`], [`Blocks::of_site`], or [`Blocks::of_bond`].
pub struct Block {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
    is_ring: bool,
}

impl Block {
    /// The sites in this block, in ascending order.
    pub fn sites(&self) -> &[SiteId] {
        &self.sites
    }

    /// The bonds in this block, in ascending order.
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
/// Partitions the bonds into maximal 2-connected subgraphs, each a ring system
/// or a single bridge bond. Isolated sites belong to no block. A site shared by
/// several blocks is an articulation point; a bond forming its own block is a
/// bridge.
///
/// Obtain via [`blocks`].
pub struct Blocks {
    blocks: Vec<Block>,
    site_index: SortedMultimap<SiteId, usize>,
    bond_index: SortedMap<BondId, usize>,
    site_count: usize,
}

impl Blocks {
    /// Number of blocks.
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    /// Returns `true` if the molecule has no bonds.
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// Returns `true` if the whole molecule is one biconnected component
    /// spanning every site.
    ///
    /// Equivalently, the molecular graph is connected and has no articulation
    /// point. A molecule with an isolated site, no bonds, or more than one block
    /// is not biconnected.
    pub fn is_biconnected(&self) -> bool {
        self.blocks.len() == 1 && self.blocks[0].sites.len() == self.site_count
    }

    /// Iterates the blocks, ordered by their sites.
    pub fn iter(&self) -> impl Iterator<Item = &Block> + '_ {
        self.blocks.iter()
    }

    /// Iterates the articulation points (cut sites), in ascending order.
    ///
    /// A site is an articulation point when removing it would disconnect the
    /// molecule; equivalently, it lies in two or more blocks.
    pub fn cuts(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.site_index
            .iter()
            .filter(|(_, blocks)| blocks.len() >= 2)
            .map(|(&site, _)| site)
    }

    /// Iterates the bridge bonds, in ascending order.
    ///
    /// A bond is a bridge when removing it would disconnect the molecule;
    /// equivalently, it forms its own single-bond block.
    pub fn bridges(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bond_index
            .iter()
            .filter(|&(_, &i)| !self.blocks[i].is_ring())
            .map(|(&bond, _)| bond)
    }

    /// Returns `true` if `site` is an articulation point.
    ///
    /// Returns `false` if `site` is absent from the molecule or lies in one
    /// block at most.
    pub fn is_cut(&self, site: SiteId) -> bool {
        self.site_index.get(&site).len() >= 2
    }

    /// Returns `true` if `bond` is a bridge.
    ///
    /// Returns `false` if `bond` is absent from the molecule or lies in a ring.
    pub fn is_bridge(&self, bond: BondId) -> bool {
        self.bond_index
            .get(&bond)
            .is_some_and(|&i| !self.blocks[i].is_ring())
    }

    /// Iterates the blocks containing `site`, ordered by their sites.
    ///
    /// Empty if `site` is absent from the molecule or isolated.
    pub fn of_site(&self, site: SiteId) -> impl Iterator<Item = &Block> + '_ {
        self.site_index.get(&site).iter().map(|&i| &self.blocks[i])
    }

    /// Returns the block containing `bond`.
    ///
    /// Returns `None` if `bond` is absent from the molecule.
    pub fn of_bond(&self, bond: BondId) -> Option<&Block> {
        self.bond_index.get(&bond).map(|&i| &self.blocks[i])
    }
}

/// Biconnected components of a molecule.
///
/// Decomposes the bonded structure into maximal 2-connected subgraphs by
/// Tarjan's algorithm, run iteratively so traversal depth is bounded by the heap
/// rather than the call stack. Each block is a ring system or a single bridge
/// bond; isolated sites belong to no block. Blocks are ordered by their sites,
/// ascending within each.
///
/// # Complexity
///
/// O(V + E log E) time and O(V + E) space.
pub fn blocks<M: HasBonds>(mol: &M) -> Blocks {
    let sites: Vec<SiteId> = mol.sites().collect();
    let bonds: Vec<BondId> = mol.bonds().collect();
    let n = sites.len();

    let pos: FxHashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let adjacency = AdjacencyList::build(
        n,
        bonds.iter().enumerate().map(|(edge, &bond)| {
            let (a, b) = mol.bond_endpoints(bond);
            (edge, pos[&a], pos[&b])
        }),
    );

    let mut disc = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut timer = 0usize;
    let mut edge_stack: Vec<(usize, usize, usize)> = Vec::new();
    let mut dfs_stack: Vec<(usize, Option<usize>, usize)> = Vec::new();
    let mut found: Vec<Block> = Vec::new();

    for root in 0..n {
        if disc[root] != usize::MAX {
            continue;
        }
        disc[root] = timer;
        low[root] = timer;
        timer += 1;
        dfs_stack.push((root, None, 0));

        while let Some(&(u, parent_edge, cursor)) = dfs_stack.last() {
            let neighbors = adjacency.neighbors(u);
            if cursor < neighbors.len() {
                let (edge, v) = neighbors[cursor];
                dfs_stack.last_mut().unwrap().2 += 1;
                if Some(edge) == parent_edge {
                    continue;
                }
                if disc[v] == usize::MAX {
                    edge_stack.push((edge, u, v));
                    disc[v] = timer;
                    low[v] = timer;
                    timer += 1;
                    dfs_stack.push((v, Some(edge), 0));
                } else if disc[v] < disc[u] {
                    edge_stack.push((edge, u, v));
                    low[u] = low[u].min(disc[v]);
                }
            } else {
                dfs_stack.pop();
                if let Some(&(parent, _, _)) = dfs_stack.last() {
                    low[parent] = low[parent].min(low[u]);
                    if low[u] >= disc[parent] {
                        let entering = parent_edge.unwrap();
                        found.push(pop_block(&mut edge_stack, entering, &sites, &bonds));
                    }
                }
            }
        }
    }

    found.sort_by(|a, b| a.sites.cmp(&b.sites));

    let bond_index = SortedMap::from_pairs(
        found
            .iter()
            .enumerate()
            .flat_map(|(i, block)| block.bonds.iter().map(move |&bond| (bond, i))),
    );
    let site_index = SortedMultimap::from_pairs(
        found
            .iter()
            .enumerate()
            .flat_map(|(i, block)| block.sites.iter().map(move |&site| (site, i))),
    );

    Blocks {
        blocks: found,
        site_index,
        bond_index,
        site_count: n,
    }
}

/// Pops the edges of one block off the stack, down to and including the
/// `entering` edge, and assembles them into a [`Block`].
fn pop_block(
    edge_stack: &mut Vec<(usize, usize, usize)>,
    entering: usize,
    sites: &[SiteId],
    bonds: &[BondId],
) -> Block {
    let mut block_sites: Vec<SiteId> = Vec::new();
    let mut block_bonds: Vec<BondId> = Vec::new();
    loop {
        let (edge, u, v) = edge_stack.pop().unwrap();
        block_bonds.push(bonds[edge]);
        block_sites.push(sites[u]);
        block_sites.push(sites[v]);
        if edge == entering {
            break;
        }
    }
    let is_ring = block_bonds.len() > 1;
    block_sites.sort_unstable();
    block_sites.dedup();
    block_bonds.sort_unstable();
    Block {
        sites: block_sites,
        bonds: block_bonds,
        is_ring,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    fn mol(sites: &[u32], bonds: &[(u32, u32, u32)]) -> Mol {
        Mol {
            sites: sites.iter().map(|&n| s(n)).collect(),
            bonds: bonds.iter().map(|&(id, _, _)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, u, v)| (s(u), s(v))).collect(),
        }
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
        }
    }

    fn empty() -> Mol {
        mol(&[], &[])
    }

    fn single() -> Mol {
        mol(&[1], &[])
    }

    fn edge() -> Mol {
        mol(&[1, 2], &[(1, 1, 2)])
    }

    fn chain() -> Mol {
        mol(&[1, 2, 3], &[(1, 1, 2), (2, 2, 3)])
    }

    fn triangle() -> Mol {
        mol(&[1, 2, 3], &[(1, 1, 2), (2, 2, 3), (3, 1, 3)])
    }

    fn lollipop() -> Mol {
        mol(&[1, 2, 3, 4], &[(1, 1, 2), (2, 2, 3), (3, 1, 3), (4, 1, 4)])
    }

    fn dumbbell() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 1, 3),
                (4, 3, 4),
                (5, 4, 5),
                (6, 5, 6),
                (7, 4, 6),
            ],
        )
    }

    #[test]
    fn empty_molecule_has_no_blocks() {
        let bl = blocks(&empty());
        assert_eq!(bl.len(), 0);
        assert!(bl.is_empty());
    }

    #[test]
    fn isolated_site_forms_no_block() {
        let bl = blocks(&single());
        assert_eq!(bl.len(), 0);
        assert_eq!(bl.of_site(s(1)).count(), 0);
    }

    #[test]
    fn molecule_with_a_bond_is_not_empty() {
        assert!(!blocks(&edge()).is_empty());
    }

    #[test]
    fn single_bond_forms_one_bridge_block() {
        let bl = blocks(&edge());
        assert_eq!(bl.len(), 1);
        assert!(!bl.iter().next().unwrap().is_ring());
    }

    #[test]
    fn single_bond_block_holds_its_sites_and_bond() {
        let bl = blocks(&edge());
        let block = bl.iter().next().unwrap();
        assert_eq!(block.sites(), &[s(1), s(2)]);
        assert_eq!(block.bonds(), &[b(1)]);
    }

    #[test]
    fn a_cycle_forms_one_ring_block() {
        let bl = blocks(&triangle());
        assert_eq!(bl.len(), 1);
        assert!(bl.iter().next().unwrap().is_ring());
    }

    #[test]
    fn block_lists_its_sites_and_bonds_in_ascending_order() {
        let bl = blocks(&triangle());
        let block = bl.iter().next().unwrap();
        assert_eq!(block.sites(), &[s(1), s(2), s(3)]);
        assert_eq!(block.bonds(), &[b(1), b(2), b(3)]);
    }

    #[test]
    fn unknown_site_is_in_no_block() {
        let bl = blocks(&triangle());
        assert_eq!(bl.of_site(s(99)).count(), 0);
        assert!(!bl.is_cut(s(99)));
    }

    #[test]
    fn unknown_bond_is_in_no_block() {
        let bl = blocks(&triangle());
        assert!(bl.of_bond(b(99)).is_none());
        assert!(!bl.is_bridge(b(99)));
    }

    #[test]
    fn chain_forms_one_block_per_bond() {
        assert_eq!(blocks(&chain()).len(), 2);
    }

    #[test]
    fn chain_middle_site_is_a_cut() {
        assert!(blocks(&chain()).is_cut(s(2)));
    }

    #[test]
    fn chain_endpoints_are_not_cuts() {
        let bl = blocks(&chain());
        assert!(!bl.is_cut(s(1)));
        assert!(!bl.is_cut(s(3)));
    }

    #[test]
    fn a_cycle_has_no_cut_sites() {
        assert_eq!(blocks(&triangle()).cuts().count(), 0);
    }

    #[test]
    fn a_bridge_bond_is_a_bridge() {
        assert!(blocks(&chain()).is_bridge(b(1)));
    }

    #[test]
    fn a_ring_bond_is_not_a_bridge() {
        assert!(!blocks(&triangle()).is_bridge(b(1)));
    }

    #[test]
    fn a_cut_site_belongs_to_every_block_it_joins() {
        assert_eq!(blocks(&lollipop()).of_site(s(1)).count(), 2);
    }

    #[test]
    fn a_non_cut_site_belongs_to_one_block() {
        assert_eq!(blocks(&triangle()).of_site(s(1)).count(), 1);
    }

    #[test]
    fn a_bond_resolves_to_its_block() {
        let bl = blocks(&lollipop());
        assert!(bl.of_bond(b(1)).is_some_and(|block| block.is_ring()));
        assert!(bl.of_bond(b(4)).is_some_and(|block| !block.is_ring()));
    }

    #[test]
    fn cuts_lists_articulation_points_in_ascending_order() {
        assert_eq!(
            blocks(&dumbbell()).cuts().collect::<Vec<_>>(),
            vec![s(3), s(4)]
        );
    }

    #[test]
    fn bridges_lists_bridge_bonds_in_ascending_order() {
        assert_eq!(
            blocks(&dumbbell()).bridges().collect::<Vec<_>>(),
            vec![b(4)]
        );
    }

    #[test]
    fn blocks_are_ordered_by_their_sites() {
        let bl = blocks(&dumbbell());
        let sites: Vec<Vec<SiteId>> = bl.iter().map(|block| block.sites().to_vec()).collect();
        assert_eq!(
            sites,
            vec![
                vec![s(1), s(2), s(3)],
                vec![s(3), s(4)],
                vec![s(4), s(5), s(6)],
            ],
        );
    }

    #[test]
    fn is_biconnected_holds_for_a_single_cycle() {
        assert!(blocks(&triangle()).is_biconnected());
    }

    #[test]
    fn is_biconnected_fails_for_a_chain() {
        assert!(!blocks(&chain()).is_biconnected());
    }

    #[test]
    fn is_biconnected_fails_when_a_site_is_isolated() {
        let m = mol(&[1, 2, 3, 4], &[(1, 1, 2), (2, 2, 3), (3, 1, 3)]);
        assert!(!blocks(&m).is_biconnected());
    }

    #[test]
    fn blocks_partition_every_bond() {
        let m = dumbbell();
        let bl = blocks(&m);
        let mut via_blocks: Vec<BondId> = bl
            .iter()
            .flat_map(|block| block.bonds().iter().copied())
            .collect();
        via_blocks.sort_unstable();
        let mut all: Vec<BondId> = m.bonds().collect();
        all.sort_unstable();
        assert_eq!(via_blocks, all);
    }

    #[test]
    fn output_is_independent_of_input_order() {
        let shape = |m: &Mol| -> (Vec<Vec<SiteId>>, Vec<SiteId>, Vec<BondId>) {
            let bl = blocks(m);
            (
                bl.iter().map(|block| block.sites().to_vec()).collect(),
                bl.cuts().collect(),
                bl.bridges().collect(),
            )
        };
        assert_eq!(shape(&dumbbell()), shape(&reversed(&dumbbell())));
    }
}

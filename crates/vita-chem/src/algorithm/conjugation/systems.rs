use vita_core::{Element, HasElements, SiteId};

use crate::algorithm::utils::{DisjointSet, FxHashMap, SortedMultimap, valence_electrons};
use crate::{BondId, BondOrder, HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// A conjugated system: a maximal set of sites whose p-orbital contributions
/// interact across intervening σ-bonds, with the bonds that carry them.
///
/// Obtain via [`ConjugatedSystems::iter`], [`ConjugatedSystems::of_site`], or
/// [`ConjugatedSystems::of_bond`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConjugatedSystem {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
    electrons: Option<u32>,
}

impl ConjugatedSystem {
    /// Returns `true` if `site` contributes to this system.
    pub fn contains_site(&self, site: SiteId) -> bool {
        self.sites.binary_search(&site).is_ok()
    }

    /// Returns `true` if `bond` lies in this system.
    pub fn contains_bond(&self, bond: BondId) -> bool {
        self.bonds.binary_search(&bond).is_ok()
    }

    /// Number of sites in the system.
    pub fn site_count(&self) -> usize {
        self.sites.len()
    }

    /// Number of bonds in the system.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Iterates the system's sites in ascending order.
    pub fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.sites.iter().copied()
    }

    /// Iterates the system's bonds in ascending order.
    ///
    /// Every bond between two member sites is a member: the multiple bonds
    /// supplying the π-electrons and the σ-bonds the interaction crosses.
    pub fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bonds.iter().copied()
    }

    /// Number of electrons the contributions supply to the system: two per
    /// π-bond and donated lone pair, one per delocalised radical, none from a
    /// vacant orbital.
    ///
    /// Returns `None` when the system holds
    /// [`Aromatic`](BondOrder::Aromatic) bonds, whose delocalised electrons
    /// have no per-system count until the ring is kekulised.
    pub fn pi_electrons(&self) -> Option<u32> {
        self.electrons
    }
}

/// The conjugated systems of a molecule.
///
/// Systems partition the π network, not the molecule: a site or bond lies in
/// up to two systems, one per orthogonal π-plane — an azide's three sites
/// carry two coincident four-electron systems.
///
/// Obtain via [`systems`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConjugatedSystems {
    groups: Vec<ConjugatedSystem>,
    site_index: SortedMultimap<SiteId, usize>,
    bond_index: SortedMultimap<BondId, usize>,
}

impl ConjugatedSystems {
    /// Number of conjugated systems.
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Returns `true` if the molecule conjugates nowhere.
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// Iterates the systems, ordered by their sites.
    pub fn iter(&self) -> impl Iterator<Item = &ConjugatedSystem> + '_ {
        self.groups.iter()
    }

    /// Iterates the systems `site` contributes to, in the order of
    /// [`iter`](Self::iter): none if it conjugates nowhere, two if it lies on
    /// orthogonal π-planes.
    ///
    /// Returns an empty iterator if `site` is absent from the molecule.
    pub fn of_site(&self, site: SiteId) -> impl Iterator<Item = &ConjugatedSystem> + '_ {
        self.site_index
            .get(&site)
            .iter()
            .map(move |&i| &self.groups[i])
    }

    /// Iterates the systems `bond` lies in, in the order of
    /// [`iter`](Self::iter).
    ///
    /// A bond between two cumulated π-planes lies in both. Returns an empty
    /// iterator if `bond` is absent from the molecule or conjugates nowhere.
    pub fn of_bond(&self, bond: BondId) -> impl Iterator<Item = &ConjugatedSystem> + '_ {
        self.bond_index
            .get(&bond)
            .iter()
            .map(move |&i| &self.groups[i])
    }

    /// Returns `true` if some system contains both `a` and `b`.
    ///
    /// Returns `false` if either site is absent from the molecule or no
    /// system holds them together.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        let a_systems = self.site_index.get(&a);
        let b_systems = self.site_index.get(&b);
        a_systems.iter().any(|i| b_systems.contains(i))
    }
}

/// Conjugated systems of a molecule.
///
/// A site contributes to the π network where its octet arithmetic — over the
/// declared bond orders, formal charge, and radical electrons — leaves a
/// p-orbital to offer: a multiple bond, a lone pair, an unpaired electron, or
/// a vacant orbital. Contributions on adjacent sites interact across the
/// intervening σ-bond, always through a multiple bond — lone pairs, radicals,
/// and vacancies extend a π-system, never one another — and every maximal
/// interacting set of two or more contributions is a system. One contribution
/// alone — an isolated double bond, a lone pair beside a saturated
/// neighbour — conjugates nothing,
/// σ-electrons themselves never participate (hyperconjugation is a separate
/// phenomenon), and a site whose arithmetic exceeds its orbitals — hypervalent
/// sulfur or phosphorus, say — neither contributes nor lets conjugation
/// cross. d- and f-block sites carry no fixed valence-electron count and
/// never participate.
///
/// Orthogonality forced by the graph itself is honoured: a cumulated sp
/// centre parts its two π-bonds into perpendicular networks — an even
/// cumulene's ends stay insulated, an odd cumulene's conjugate — and a triple
/// bond runs one network per plane, so a site or bond can lie in two systems
/// at once. Everything else is taken maximally planar: the perception reports
/// what idealised geometry admits, leaving the overlap any one conformer
/// realises to geometric analysis; a cumulated cycle admitting no consistent
/// alternation collapses into a single system.
///
/// [`Aromatic`](BondOrder::Aromatic) bonds join their ring as one
/// declared-delocalised block — aromatic implies conjugated — whose electron
/// count stays indeterminate until kekulised (see
/// [`ConjugatedSystem::pi_electrons`]).
///
/// # Complexity
///
/// O((V + E) · log (V + E)) time and O(V + E) space, over the molecule's `V`
/// sites and `E` bonds; the log factor orders the systems, their sites, and
/// their bonds canonically.
pub fn systems<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>(
    mol: &M,
) -> ConjugatedSystems {
    let mut sites: Vec<SiteId> = mol.sites().collect();
    sites.sort_unstable();
    let n = sites.len();

    let mut rows: Vec<(BondId, usize, usize, BondOrder)> = mol
        .bonds()
        .map(|bond| {
            let (a, b) = mol.bond_endpoints(bond);
            let i = sites.binary_search(&a).unwrap();
            let j = sites.binary_search(&b).unwrap();
            (bond, i.min(j), i.max(j), mol.bond_order(bond))
        })
        .collect();
    rows.sort_unstable_by_key(|&(_, lo, hi, _)| (lo, hi));

    // Per-site bond arithmetic; aromatic declarations bypass it entirely.
    let mut sigma = vec![0u32; n];
    let mut unsaturation = vec![0u32; n];
    let mut aromatic = vec![false; n];
    for &(_, lo, hi, order) in &rows {
        let Some(beyond) = beyond_sigma(order) else {
            aromatic[lo] = true;
            aromatic[hi] = true;
            continue;
        };
        for end in [lo, hi] {
            sigma[end] += 1;
            unsaturation[end] += beyond;
        }
    }

    // Octet arithmetic: who is a Lewis participant, and what p-orbital slots
    // remain beyond the bond system — filled vacancy first, then radical,
    // then lone pair, in reverse order of their pull toward s-character.
    let mut lewis = vec![false; n];
    let mut slots: Vec<[Option<u32>; 2]> = vec![[None; 2]; n];
    for i in 0..n {
        if aromatic[i] {
            continue;
        }
        let element = mol.element(sites[i]);
        let Some(electrons) = valence_electrons(element) else {
            continue;
        };
        let charge = i32::from(mol.formal_charge(sites[i]));
        let radicals = u32::from(mol.radical_electron(sites[i]));
        let bonding = sigma[i] + unsaturation[i];
        let free = i32::from(electrons) - charge - bonding as i32 - radicals as i32;
        if free < 0 {
            continue;
        }
        let pairs = free as u32 / 2;
        let used = bonding + pairs + radicals;
        let budget = orbitals(element);
        if used > budget {
            continue;
        }
        lewis[i] = true;

        let perpendicular: u32 = match sigma[i] {
            0 => 0,
            1 => 2,
            2 if unsaturation[i] >= 2 => 2,
            2 | 3 => 1,
            _ => 0,
        };
        let extra = perpendicular.saturating_sub(unsaturation[i]) as usize;
        let mut filled = 0;
        for (count, contribution) in [(budget - used, 0), (radicals, 1), (pairs, 2)] {
            for _ in 0..count {
                if filled == extra {
                    break;
                }
                slots[i][filled] = Some(contribution);
                filled += 1;
            }
        }
    }

    // π units — one node per π-bond — and which units each site hosts.
    let mut units: Vec<(usize, usize)> = Vec::new();
    let mut row_units = vec![false; rows.len()];
    let mut hosted: Vec<[Option<usize>; 2]> = vec![[None; 2]; n];
    for (row, &(_, lo, hi, order)) in rows.iter().enumerate() {
        let pi = beyond_sigma(order).map_or(0, |beyond| beyond.min(2));
        if pi == 0 || !lewis[lo] || !lewis[hi] {
            continue;
        }
        row_units[row] = true;
        for _ in 0..pi {
            let unit = units.len();
            units.push((lo, hi));
            for end in [lo, hi] {
                let host = &mut hosted[end];
                if host[0].is_none() {
                    host[0] = Some(unit);
                } else {
                    host[1] = Some(unit);
                }
            }
        }
    }

    // Rigid clusters: units sharing a site sit on perpendicular planes, so
    // 2-colour each cluster; a cluster its cycle leaves uncolourable is
    // geometrically frustrated and collapses into one channel.
    let mut adjacent: Vec<Vec<usize>> = vec![Vec::new(); units.len()];
    for host in &hosted {
        if let [Some(a), Some(b)] = *host {
            adjacent[a].push(b);
            adjacent[b].push(a);
        }
    }
    let mut cluster = vec![usize::MAX; units.len()];
    let mut color = vec![0u8; units.len()];
    let mut frustrated = Vec::new();
    for root in 0..units.len() {
        if cluster[root] != usize::MAX {
            continue;
        }
        let c = frustrated.len();
        frustrated.push(false);
        cluster[root] = c;
        let mut stack = vec![root];
        while let Some(u) = stack.pop() {
            for &v in &adjacent[u] {
                if cluster[v] == usize::MAX {
                    cluster[v] = c;
                    color[v] = 1 - color[u];
                    stack.push(v);
                } else if color[v] == color[u] {
                    frustrated[c] = true;
                }
            }
        }
    }

    // Order each pair of hosted units principal plane first, so junctions and
    // slots meet the two planes of both sides in matching order.
    for host in &mut hosted {
        if let [Some(a), Some(b)] = *host
            && color[a] > color[b]
        {
            host.swap(0, 1);
        }
    }

    // Declared-aromatic blocks: connected components over aromatic bonds.
    let mut block_dsu = DisjointSet::new(n);
    for &(_, lo, hi, order) in &rows {
        if order == BondOrder::Aromatic {
            block_dsu.union(lo, hi);
        }
    }
    let mut block_of: Vec<Option<usize>> = vec![None; n];
    let mut block_sites: Vec<Vec<usize>> = Vec::new();
    let mut block_index: FxHashMap<usize, usize> = FxHashMap::default();
    for i in 0..n {
        if !aromatic[i] {
            continue;
        }
        let root = block_dsu.find(i);
        let block = *block_index.entry(root).or_insert_with(|| {
            block_sites.push(Vec::new());
            block_sites.len() - 1
        });
        block_of[i] = Some(block);
        block_sites[block].push(i);
    }

    // Channel merging: nodes are the π units plus one node per aromatic
    // block; the channels present at a site, principal plane first.
    let mut channels = DisjointSet::new(units.len() + block_sites.len());
    let presence = |i: usize| -> [Option<usize>; 2] {
        match block_of[i] {
            Some(block) => [Some(units.len() + block), None],
            None => hosted[i],
        }
    };

    // Units of one cluster and colour share a channel; frustration merges
    // both colours.
    let mut leaders: FxHashMap<(usize, u8), usize> = FxHashMap::default();
    for unit in 0..units.len() {
        let key = (
            cluster[unit],
            if frustrated[cluster[unit]] {
                0
            } else {
                color[unit]
            },
        );
        match leaders.get(&key) {
            Some(&leader) => {
                channels.union(leader, unit);
            }
            None => {
                leaders.insert(key, unit);
            }
        }
    }

    // Junctions: a σ-bond between two channel-bearing sites aligns their
    // planes pairwise under the maximally planar geometry.
    for (row, &(_, lo, hi, order)) in rows.iter().enumerate() {
        if order == BondOrder::Aromatic || row_units[row] {
            continue;
        }
        let (at_lo, at_hi) = (presence(lo), presence(hi));
        if let (Some(a), Some(b)) = (at_lo[0], at_hi[0]) {
            channels.union(a, b);
            if let (Some(a), Some(b)) = (at_lo[1], at_hi[1]) {
                channels.union(a, b);
            }
        }
    }

    // Slot attachment: a pure donor offers its slots to the neighbouring
    // planes in order; the leftover slot of a π-bearing terminus reaches only
    // the plane perpendicular to its own π-bond.
    let mut slot_targets: Vec<[Vec<usize>; 2]> = (0..n).map(|_| Default::default()).collect();
    for &(_, lo, hi, order) in &rows {
        if order == BondOrder::Aromatic {
            continue;
        }
        for (from, to) in [(lo, hi), (hi, lo)] {
            if slots[from][0].is_none() {
                continue;
            }
            let reachable = presence(to);
            if hosted[from][0].is_none() {
                for (k, target) in reachable.into_iter().enumerate() {
                    if slots[from][k].is_some()
                        && let Some(node) = target
                    {
                        slot_targets[from][k].push(node);
                    }
                }
            } else {
                for node in reachable.into_iter().flatten() {
                    if hosted[from].iter().flatten().all(|&own| own != node) {
                        slot_targets[from][0].push(node);
                    }
                }
            }
        }
    }
    let mut attached: Vec<(usize, u32, usize)> = Vec::new();
    for i in 0..n {
        for k in 0..2 {
            let targets = &slot_targets[i][k];
            let (Some(&anchor), Some(contribution)) = (targets.first(), slots[i][k]) else {
                continue;
            };
            for pair in targets.windows(2) {
                channels.union(pair[0], pair[1]);
            }
            attached.push((i, contribution, anchor));
        }
    }

    // Gather channels, keep those with at least two contributions, and
    // materialise canonically ordered systems.
    let mut bucket_index: FxHashMap<usize, usize> = FxHashMap::default();
    let mut buckets: Vec<(u32, Vec<usize>, Option<u32>)> = Vec::new();
    let mut bucket = |root: usize, buckets: &mut Vec<(u32, Vec<usize>, Option<u32>)>| -> usize {
        *bucket_index.entry(root).or_insert_with(|| {
            buckets.push((0, Vec::new(), Some(0)));
            buckets.len() - 1
        })
    };
    for (unit, &(lo, hi)) in units.iter().enumerate() {
        let at = bucket(channels.find(unit), &mut buckets);
        let (count, members, electrons) = &mut buckets[at];
        *count += 1;
        members.extend([lo, hi]);
        *electrons = electrons.map(|total| total + 2);
    }
    for (block, members_of_block) in block_sites.iter().enumerate() {
        let at = bucket(channels.find(units.len() + block), &mut buckets);
        let (count, members, electrons) = &mut buckets[at];
        *count += 2;
        members.extend(members_of_block);
        *electrons = None;
    }
    for &(site, contribution, anchor) in &attached {
        let at = bucket(channels.find(anchor), &mut buckets);
        let (count, members, electrons) = &mut buckets[at];
        *count += 1;
        members.push(site);
        *electrons = electrons.map(|total| total + contribution);
    }

    let mut kept: Vec<(Vec<usize>, Option<u32>)> = buckets
        .into_iter()
        .filter(|&(count, ..)| count >= 2)
        .map(|(_, mut members, electrons)| {
            members.sort_unstable();
            members.dedup();
            (members, electrons)
        })
        .collect();
    kept.sort_unstable();

    let mut memberships: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (g, (members, _)) in kept.iter().enumerate() {
        for &i in members {
            memberships[i].push(g);
        }
    }
    let mut group_bonds: Vec<Vec<BondId>> = vec![Vec::new(); kept.len()];
    for &(bond, lo, hi, _) in &rows {
        for &g in &memberships[lo] {
            if memberships[hi].contains(&g) {
                group_bonds[g].push(bond);
            }
        }
    }

    let mut groups: Vec<ConjugatedSystem> = kept
        .into_iter()
        .zip(group_bonds)
        .map(|((members, electrons), mut bonds)| {
            bonds.sort_unstable();
            ConjugatedSystem {
                sites: members.into_iter().map(|i| sites[i]).collect(),
                bonds,
                electrons,
            }
        })
        .collect();
    groups.sort_by(|a, b| {
        a.sites
            .cmp(&b.sites)
            .then_with(|| a.bonds.cmp(&b.bonds))
            .then_with(|| a.electrons.cmp(&b.electrons))
    });

    let site_index = SortedMultimap::from_pairs(
        groups
            .iter()
            .enumerate()
            .flat_map(|(g, system)| system.sites.iter().map(move |&site| (site, g))),
    );
    let bond_index = SortedMultimap::from_pairs(
        groups
            .iter()
            .enumerate()
            .flat_map(|(g, system)| system.bonds.iter().map(move |&bond| (bond, g))),
    );

    ConjugatedSystems {
        groups,
        site_index,
        bond_index,
    }
}

/// Bonding orbitals an order engages at each end beyond its σ — π, δ, and φ
/// components alike; `None` for the delocalised aromatic order.
fn beyond_sigma(order: BondOrder) -> Option<u32> {
    match order {
        BondOrder::Single => Some(0),
        BondOrder::Double => Some(1),
        BondOrder::Triple => Some(2),
        BondOrder::Quadruple => Some(3),
        BondOrder::Quintuple => Some(4),
        BondOrder::Hextuple => Some(5),
        BondOrder::Aromatic => None,
    }
}

/// Valence orbitals a main-group element bonds and holds electrons with: one
/// in the first shell, four beyond.
fn orbitals(element: Element) -> u32 {
    if element.atomic_number() <= 2 { 1 } else { 4 }
}

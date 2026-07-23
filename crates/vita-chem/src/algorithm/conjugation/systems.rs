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
    donations: Vec<(SiteId, u32)>,
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
    /// π-bond and donated lone pair, one per delocalized radical, none from a
    /// vacant orbital.
    ///
    /// Returns `None` when the system holds
    /// [`Aromatic`](BondOrder::Aromatic) bonds, whose delocalized electrons
    /// have no per-system count until the ring is kekulized.
    pub fn pi_electrons(&self) -> Option<u32> {
        self.electrons
    }

    /// Number of lone pairs `site` donates into this system's π network.
    ///
    /// Zero when `site` donates nothing: it contributes through its π-bonds
    /// alone, lies outside the system, or is absent from the molecule. The
    /// attribution follows the declared Lewis form — an amide's neutral form
    /// donates from the nitrogen, its zwitterionic form from the oxygen —
    /// while the partition and its electron count stay form-independent.
    pub fn donated_pairs(&self, site: SiteId) -> u32 {
        self.donations
            .binary_search_by_key(&site, |&(donor, _)| donor)
            .map_or(0, |i| self.donations[i].1)
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
/// neighbor — conjugates nothing,
/// σ-electrons themselves never participate (hyperconjugation is a separate
/// phenomenon), and a site whose arithmetic exceeds its orbitals — hypervalent
/// sulfur or phosphorus, say — neither contributes nor lets conjugation
/// cross. d- and f-block sites carry no fixed valence-electron count and
/// never participate.
///
/// Orthogonality forced by the graph itself is honored: a cumulated sp
/// center parts its two π-bonds into perpendicular networks — an even
/// cumulene's ends stay insulated, an odd cumulene's conjugate — and a triple
/// bond runs one network per plane, so a site or bond can lie in two systems
/// at once. Everything else is taken maximally planar: the perception reports
/// what idealised geometry admits, leaving the overlap any one conformer
/// realizes to geometric analysis; a cumulated cycle admitting no consistent
/// alternation collapses into a single system.
///
/// [`Aromatic`](BondOrder::Aromatic) bonds join their ring as one
/// declared-delocalized block — aromatic implies conjugated — whose electron
/// count stays indeterminate until kekulized (see
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
    // 2-color each cluster; a cluster its cycle leaves uncolorable is
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

    // Units of one cluster and color share a channel; frustration merges
    // both colors.
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

    // Slot attachment: a pure donor offers its slots to the neighboring
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
    // materialize canonically ordered systems.
    let mut bucket_index: FxHashMap<usize, usize> = FxHashMap::default();
    let mut buckets: Vec<Bucket> = Vec::new();
    let mut bucket = |root: usize, buckets: &mut Vec<Bucket>| -> usize {
        *bucket_index.entry(root).or_insert_with(|| {
            buckets.push((0, Vec::new(), Some(0), Vec::new()));
            buckets.len() - 1
        })
    };
    for (unit, &(lo, hi)) in units.iter().enumerate() {
        let at = bucket(channels.find(unit), &mut buckets);
        let (count, members, electrons, _) = &mut buckets[at];
        *count += 1;
        members.extend([lo, hi]);
        *electrons = electrons.map(|total| total + 2);
    }
    for (block, members_of_block) in block_sites.iter().enumerate() {
        let at = bucket(channels.find(units.len() + block), &mut buckets);
        let (count, members, electrons, _) = &mut buckets[at];
        *count += 2;
        members.extend(members_of_block);
        *electrons = None;
    }
    for &(site, contribution, anchor) in &attached {
        let at = bucket(channels.find(anchor), &mut buckets);
        let (count, members, electrons, donors) = &mut buckets[at];
        *count += 1;
        members.push(site);
        *electrons = electrons.map(|total| total + contribution);
        if contribution == 2 {
            donors.push(site);
        }
    }

    let mut kept: Vec<Channel> = buckets
        .into_iter()
        .filter(|&(count, ..)| count >= 2)
        .map(|(_, mut members, electrons, mut donors)| {
            members.sort_unstable();
            members.dedup();
            donors.sort_unstable();
            let mut donations: Vec<(usize, u32)> = Vec::new();
            for donor in donors {
                match donations.last_mut() {
                    Some((site, pairs)) if *site == donor => *pairs += 1,
                    _ => donations.push((donor, 1)),
                }
            }
            (members, electrons, donations)
        })
        .collect();
    kept.sort_unstable();

    let mut memberships: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (g, (members, ..)) in kept.iter().enumerate() {
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
        .map(|((members, electrons, donations), mut bonds)| {
            bonds.sort_unstable();
            ConjugatedSystem {
                sites: members.into_iter().map(|i| sites[i]).collect(),
                bonds,
                electrons,
                donations: donations
                    .into_iter()
                    .map(|(i, pairs)| (sites[i], pairs))
                    .collect(),
            }
        })
        .collect();
    groups.sort_by(|a, b| {
        a.sites
            .cmp(&b.sites)
            .then_with(|| a.bonds.cmp(&b.bonds))
            .then_with(|| a.electrons.cmp(&b.electrons))
            .then_with(|| a.donations.cmp(&b.donations))
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

/// One channel being gathered: its contribution count, member site indices,
/// running electron count, and lone-pair donor site indices.
type Bucket = (u32, Vec<usize>, Option<u32>, Vec<usize>);

/// A channel past the two-contribution threshold: its member site indices,
/// electron count, and donated pairs per donor site index.
type Channel = (Vec<usize>, Option<u32>, Vec<(usize, u32)>);

/// Bonding orbitals an order engages at each end beyond its σ — π, δ, and φ
/// components alike; `None` for the delocalized aromatic order.
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

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::HasBonds;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
        formal_charges: Vec<i8>,
        radicals: Vec<u8>,
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasElements for Mol {
        fn element(&self, site: SiteId) -> Element {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.elements[i]
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

    impl HasBondOrders for Mol {
        fn bond_order(&self, bond: BondId) -> BondOrder {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.orders[i]
        }
    }

    impl HasFormalCharges for Mol {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.formal_charges[i]
        }
    }

    impl HasRadicalElectrons for Mol {
        fn radical_electron(&self, site: SiteId) -> u8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.radicals[i]
        }
    }

    fn molecule(atoms: &[(u32, &str, i8, u8)], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, ..)| s(id)).collect(),
            elements: atoms.iter().map(|&(_, symbol, ..)| elem(symbol)).collect(),
            formal_charges: atoms.iter().map(|&(_, _, charge, _)| charge).collect(),
            radicals: atoms.iter().map(|&(.., radicals)| radicals).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, u, v, _)| (s(u), s(v))).collect(),
            orders: bonds.iter().map(|&(.., order)| order).collect(),
        }
    }

    fn shape(perceived: &ConjugatedSystems) -> Vec<(Vec<SiteId>, Vec<BondId>, Option<u32>)> {
        perceived
            .iter()
            .map(|system| {
                (
                    system.sites().collect(),
                    system.bonds().collect(),
                    system.pi_electrons(),
                )
            })
            .collect()
    }

    fn empty() -> Mol {
        molecule(&[], &[])
    }

    fn ethane() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 2, 6, BondOrder::Single),
                (6, 2, 7, BondOrder::Single),
                (7, 2, 8, BondOrder::Single),
            ],
        )
    }

    fn ethylene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 2, 5, BondOrder::Single),
                (5, 2, 6, BondOrder::Single),
            ],
        )
    }

    fn formaldehyde() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "O", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
            ],
        )
    }

    fn penta_1_4_diene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 3, 4, BondOrder::Single),
                (4, 4, 5, BondOrder::Double),
                (5, 1, 6, BondOrder::Single),
                (6, 1, 7, BondOrder::Single),
                (7, 2, 8, BondOrder::Single),
                (8, 3, 9, BondOrder::Single),
                (9, 3, 10, BondOrder::Single),
                (10, 4, 11, BondOrder::Single),
                (11, 5, 12, BondOrder::Single),
                (12, 5, 13, BondOrder::Single),
            ],
        )
    }

    fn butadiene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 3, 4, BondOrder::Double),
                (4, 1, 5, BondOrder::Single),
                (5, 1, 6, BondOrder::Single),
                (6, 2, 7, BondOrder::Single),
                (7, 3, 8, BondOrder::Single),
                (8, 4, 9, BondOrder::Single),
                (9, 4, 10, BondOrder::Single),
            ],
        )
    }

    fn vinyl_chloride() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "Cl", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 2, 6, BondOrder::Single),
            ],
        )
    }

    fn hydrazine() -> Mol {
        molecule(
            &[
                (1, "N", 0, 0),
                (2, "N", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 2, 5, BondOrder::Single),
                (5, 2, 6, BondOrder::Single),
            ],
        )
    }

    fn allyl(charge: i8, radicals: u8) -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", charge, radicals),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 2, 6, BondOrder::Single),
                (6, 3, 7, BondOrder::Single),
                (7, 3, 8, BondOrder::Single),
            ],
        )
    }

    fn acetate() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "O", 0, 0),
                (4, "O", -1, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Double),
                (3, 2, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 1, 6, BondOrder::Single),
                (6, 1, 7, BondOrder::Single),
            ],
        )
    }

    fn divinyl_ketone() -> Mol {
        molecule(
            &[
                (1, "O", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 3, 4, BondOrder::Double),
                (4, 2, 5, BondOrder::Single),
                (5, 5, 6, BondOrder::Double),
                (6, 3, 7, BondOrder::Single),
                (7, 4, 8, BondOrder::Single),
                (8, 4, 9, BondOrder::Single),
                (9, 5, 10, BondOrder::Single),
                (10, 6, 11, BondOrder::Single),
                (11, 6, 12, BondOrder::Single),
            ],
        )
    }

    fn kekule_benzene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 3, 4, BondOrder::Double),
                (4, 4, 5, BondOrder::Single),
                (5, 5, 6, BondOrder::Double),
                (6, 6, 1, BondOrder::Single),
                (7, 1, 7, BondOrder::Single),
                (8, 2, 8, BondOrder::Single),
                (9, 3, 9, BondOrder::Single),
                (10, 4, 10, BondOrder::Single),
                (11, 5, 11, BondOrder::Single),
                (12, 6, 12, BondOrder::Single),
            ],
        )
    }

    fn kekule_furan() -> Mol {
        molecule(
            &[
                (1, "O", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Double),
                (3, 3, 4, BondOrder::Single),
                (4, 4, 5, BondOrder::Double),
                (5, 5, 1, BondOrder::Single),
                (6, 2, 6, BondOrder::Single),
                (7, 3, 7, BondOrder::Single),
                (8, 4, 8, BondOrder::Single),
                (9, 5, 9, BondOrder::Single),
            ],
        )
    }

    fn divinyl_ether() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "O", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 3, 4, BondOrder::Single),
                (4, 4, 5, BondOrder::Double),
                (5, 1, 6, BondOrder::Single),
                (6, 1, 7, BondOrder::Single),
                (7, 2, 8, BondOrder::Single),
                (8, 4, 9, BondOrder::Single),
                (9, 5, 10, BondOrder::Single),
                (10, 5, 11, BondOrder::Single),
            ],
        )
    }

    fn allene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Double),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
                (6, 3, 7, BondOrder::Single),
            ],
        )
    }

    fn butatriene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Double),
                (3, 3, 4, BondOrder::Double),
                (4, 1, 5, BondOrder::Single),
                (5, 1, 6, BondOrder::Single),
                (6, 4, 7, BondOrder::Single),
                (7, 4, 8, BondOrder::Single),
            ],
        )
    }

    fn pentatetraene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Double),
                (3, 3, 4, BondOrder::Double),
                (4, 4, 5, BondOrder::Double),
                (5, 1, 6, BondOrder::Single),
                (6, 1, 7, BondOrder::Single),
                (7, 5, 8, BondOrder::Single),
                (8, 5, 9, BondOrder::Single),
            ],
        )
    }

    fn vinylacetylene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 3, 4, BondOrder::Triple),
                (4, 1, 5, BondOrder::Single),
                (5, 1, 6, BondOrder::Single),
                (6, 2, 7, BondOrder::Single),
                (7, 4, 8, BondOrder::Single),
            ],
        )
    }

    fn diacetylene() -> Mol {
        molecule(
            &[
                (1, "H", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 2, 3, BondOrder::Triple),
                (2, 3, 4, BondOrder::Single),
                (3, 4, 5, BondOrder::Triple),
                (4, 1, 2, BondOrder::Single),
                (5, 5, 6, BondOrder::Single),
            ],
        )
    }

    fn azide() -> Mol {
        molecule(
            &[(1, "N", -1, 0), (2, "N", 1, 0), (3, "N", -1, 0)],
            &[(1, 1, 2, BondOrder::Double), (2, 2, 3, BondOrder::Double)],
        )
    }

    fn chloroacetylene() -> Mol {
        molecule(
            &[
                (1, "Cl", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Triple),
                (3, 3, 4, BondOrder::Single),
            ],
        )
    }

    fn acrylonitrile() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "N", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 3, 4, BondOrder::Triple),
                (4, 1, 5, BondOrder::Single),
                (5, 1, 6, BondOrder::Single),
                (6, 2, 7, BondOrder::Single),
            ],
        )
    }

    fn cumulated_cycle() -> Mol {
        molecule(
            &[(1, "C", 0, 0), (2, "C", 0, 0), (3, "C", 0, 0)],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Double),
                (3, 3, 1, BondOrder::Double),
            ],
        )
    }

    fn dimethyl_sulfone() -> Mol {
        molecule(
            &[
                (1, "S", 0, 0),
                (2, "O", 0, 0),
                (3, "O", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "H", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 1, 3, BondOrder::Double),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 4, 6, BondOrder::Single),
                (6, 4, 7, BondOrder::Single),
                (7, 4, 8, BondOrder::Single),
                (8, 5, 9, BondOrder::Single),
                (9, 5, 10, BondOrder::Single),
                (10, 5, 11, BondOrder::Single),
            ],
        )
    }

    fn iron_vinyl() -> Mol {
        molecule(
            &[
                (1, "Fe", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Double),
                (3, 2, 4, BondOrder::Single),
                (4, 3, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
            ],
        )
    }

    fn aromatic_benzene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Aromatic),
                (2, 2, 3, BondOrder::Aromatic),
                (3, 3, 4, BondOrder::Aromatic),
                (4, 4, 5, BondOrder::Aromatic),
                (5, 5, 6, BondOrder::Aromatic),
                (6, 6, 1, BondOrder::Aromatic),
                (7, 1, 7, BondOrder::Single),
                (8, 2, 8, BondOrder::Single),
                (9, 3, 9, BondOrder::Single),
                (10, 4, 10, BondOrder::Single),
                (11, 5, 11, BondOrder::Single),
                (12, 6, 12, BondOrder::Single),
            ],
        )
    }

    fn aromatic_styrene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "C", 0, 0),
                (8, "C", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
                (15, "H", 0, 0),
                (16, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Aromatic),
                (2, 2, 3, BondOrder::Aromatic),
                (3, 3, 4, BondOrder::Aromatic),
                (4, 4, 5, BondOrder::Aromatic),
                (5, 5, 6, BondOrder::Aromatic),
                (6, 6, 1, BondOrder::Aromatic),
                (7, 1, 7, BondOrder::Single),
                (8, 7, 8, BondOrder::Double),
                (9, 2, 9, BondOrder::Single),
                (10, 3, 10, BondOrder::Single),
                (11, 4, 11, BondOrder::Single),
                (12, 5, 12, BondOrder::Single),
                (13, 6, 13, BondOrder::Single),
                (14, 7, 14, BondOrder::Single),
                (15, 8, 15, BondOrder::Single),
                (16, 8, 16, BondOrder::Single),
            ],
        )
    }

    fn aniline() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "N", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
                (13, "H", 0, 0),
                (14, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Aromatic),
                (2, 2, 3, BondOrder::Aromatic),
                (3, 3, 4, BondOrder::Aromatic),
                (4, 4, 5, BondOrder::Aromatic),
                (5, 5, 6, BondOrder::Aromatic),
                (6, 6, 1, BondOrder::Aromatic),
                (7, 1, 7, BondOrder::Single),
                (8, 7, 8, BondOrder::Single),
                (9, 7, 9, BondOrder::Single),
                (10, 2, 10, BondOrder::Single),
                (11, 3, 11, BondOrder::Single),
                (12, 4, 12, BondOrder::Single),
                (13, 5, 13, BondOrder::Single),
                (14, 6, 14, BondOrder::Single),
            ],
        )
    }

    fn formamide() -> Mol {
        molecule(
            &[
                (1, "O", 0, 0),
                (2, "C", 0, 0),
                (3, "N", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 2, 4, BondOrder::Single),
                (4, 3, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
            ],
        )
    }

    fn formamide_polar() -> Mol {
        molecule(
            &[
                (1, "O", -1, 0),
                (2, "C", 0, 0),
                (3, "N", 1, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Double),
                (3, 2, 4, BondOrder::Single),
                (4, 3, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
            ],
        )
    }

    fn two_fragments() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "Cl", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
                (11, "C", 0, 0),
                (12, "C", 0, 0),
                (13, "C", 0, 0),
                (14, "C", 0, 0),
                (15, "H", 0, 0),
                (16, "H", 0, 0),
                (17, "H", 0, 0),
                (18, "H", 0, 0),
                (19, "H", 0, 0),
                (20, "H", 0, 0),
            ],
            &[
                (1, 11, 12, BondOrder::Double),
                (2, 12, 13, BondOrder::Single),
                (3, 13, 14, BondOrder::Double),
                (4, 1, 2, BondOrder::Double),
                (5, 2, 3, BondOrder::Single),
                (6, 1, 4, BondOrder::Single),
                (7, 1, 5, BondOrder::Single),
                (8, 2, 6, BondOrder::Single),
                (9, 11, 15, BondOrder::Single),
                (10, 11, 16, BondOrder::Single),
                (11, 12, 17, BondOrder::Single),
                (12, 13, 18, BondOrder::Single),
                (13, 14, 19, BondOrder::Single),
                (14, 14, 20, BondOrder::Single),
            ],
        )
    }

    #[test]
    fn an_empty_molecule_has_no_systems() {
        let perceived = systems(&empty());
        assert_eq!(perceived.len(), 0);
        assert!(perceived.is_empty());
        assert_eq!(perceived.iter().count(), 0);
    }

    #[test]
    fn a_saturated_molecule_has_no_systems() {
        assert!(systems(&ethane()).is_empty());
    }

    #[test]
    fn an_isolated_multiple_bond_is_not_a_system() {
        assert!(systems(&ethylene()).is_empty());
    }

    #[test]
    fn an_isolated_carbonyl_is_not_a_system() {
        assert!(systems(&formaldehyde()).is_empty());
    }

    #[test]
    fn a_saturated_site_does_not_bridge() {
        assert!(systems(&penta_1_4_diene()).is_empty());
    }

    #[test]
    fn alternating_multiple_bonds_conjugate() {
        let perceived = systems(&butadiene());
        assert_eq!(
            shape(&perceived),
            vec![(
                vec![s(1), s(2), s(3), s(4)],
                vec![b(1), b(2), b(3)],
                Some(4)
            )]
        );
        let system = perceived.iter().next().unwrap();
        assert_eq!(system.site_count(), 4);
        assert_eq!(system.bond_count(), 3);
        assert!(system.contains_site(s(1)));
        assert!(!system.contains_site(s(5)));
        assert!(system.contains_bond(b(2)));
        assert!(!system.contains_bond(b(4)));
    }

    #[test]
    fn a_lone_pair_joins_an_adjacent_pi_bond() {
        assert_eq!(
            shape(&systems(&vinyl_chloride())),
            vec![(vec![s(1), s(2), s(3)], vec![b(1), b(2)], Some(4))]
        );
    }

    #[test]
    fn adjacent_lone_pairs_alone_do_not_conjugate() {
        assert!(systems(&hydrazine()).is_empty());
    }

    #[test]
    fn a_radical_joins_an_adjacent_pi_bond() {
        assert_eq!(
            shape(&systems(&allyl(0, 1))),
            vec![(vec![s(1), s(2), s(3)], vec![b(1), b(2)], Some(3))]
        );
    }

    #[test]
    fn a_vacant_orbital_joins_an_adjacent_pi_bond() {
        assert_eq!(
            shape(&systems(&allyl(1, 0))),
            vec![(vec![s(1), s(2), s(3)], vec![b(1), b(2)], Some(2))]
        );
    }

    #[test]
    fn a_charged_donor_joins_an_adjacent_pi_bond() {
        let perceived = systems(&acetate());
        assert_eq!(
            shape(&perceived),
            vec![(vec![s(2), s(3), s(4)], vec![b(2), b(3)], Some(4))]
        );
        assert_eq!(perceived.iter().next().unwrap().donated_pairs(s(4)), 1);
    }

    #[test]
    fn cross_conjugation_merges_into_one_system() {
        assert_eq!(
            shape(&systems(&divinyl_ketone())),
            vec![(
                vec![s(1), s(2), s(3), s(4), s(5), s(6)],
                vec![b(1), b(2), b(3), b(4), b(5)],
                Some(6)
            )]
        );
    }

    #[test]
    fn a_kekule_ring_is_one_system() {
        assert_eq!(
            shape(&systems(&kekule_benzene())),
            vec![(
                vec![s(1), s(2), s(3), s(4), s(5), s(6)],
                vec![b(1), b(2), b(3), b(4), b(5), b(6)],
                Some(6)
            )]
        );
    }

    #[test]
    fn a_kekule_heteroring_counts_its_donated_pair() {
        let perceived = systems(&kekule_furan());
        assert_eq!(
            shape(&perceived),
            vec![(
                vec![s(1), s(2), s(3), s(4), s(5)],
                vec![b(1), b(2), b(3), b(4), b(5)],
                Some(6)
            )]
        );
        assert_eq!(perceived.iter().next().unwrap().donated_pairs(s(1)), 1);
    }

    #[test]
    fn a_donor_bridges_disjoint_pi_bonds() {
        assert_eq!(
            shape(&systems(&divinyl_ether())),
            vec![(
                vec![s(1), s(2), s(3), s(4), s(5)],
                vec![b(1), b(2), b(3), b(4)],
                Some(6)
            )]
        );
    }

    #[test]
    fn an_even_cumulene_does_not_conjugate_through() {
        assert!(systems(&allene()).is_empty());
    }

    #[test]
    fn an_odd_cumulene_conjugates_end_to_end() {
        assert_eq!(
            shape(&systems(&butatriene())),
            vec![(
                vec![s(1), s(2), s(3), s(4)],
                vec![b(1), b(2), b(3)],
                Some(4)
            )]
        );
    }

    #[test]
    fn a_longer_even_cumulene_conjugates_within_each_plane() {
        let perceived = systems(&pentatetraene());
        assert_eq!(
            shape(&perceived),
            vec![
                (
                    vec![s(1), s(2), s(3), s(4)],
                    vec![b(1), b(2), b(3)],
                    Some(4)
                ),
                (
                    vec![s(2), s(3), s(4), s(5)],
                    vec![b(2), b(3), b(4)],
                    Some(4)
                ),
            ]
        );
        assert_eq!(perceived.of_bond(b(2)).count(), 2);
    }

    #[test]
    fn a_triple_bond_conjugates_in_one_plane() {
        assert_eq!(
            shape(&systems(&vinylacetylene())),
            vec![(
                vec![s(1), s(2), s(3), s(4)],
                vec![b(1), b(2), b(3)],
                Some(4)
            )]
        );
    }

    #[test]
    fn adjacent_triple_bonds_conjugate_in_both_planes() {
        assert_eq!(
            shape(&systems(&diacetylene())),
            vec![
                (
                    vec![s(2), s(3), s(4), s(5)],
                    vec![b(1), b(2), b(3)],
                    Some(4)
                ),
                (
                    vec![s(2), s(3), s(4), s(5)],
                    vec![b(1), b(2), b(3)],
                    Some(4)
                ),
            ]
        );
    }

    #[test]
    fn orthogonal_systems_coincide_on_their_sites() {
        let perceived = systems(&azide());
        assert_eq!(
            shape(&perceived),
            vec![
                (vec![s(1), s(2), s(3)], vec![b(1), b(2)], Some(4)),
                (vec![s(1), s(2), s(3)], vec![b(1), b(2)], Some(4)),
            ]
        );
        assert_eq!(perceived.of_site(s(2)).count(), 2);
        assert_eq!(perceived.of_bond(b(1)).count(), 2);
    }

    #[test]
    fn a_terminal_donor_reaches_both_orthogonal_planes() {
        assert_eq!(
            shape(&systems(&chloroacetylene())),
            vec![
                (vec![s(1), s(2), s(3)], vec![b(1), b(2)], Some(4)),
                (vec![s(1), s(2), s(3)], vec![b(1), b(2)], Some(4)),
            ]
        );
    }

    #[test]
    fn a_nitrile_lone_pair_points_along_the_axis() {
        let perceived = systems(&acrylonitrile());
        assert_eq!(
            shape(&perceived),
            vec![(
                vec![s(1), s(2), s(3), s(4)],
                vec![b(1), b(2), b(3)],
                Some(4)
            )]
        );
        assert_eq!(perceived.iter().next().unwrap().donated_pairs(s(4)), 0);
    }

    #[test]
    fn a_frustrated_cumulated_cycle_merges_into_one_system() {
        assert_eq!(
            shape(&systems(&cumulated_cycle())),
            vec![(vec![s(1), s(2), s(3)], vec![b(1), b(2), b(3)], Some(6))]
        );
    }

    #[test]
    fn conjugation_does_not_cross_a_hypervalent_center() {
        assert!(systems(&dimethyl_sulfone()).is_empty());
    }

    #[test]
    fn a_d_block_site_never_participates() {
        assert!(systems(&iron_vinyl()).is_empty());
    }

    #[test]
    fn a_declared_aromatic_ring_is_a_system_with_indeterminate_electrons() {
        assert_eq!(
            shape(&systems(&aromatic_benzene())),
            vec![(
                vec![s(1), s(2), s(3), s(4), s(5), s(6)],
                vec![b(1), b(2), b(3), b(4), b(5), b(6)],
                None
            )]
        );
    }

    #[test]
    fn a_declared_aromatic_ring_conjugates_with_its_substituent() {
        assert_eq!(
            shape(&systems(&aromatic_styrene())),
            vec![(
                vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8)],
                vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8)],
                None
            )]
        );
    }

    #[test]
    fn a_donor_joins_a_declared_aromatic_ring() {
        assert_eq!(
            shape(&systems(&aniline())),
            vec![(
                vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7)],
                vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7)],
                None
            )]
        );
    }

    #[test]
    fn systems_iterate_in_canonical_order() {
        let perceived = systems(&two_fragments());
        assert_eq!(
            shape(&perceived),
            vec![
                (vec![s(1), s(2), s(3)], vec![b(4), b(5)], Some(4)),
                (
                    vec![s(11), s(12), s(13), s(14)],
                    vec![b(1), b(2), b(3)],
                    Some(4)
                ),
            ]
        );
    }

    #[test]
    fn same_relates_sites_sharing_a_system() {
        let perceived = systems(&butadiene());
        assert!(perceived.same(s(1), s(4)));
        assert!(!perceived.same(s(1), s(5)));
        assert!(!perceived.same(s(1), s(99)));
    }

    #[test]
    fn same_is_false_across_orthogonal_planes() {
        let perceived = systems(&pentatetraene());
        assert!(perceived.same(s(1), s(4)));
        assert!(!perceived.same(s(1), s(5)));
    }

    #[test]
    fn of_site_of_an_absent_site_is_empty() {
        assert_eq!(systems(&butadiene()).of_site(s(99)).count(), 0);
    }

    #[test]
    fn of_site_of_a_saturated_site_is_empty() {
        assert_eq!(systems(&butadiene()).of_site(s(5)).count(), 0);
    }

    #[test]
    fn of_bond_of_an_unconjugated_bond_is_empty() {
        assert_eq!(systems(&butadiene()).of_bond(b(4)).count(), 0);
    }

    #[test]
    fn donated_pairs_counts_lone_pairs_given_to_the_system() {
        let perceived = systems(&vinyl_chloride());
        let system = perceived.iter().next().unwrap();
        assert_eq!(system.donated_pairs(s(3)), 1);
        assert_eq!(system.donated_pairs(s(1)), 0);
        assert_eq!(system.donated_pairs(s(4)), 0);
        assert_eq!(system.donated_pairs(s(99)), 0);
    }

    #[test]
    fn a_radical_or_vacant_donor_donates_no_pairs() {
        let radical = systems(&allyl(0, 1));
        let cation = systems(&allyl(1, 0));
        assert_eq!(radical.iter().next().unwrap().donated_pairs(s(3)), 0);
        assert_eq!(cation.iter().next().unwrap().donated_pairs(s(3)), 0);
    }

    #[test]
    fn a_terminal_donor_donates_one_pair_per_plane() {
        let perceived = systems(&chloroacetylene());
        for system in perceived.iter() {
            assert_eq!(system.donated_pairs(s(1)), 1);
        }
    }

    #[test]
    fn orthogonal_planes_split_their_donations() {
        let perceived = systems(&azide());
        let split: Vec<(u32, u32)> = perceived
            .iter()
            .map(|system| (system.donated_pairs(s(1)), system.donated_pairs(s(3))))
            .collect();
        assert_eq!(split, vec![(1, 0), (0, 1)]);
    }

    #[test]
    fn donation_stays_exact_in_a_declared_aromatic_system() {
        let perceived = systems(&aniline());
        let system = perceived.iter().next().unwrap();
        assert_eq!(system.donated_pairs(s(7)), 1);
        assert_eq!(system.pi_electrons(), None);
    }

    #[test]
    fn the_systems_are_independent_of_the_resonance_form() {
        let neutral = systems(&formamide());
        let polar = systems(&formamide_polar());
        assert_eq!(shape(&neutral), shape(&polar));
        assert_eq!(
            shape(&neutral),
            vec![(vec![s(1), s(2), s(3)], vec![b(1), b(2)], Some(4))]
        );
    }

    #[test]
    fn donation_attribution_follows_the_declared_form() {
        let neutral = systems(&formamide());
        let polar = systems(&formamide_polar());
        let from_neutral = neutral.iter().next().unwrap();
        let from_polar = polar.iter().next().unwrap();
        assert_eq!(from_neutral.donated_pairs(s(3)), 1);
        assert_eq!(from_neutral.donated_pairs(s(1)), 0);
        assert_eq!(from_polar.donated_pairs(s(1)), 1);
        assert_eq!(from_polar.donated_pairs(s(3)), 0);
    }

    #[test]
    fn the_systems_are_independent_of_input_order() {
        let shuffled = molecule(
            &[
                (8, "H", 0, 0),
                (3, "C", 0, 0),
                (10, "H", 0, 0),
                (1, "C", 0, 0),
                (6, "H", 0, 0),
                (4, "C", 0, 0),
                (9, "H", 0, 0),
                (2, "C", 0, 0),
                (5, "H", 0, 0),
                (7, "H", 0, 0),
            ],
            &[
                (7, 3, 8, BondOrder::Single),
                (3, 3, 4, BondOrder::Double),
                (9, 4, 10, BondOrder::Single),
                (1, 1, 2, BondOrder::Double),
                (5, 1, 6, BondOrder::Single),
                (8, 4, 9, BondOrder::Single),
                (2, 2, 3, BondOrder::Single),
                (6, 2, 7, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
            ],
        );
        assert_eq!(systems(&butadiene()), systems(&shuffled));
    }
}

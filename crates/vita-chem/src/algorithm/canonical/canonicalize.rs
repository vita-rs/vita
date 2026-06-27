use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use vita_core::SiteId;

use crate::utils::refine;
use crate::{BondId, HasBonds};

/// The canonical labelling and form of a molecule's sites.
///
/// [`canonicalize`] assigns every site a rank — its place in a total order fixed
/// by the molecular graph and the colouring it was built with, not by the order
/// the sites were given in — and records the labelled graph in that order as a
/// canonical form. Two molecules the colouring makes isomorphic share that form,
/// so a `Canonical` is a portable identity: compare it to test sameness, hash it
/// to key a registry, sort molecules into a stable order. Equivalent atoms of a
/// symmetric molecule are interchangeable and may take their shared ranks in
/// either arrangement, but the form — and with it the identity — is the same.
///
/// Obtain via [`canonicalize`].
#[derive(Debug)]
pub struct Canonical<VK, EK> {
    order: Vec<SiteId>,
    ranks: HashMap<SiteId, usize>,
    // The canonical form: each site's key in rank order, then every bond as its
    // ordered endpoint ranks and key. Equal exactly for molecules the keys make
    // isomorphic — the basis of `Eq`, `Ord`, and `Hash`.
    form: (Vec<VK>, Vec<(usize, usize, EK)>),
}

impl<VK, EK> Canonical<VK, EK> {
    /// Number of ranked sites.
    pub fn len(&self) -> usize {
        self.order.len()
    }

    /// Returns `true` if the molecule contains no sites.
    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// Returns the canonical rank of `site`.
    ///
    /// Returns `None` if `site` is not present in the molecule.
    pub fn rank(&self, site: SiteId) -> Option<usize> {
        self.ranks.get(&site).copied()
    }

    /// Iterates the sites in canonical order, rank `0` first.
    pub fn order(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.order.iter().copied()
    }
}

impl<VK: PartialEq, EK: PartialEq> PartialEq for Canonical<VK, EK> {
    fn eq(&self, other: &Self) -> bool {
        self.form == other.form
    }
}

impl<VK: Eq, EK: Eq> Eq for Canonical<VK, EK> {}

impl<VK: Ord, EK: Ord> PartialOrd for Canonical<VK, EK> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<VK: Ord, EK: Ord> Ord for Canonical<VK, EK> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.form.cmp(&other.form)
    }
}

impl<VK: Hash, EK: Hash> Hash for Canonical<VK, EK> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.form.hash(state);
    }
}

/// Canonically ranks the sites of a molecule.
///
/// The colouring is the caller's: `site_key` and `bond_key` map each site and
/// bond to an ordered key, and two atoms or bonds count as the same exactly when
/// their keys are equal. That choice is the definition of identity — pass the
/// element and bond order to rank by constitution, fold in the formal charge to
/// tell charge states apart — and the library imposes no default. The sites are
/// ranked, and the labelled graph recorded, into the one canonical form their
/// structure and keys dictate: molecules the keys make isomorphic yield equal
/// [`Canonical`]s, whatever sequence the sites and bonds arrived in.
///
/// The ranking is exact. Colour refinement settles the sites into classes by
/// their coloured neighbourhoods; where symmetry leaves a class unsplit the
/// search individualises each of its members in turn and keeps the labelling of
/// least certificate. Taking the least over every branch — rather than
/// committing to a greedy tie-break — is what frees the result from the input
/// order.
///
/// # Complexity
///
/// O(V · (V + E) · log V) per refinement, run once for an asymmetric molecule
/// but once per node of the search where symmetry forces a choice: near-linear in
/// the node count in practice, yet exponential in the worst case, as canonical
/// labelling is not known to be polynomial.
pub fn canonicalize<M, VK, EK>(
    mol: &M,
    site_key: impl Fn(SiteId) -> VK,
    bond_key: impl Fn(BondId) -> EK,
) -> Canonical<VK, EK>
where
    M: HasBonds,
    VK: Ord,
    EK: Ord,
{
    let sites: Vec<SiteId> = mol.sites().collect();
    let n = sites.len();
    if n == 0 {
        return Canonical {
            order: Vec::new(),
            ranks: HashMap::new(),
            form: (Vec::new(), Vec::new()),
        };
    }

    let pos: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let site_keys: Vec<VK> = sites.iter().map(|&s| site_key(s)).collect();
    let seed = dense(&site_keys);

    let bonds: Vec<BondId> = mol.bonds().collect();
    let bond_keys: Vec<EK> = bonds.iter().map(|&b| bond_key(b)).collect();
    let edge = dense(&bond_keys);
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (i, &bond) in bonds.iter().enumerate() {
        let (a, b) = mol.bond_endpoints(bond);
        adjacency[pos[&a]].push((pos[&b], edge[i]));
        adjacency[pos[&b]].push((pos[&a], edge[i]));
    }

    let mut search = Search {
        adjacency: &adjacency,
        seed: &seed,
        n,
        best: None,
    };
    search.run(seed.clone());
    let labelling = search.best.expect("the search visits at least one leaf").1;

    let mut order = vec![sites[0]; n];
    let mut ranks = HashMap::with_capacity(n);
    for (i, &site) in sites.iter().enumerate() {
        order[labelling[i]] = site;
        ranks.insert(site, labelling[i]);
    }

    let mut keyed: Vec<(usize, VK)> = site_keys
        .into_iter()
        .enumerate()
        .map(|(i, key)| (labelling[i], key))
        .collect();
    keyed.sort_by_key(|&(rank, _)| rank);
    let form_sites: Vec<VK> = keyed.into_iter().map(|(_, key)| key).collect();

    let mut form_bonds: Vec<(usize, usize, EK)> = bonds
        .iter()
        .zip(bond_keys)
        .map(|(&bond, key)| {
            let (a, b) = mol.bond_endpoints(bond);
            let (ra, rb) = (labelling[pos[&a]], labelling[pos[&b]]);
            (ra.min(rb), ra.max(rb), key)
        })
        .collect();
    form_bonds.sort();

    Canonical {
        order,
        ranks,
        form: (form_sites, form_bonds),
    }
}

/// Dense ranks of values: equal values share a rank, distinct values rank by
/// `Ord`, so the result depends only on the values and never on their position.
fn dense<K: Ord>(keys: &[K]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..keys.len()).collect();
    order.sort_by(|&a, &b| keys[a].cmp(&keys[b]));
    let mut ranks = vec![0; keys.len()];
    let mut rank = 0;
    for i in 1..order.len() {
        if keys[order[i]] != keys[order[i - 1]] {
            rank += 1;
        }
        ranks[order[i]] = rank;
    }
    ranks
}

/// The individualisation–refinement search for the least certificate.
struct Search<'a> {
    adjacency: &'a [Vec<(usize, usize)>],
    seed: &'a [usize],
    n: usize,
    best: Option<(Vec<usize>, Vec<usize>)>,
}

impl Search<'_> {
    /// Refines, then records a discrete colouring as a leaf or individualises a
    /// target cell vertex by vertex.
    fn run(&mut self, mut colours: Vec<usize>) {
        let count = refine(self.adjacency, &mut colours);
        if count == self.n {
            let certificate = self.certificate(&colours);
            let better = match &self.best {
                Some((best, _)) => certificate < *best,
                None => true,
            };
            if better {
                self.best = Some((certificate, colours));
            }
            return;
        }
        for v in self.target(&colours, count) {
            let mut next: Vec<usize> = colours.iter().map(|&c| c * 2).collect();
            next[v] += 1;
            self.run(next);
        }
    }

    /// The vertices of the smallest non-singleton colour class. Ties fall to the
    /// lowest colour, itself canonical, so the branch taken — and the least
    /// certificate it leads to — does not depend on input order.
    fn target(&self, colours: &[usize], count: usize) -> Vec<usize> {
        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); count];
        for (v, &c) in colours.iter().enumerate() {
            cells[c].push(v);
        }
        cells
            .into_iter()
            .filter(|cell| cell.len() > 1)
            .min_by_key(Vec::len)
            .expect("a non-discrete colouring has a non-singleton class")
    }

    /// The certificate of a discrete colouring: the seeded, labelled graph
    /// written out in rank order. The least certificate over all leaves picks the
    /// canonical labelling.
    fn certificate(&self, colours: &[usize]) -> Vec<usize> {
        let mut by_rank = vec![0; self.n];
        for (v, &rank) in colours.iter().enumerate() {
            by_rank[rank] = v;
        }
        let mut certificate = Vec::new();
        for &v in &by_rank {
            certificate.push(self.seed[v]);
            let mut incident: Vec<(usize, usize)> = self.adjacency[v]
                .iter()
                .map(|&(u, e)| (e, colours[u]))
                .collect();
            incident.sort_unstable();
            certificate.push(incident.len());
            for (edge, neighbour) in incident {
                certificate.push(edge);
                certificate.push(neighbour);
            }
        }
        certificate
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondOrder, HasBondOrders};
    use std::collections::HashSet;
    use vita_core::{Element, HasElements, HasSites};

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

    fn canon(mol: &Mol) -> Canonical<Element, u8> {
        canonicalize(
            mol,
            |site| mol.element(site),
            |bond| mol.bond_order(bond) as u8,
        )
    }

    fn empty() -> Mol {
        Mol {
            sites: vec![],
            elements: vec![],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        }
    }

    fn methane() -> Mol {
        Mol {
            sites: vec![s(1)],
            elements: vec![elem("C")],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        }
    }

    fn cyanate() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("N"), elem("C"), elem("O")],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    fn cyanate_shuffled() -> Mol {
        Mol {
            sites: vec![s(3), s(1), s(2)],
            elements: vec![elem("O"), elem("N"), elem("C")],
            bonds: vec![b(2), b(1)],
            endpoints: vec![(s(2), s(3)), (s(1), s(2))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    fn cyanate_relabelled() -> Mol {
        Mol {
            sites: vec![s(4), s(5), s(6)],
            elements: vec![elem("N"), elem("C"), elem("O")],
            bonds: vec![b(3), b(4)],
            endpoints: vec![(s(4), s(5)), (s(5), s(6))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    fn ring6() -> Mol {
        Mol {
            sites: (1..=6).map(s).collect(),
            elements: vec![elem("C"); 6],
            bonds: (1..=6).map(b).collect(),
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
            ],
            orders: vec![BondOrder::Single; 6],
        }
    }

    fn ring6_rotated() -> Mol {
        Mol {
            sites: vec![s(4), s(5), s(6), s(1), s(2), s(3)],
            elements: vec![elem("C"); 6],
            bonds: vec![b(4), b(5), b(6), b(1), b(2), b(3)],
            endpoints: vec![
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
            ],
            orders: vec![BondOrder::Single; 6],
        }
    }

    fn butane() -> Mol {
        Mol {
            sites: (1..=4).map(s).collect(),
            elements: vec![elem("C"); 4],
            bonds: (1..=3).map(b).collect(),
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4))],
            orders: vec![BondOrder::Single; 3],
        }
    }

    fn isobutane() -> Mol {
        Mol {
            sites: (1..=4).map(s).collect(),
            elements: vec![elem("C"); 4],
            bonds: (1..=3).map(b).collect(),
            endpoints: vec![(s(1), s(2)), (s(1), s(3)), (s(1), s(4))],
            orders: vec![BondOrder::Single; 3],
        }
    }

    fn fragments() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("C"), elem("C"), elem("O")],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
            orders: vec![BondOrder::Single],
        }
    }

    #[test]
    fn empty_molecule_is_empty() {
        let canonical = canon(&empty());
        assert_eq!(canonical.len(), 0);
        assert!(canonical.is_empty());
    }

    #[test]
    fn single_site_ranks_zero() {
        let canonical = canon(&methane());
        assert_eq!(canonical.len(), 1);
        assert_eq!(canonical.rank(s(1)), Some(0));
    }

    #[test]
    fn ranks_are_a_permutation() {
        let canonical = canon(&cyanate());
        let ranks: HashSet<usize> = cyanate()
            .sites
            .iter()
            .map(|&x| canonical.rank(x).unwrap())
            .collect();
        assert_eq!(ranks, (0..3).collect());
    }

    #[test]
    fn order_inverts_rank() {
        let canonical = canon(&cyanate());
        for (r, site) in canonical.order().enumerate() {
            assert_eq!(canonical.rank(site), Some(r));
        }
    }

    #[test]
    fn unknown_site_has_no_rank() {
        assert_eq!(canon(&cyanate()).rank(s(99)), None);
    }

    #[test]
    fn asymmetric_labelling_is_independent_of_input_order() {
        let plain = canon(&cyanate());
        let shuffled = canon(&cyanate_shuffled());
        for site in [s(1), s(2), s(3)] {
            assert_eq!(plain.rank(site), shuffled.rank(site));
        }
    }

    #[test]
    fn least_atomic_number_ranks_first() {
        let canonical = canon(&cyanate());
        assert_eq!(canonical.rank(s(2)), Some(0));
        assert_eq!(canonical.rank(s(1)), Some(1));
        assert_eq!(canonical.rank(s(3)), Some(2));
    }

    #[test]
    fn symmetric_ring_ranks_every_atom() {
        let canonical = canon(&ring6());
        let ranks: HashSet<usize> = (1..=6).map(|i| canonical.rank(s(i)).unwrap()).collect();
        assert_eq!(ranks, (0..6).collect());
    }

    #[test]
    fn disconnected_molecule_ranks_every_atom() {
        let canonical = canon(&fragments());
        let ranks: HashSet<usize> = (1..=3).map(|i| canonical.rank(s(i)).unwrap()).collect();
        assert_eq!(ranks, (0..3).collect());
    }

    #[test]
    fn reordered_molecule_is_equal() {
        assert_eq!(canon(&cyanate()), canon(&cyanate_shuffled()));
    }

    #[test]
    fn relabelled_molecule_is_equal() {
        assert_eq!(canon(&cyanate()), canon(&cyanate_relabelled()));
    }

    #[test]
    fn symmetric_molecule_identity_is_independent_of_input_order() {
        assert_eq!(canon(&ring6()), canon(&ring6_rotated()));
    }

    #[test]
    fn different_elements_differ() {
        let mut thiocyanate = cyanate();
        thiocyanate.elements = vec![elem("O"), elem("C"), elem("S")];
        assert_ne!(canon(&cyanate()), canon(&thiocyanate));
    }

    #[test]
    fn different_bond_orders_differ() {
        let mut doubled = cyanate();
        doubled.orders = vec![BondOrder::Double, BondOrder::Single];
        assert_ne!(canon(&cyanate()), canon(&doubled));
    }

    #[test]
    fn topology_is_part_of_identity() {
        assert_ne!(canon(&butane()), canon(&isobutane()));
    }

    #[test]
    fn equal_canonicals_hash_and_order_alike() {
        let set: HashSet<Canonical<Element, u8>> =
            [canon(&cyanate()), canon(&cyanate_relabelled())]
                .into_iter()
                .collect();
        assert_eq!(set.len(), 1);
        assert_eq!(
            canon(&cyanate()).cmp(&canon(&cyanate_shuffled())),
            Ordering::Equal
        );
    }
}

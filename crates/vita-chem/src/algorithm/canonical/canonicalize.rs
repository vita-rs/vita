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

use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

use vita_core::SiteId;

use crate::algorithm::utils::{FxHashMap, SortedMap, labeling};
use crate::{BondId, HasBonds};

/// The canonical labeling of a molecule's sites: a portable identity, with its
/// symmetry classes.
///
/// [`canonicalize`] assigns every site a rank — its place in a total order fixed
/// by the molecular graph and the coloring it was built with, not by the order
/// the sites were given in — and records the labeled graph in that order as a
/// canonical form. Two molecules the coloring makes isomorphic share that form,
/// so a `Canonical` is a portable identity: compare it to test sameness, hash it
/// to key a registry, sort molecules into a stable order.
///
/// The search that fixes the order also uncovers the molecule's symmetry. Sites
/// an automorphism can interchange form an orbit — the equivalent atoms of the
/// structure — and may take their shared ranks in either arrangement; the form,
/// and with it the identity, is the same however they fall.
///
/// Obtain via [`canonicalize`].
pub struct Canonical<VK, EK> {
    order: Vec<SiteId>,
    ranks: SortedMap<SiteId, usize>,
    classes: Vec<Vec<SiteId>>,
    class_of: SortedMap<SiteId, usize>,
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
    /// Returns `None` if `site` is absent from the molecule.
    pub fn rank(&self, site: SiteId) -> Option<usize> {
        self.ranks.get(&site).copied()
    }

    /// Iterates the sites in canonical order, rank `0` first.
    pub fn order(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.order.iter().copied()
    }

    /// Iterates the symmetry classes, each a set of sites an automorphism can
    /// interchange. Classes are ordered by their sites, ascending within each.
    pub fn orbits(&self) -> impl Iterator<Item = &[SiteId]> + '_ {
        self.classes.iter().map(Vec::as_slice)
    }

    /// Returns the symmetry class containing `site`.
    ///
    /// Returns `None` if `site` is absent from the molecule.
    pub fn orbit(&self, site: SiteId) -> Option<&[SiteId]> {
        self.class_of
            .get(&site)
            .map(|&class| self.classes[class].as_slice())
    }

    /// Returns `true` if `a` and `b` are interchangeable by a symmetry of the
    /// molecule.
    ///
    /// Returns `false` if either site is absent from the molecule.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.class_of.get(&a), self.class_of.get(&b)) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        }
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
/// The coloring is the caller's: `site_key` and `bond_key` map each site and
/// bond to an ordered key, and two atoms or bonds count as the same exactly when
/// their keys are equal. That choice is the definition of identity — pass the
/// element and bond order to rank by constitution, fold in the formal charge to
/// tell charge states apart — and the library imposes no default. The sites are
/// ranked, and the labeled graph recorded, into the one canonical form their
/// structure and keys dictate: molecules the keys make isomorphic yield equal
/// [`Canonical`]s, whatever sequence the sites and bonds arrived in.
///
/// The ranking is exact. Color refinement settles the sites into classes by
/// their colored neighborhoods; where symmetry leaves a class unsplit the search
/// individualizes each of its members in turn and keeps the labeling of least
/// certificate. Taking the least over every branch — rather than committing to a
/// greedy tie-break — is what frees the result from the input order, and the
/// automorphisms it meets along the way are the molecule's
/// [`orbits`](Canonical::orbits).
///
/// # Complexity
///
/// O(V · (V + E) · log V) time per refinement and O(V + E) space, over the
/// molecule's `V` sites and `E` bonds — one refinement for a rigid molecule, one
/// per search node where symmetry forces a branch; near-linear in practice,
/// exponential in the worst case.
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
    let pos: FxHashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

    let site_keys: Vec<VK> = sites.iter().map(|&site| site_key(site)).collect();
    let seed = dense(&site_keys);

    let bonds: Vec<BondId> = mol.bonds().collect();
    let bond_keys: Vec<EK> = bonds.iter().map(|&bond| bond_key(bond)).collect();
    let colors = dense(&bond_keys);

    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (i, &bond) in bonds.iter().enumerate() {
        let (a, b) = mol.bond_endpoints(bond);
        adjacency[pos[&a]].push((pos[&b], colors[i]));
        adjacency[pos[&b]].push((pos[&a], colors[i]));
    }

    let labeled = labeling(&adjacency, &seed);
    let rank = labeled.ranks();
    let orbit = labeled.orbits();

    let mut order = sites.clone();
    for (vertex, &site) in sites.iter().enumerate() {
        order[rank[vertex]] = site;
    }

    let ranks = SortedMap::from_pairs(
        sites
            .iter()
            .enumerate()
            .map(|(vertex, &site)| (site, rank[vertex])),
    );

    let mut classes: Vec<Vec<SiteId>> = vec![Vec::new(); n];
    for (vertex, &site) in sites.iter().enumerate() {
        classes[orbit[vertex]].push(site);
    }
    classes.retain(|class| !class.is_empty());
    for class in &mut classes {
        class.sort_unstable();
    }
    classes.sort_unstable();

    let class_of = SortedMap::from_pairs(
        classes
            .iter()
            .enumerate()
            .flat_map(|(i, class)| class.iter().map(move |&site| (site, i))),
    );

    let mut keyed: Vec<(usize, VK)> = site_keys
        .into_iter()
        .enumerate()
        .map(|(vertex, key)| (rank[vertex], key))
        .collect();
    keyed.sort_unstable_by_key(|entry| entry.0);
    let form_sites: Vec<VK> = keyed.into_iter().map(|(_, key)| key).collect();

    let mut form_bonds: Vec<(usize, usize, EK)> = bonds
        .iter()
        .zip(bond_keys)
        .map(|(&bond, key)| {
            let (a, b) = mol.bond_endpoints(bond);
            let (a, b) = (rank[pos[&a]], rank[pos[&b]]);
            (a.min(b), a.max(b), key)
        })
        .collect();
    form_bonds.sort();

    Canonical {
        order,
        ranks,
        classes,
        class_of,
        form: (form_sites, form_bonds),
    }
}

/// Dense ranks of values: equal values share a rank, distinct values rank by
/// `Ord`, so the result depends only on the values, never on their position.
fn dense<K: Ord>(keys: &[K]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..keys.len()).collect();
    order.sort_unstable_by(|&a, &b| keys[a].cmp(&keys[b]));
    let mut ranks = vec![0; keys.len()];
    let mut rank = 0;
    for window in order.windows(2) {
        if keys[window[1]] != keys[window[0]] {
            rank += 1;
        }
        ranks[window[1]] = rank;
    }
    ranks
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

    fn canon(mol: &Mol) -> Canonical<u32, u32> {
        canonicalize(mol, |_| 0u32, |_| 0u32)
    }

    fn hash_of(canonical: &Canonical<u32, u32>) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        canonical.hash(&mut hasher);
        hasher.finish()
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

    fn path() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
        }
    }

    fn shifted_path() -> Mol {
        Mol {
            sites: vec![s(5), s(6), s(7)],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(5), s(6)), (s(6), s(7))],
        }
    }

    fn reordered_path() -> Mol {
        Mol {
            sites: vec![s(3), s(1), s(2)],
            bonds: vec![b(2), b(1)],
            endpoints: vec![(s(2), s(3)), (s(1), s(2))],
        }
    }

    fn triangle() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3))],
        }
    }

    fn star() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(1), s(3)), (s(1), s(4))],
        }
    }

    #[test]
    fn empty_molecule_is_empty() {
        let c = canon(&empty());
        assert_eq!(c.len(), 0);
        assert!(c.is_empty());
    }

    #[test]
    fn single_site_has_length_one() {
        let c = canon(&single());
        assert_eq!(c.len(), 1);
        assert!(!c.is_empty());
    }

    #[test]
    fn single_site_ranks_zero() {
        assert_eq!(canon(&single()).rank(s(1)), Some(0));
    }

    #[test]
    fn single_site_forms_one_orbit() {
        let c = canon(&single());
        let orbits: Vec<Vec<SiteId>> = c.orbits().map(<[SiteId]>::to_vec).collect();
        assert_eq!(orbits, vec![vec![s(1)]]);
    }

    #[test]
    fn order_enumerates_every_site() {
        let c = canon(&star());
        assert_eq!(c.len(), 4);
        let mut sites: Vec<SiteId> = c.order().collect();
        sites.sort_unstable();
        assert_eq!(sites, vec![s(1), s(2), s(3), s(4)]);
    }

    #[test]
    fn rank_is_the_position_in_the_canonical_order() {
        let c = canon(&star());
        for (position, site) in c.order().enumerate() {
            assert_eq!(c.rank(site), Some(position));
        }
    }

    #[test]
    fn symmetric_sites_share_an_orbit() {
        assert!(canon(&path()).same(s(1), s(3)));
    }

    #[test]
    fn orbits_group_the_symmetric_sites() {
        let c = canon(&star());
        let orbits: Vec<Vec<SiteId>> = c.orbits().map(<[SiteId]>::to_vec).collect();
        assert_eq!(orbits, vec![vec![s(1)], vec![s(2), s(3), s(4)]]);
    }

    #[test]
    fn orbit_of_a_site_is_its_symmetry_class() {
        let c = canon(&star());
        assert_eq!(c.orbit(s(2)), Some([s(2), s(3), s(4)].as_slice()));
        assert_eq!(c.orbit(s(1)), Some([s(1)].as_slice()));
    }

    #[test]
    fn rank_of_an_absent_site_is_none() {
        assert_eq!(canon(&path()).rank(s(99)), None);
    }

    #[test]
    fn same_is_false_for_an_absent_site() {
        assert!(!canon(&path()).same(s(1), s(99)));
    }

    #[test]
    fn orbit_of_an_absent_site_is_none() {
        assert_eq!(canon(&path()).orbit(s(99)), None);
    }

    #[test]
    fn asymmetric_sites_lie_in_different_orbits() {
        assert!(!canon(&path()).same(s(1), s(2)));
    }

    #[test]
    fn a_symmetric_ring_is_a_single_orbit() {
        let c = canon(&triangle());
        assert_eq!(c.orbits().count(), 1);
        assert!(c.same(s(1), s(2)));
        assert!(c.same(s(2), s(3)));
    }

    #[test]
    fn a_site_key_splits_a_symmetric_pair() {
        let c = canonicalize(&path(), |site| (site == s(1)) as u32, |_| 0u32);
        assert!(!c.same(s(1), s(3)));
    }

    #[test]
    fn a_bond_key_splits_a_symmetric_pair() {
        let c = canonicalize(&path(), |_| 0u32, |bond| (bond == b(1)) as u32);
        assert!(!c.same(s(1), s(3)));
    }

    #[test]
    fn isomorphic_molecules_share_a_canonical_form() {
        assert!(canon(&path()) == canon(&shifted_path()));
    }

    #[test]
    fn distinct_structures_are_not_equal() {
        assert!(canon(&path()) != canon(&triangle()));
    }

    #[test]
    fn a_site_key_distinguishes_identity() {
        let plain = canon(&path());
        let colored = canonicalize(&path(), |site| (site == s(1)) as u32, |_| 0u32);
        assert!(plain != colored);
    }

    #[test]
    fn a_bond_key_distinguishes_identity() {
        let plain = canon(&path());
        let colored = canonicalize(&path(), |_| 0u32, |bond| (bond == b(1)) as u32);
        assert!(plain != colored);
    }

    #[test]
    fn equal_canonicals_hash_equally() {
        assert_eq!(hash_of(&canon(&path())), hash_of(&canon(&shifted_path())));
    }

    #[test]
    fn ordering_agrees_with_equality() {
        let a = canon(&path());
        let b = canon(&shifted_path());
        let c = canon(&triangle());
        assert_eq!(a.cmp(&b), Ordering::Equal);
        assert_ne!(a.cmp(&c), Ordering::Equal);
    }

    #[test]
    fn canonicalization_is_independent_of_input_order() {
        let plain = canon(&path());
        let reordered = canon(&reordered_path());
        assert!(plain == reordered);
        let orbits = |c: &Canonical<u32, u32>| -> Vec<Vec<SiteId>> {
            c.orbits().map(<[SiteId]>::to_vec).collect()
        };
        assert_eq!(orbits(&plain), orbits(&reordered));
    }
}

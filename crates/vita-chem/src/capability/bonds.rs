use vita_core::{HasSites, SiteId};

use crate::BondId;

/// The bonding skeleton: the chemical bonds connecting sites.
///
/// `HasBonds` is the supertrait of all per-bond capabilities, just as
/// [`HasSites`] is the supertrait of all per-site capabilities. A type
/// cannot expose data *about* bonds without first declaring *which* bonds
/// exist and which sites they connect.
///
/// # Contract
///
/// - [`bonds`](HasBonds::bonds) yields each [`BondId`] exactly once, with no
///   duplicates.
/// - [`bond_endpoints`](HasBonds::bond_endpoints) is total over
///   [`bonds`](HasBonds::bonds): every bond has exactly two endpoints, both
///   of which are in [`sites`](HasSites::sites).
/// - At most one bond exists per unordered site pair (simple graph).
pub trait HasBonds: HasSites {
    /// Returns an iterator over the identifier of every bond.
    fn bonds(&self) -> impl Iterator<Item = BondId> + '_;

    /// Returns the two sites connected by `bond`.
    ///
    /// # Panics
    ///
    /// Panics if `bond` is not in [`bonds`](HasBonds::bonds).
    fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId);

    /// Returns the number of bonds.
    ///
    /// The default implementation consumes [`bonds`](HasBonds::bonds);
    /// override it when the count is known in `O(1)`.
    #[inline]
    fn bond_count(&self) -> usize {
        self.bonds().count()
    }

    /// Returns whether `bond` is in [`bonds`](HasBonds::bonds).
    ///
    /// The default implementation scans [`bonds`](HasBonds::bonds); override
    /// it when membership can be decided in better than `O(n)`.
    #[inline]
    fn contains_bond(&self, bond: BondId) -> bool {
        self.bonds().any(|b| b == bond)
    }

    /// Returns the bond connecting `a` and `b`, if one exists.
    ///
    /// Returns `None` when no bond connects the two sites. The check is
    /// symmetric: `bond_between(a, b) == bond_between(b, a)`. The default
    /// implementation scans [`bonds`](HasBonds::bonds) in `O(n)`; override
    /// it for `O(1)` lookup via an adjacency map.
    #[inline]
    fn bond_between(&self, a: SiteId, b: SiteId) -> Option<BondId> {
        self.bonds().find(|&bond| {
            let (u, v) = self.bond_endpoints(bond);
            (u == a && v == b) || (u == b && v == a)
        })
    }

    /// Returns an iterator over `(bond, other_site)` pairs for every bond
    /// incident to `site`.
    ///
    /// `other_site` is the far endpoint of the bond from `site`. The default
    /// implementation scans all bonds in `O(n)`; override it with an
    /// adjacency list for `O(degree)`.
    #[inline]
    fn bonds_of(&self, site: SiteId) -> impl Iterator<Item = (BondId, SiteId)> + '_ {
        self.bonds().filter_map(move |bond| {
            let (a, b) = self.bond_endpoints(bond);
            if a == site {
                Some((bond, b))
            } else if b == site {
                Some((bond, a))
            } else {
                None
            }
        })
    }

    /// Returns an iterator over the site identifiers of every neighbour of
    /// `site`.
    ///
    /// The default implementation delegates to
    /// [`bonds_of`](HasBonds::bonds_of).
    #[inline]
    fn neighbors(&self, site: SiteId) -> impl Iterator<Item = SiteId> + '_ {
        self.bonds_of(site).map(|(_, nb)| nb)
    }

    /// Returns the degree (bond count) of `site`.
    ///
    /// The default implementation counts [`bonds_of`](HasBonds::bonds_of);
    /// override it when the degree is known in `O(1)`.
    #[inline]
    fn degree(&self, site: SiteId) -> usize {
        self.bonds_of(site).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn bond(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Bare {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasBonds for Bare {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, b: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == b).unwrap();
            self.endpoints[i]
        }
    }

    struct Adjacency {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
    }
    impl HasSites for Adjacency {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasBonds for Adjacency {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, b: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == b).unwrap();
            self.endpoints[i]
        }

        fn bond_count(&self) -> usize {
            self.bonds.len()
        }

        fn contains_bond(&self, b: BondId) -> bool {
            self.bonds.contains(&b)
        }

        fn bond_between(&self, a: SiteId, b: SiteId) -> Option<BondId> {
            self.bonds
                .iter()
                .copied()
                .zip(self.endpoints.iter().copied())
                .find_map(|(bond, (u, v))| {
                    if (u == a && v == b) || (u == b && v == a) {
                        Some(bond)
                    } else {
                        None
                    }
                })
        }
    }

    fn chain() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            bonds: vec![bond(1), bond(2)],
            endpoints: vec![(site(1), site(2)), (site(2), site(3))],
        }
    }

    #[test]
    fn bonds() {
        let mol = chain();
        assert_eq!(mol.bonds().collect::<Vec<_>>(), vec![bond(1), bond(2)]);
    }

    #[test]
    fn bonds_empty() {
        let mol = Bare {
            sites: vec![],
            bonds: vec![],
            endpoints: vec![],
        };
        assert_eq!(mol.bonds().count(), 0);
    }

    #[test]
    fn bond_endpoints() {
        let mol = chain();
        assert_eq!(mol.bond_endpoints(bond(1)), (site(1), site(2)));
        assert_eq!(mol.bond_endpoints(bond(2)), (site(2), site(3)));
    }

    #[test]
    fn bond_count() {
        let mol = chain();
        assert_eq!(mol.bond_count(), 2);
    }

    #[test]
    fn bond_count_empty_is_zero() {
        let mol = Bare {
            sites: vec![],
            bonds: vec![],
            endpoints: vec![],
        };
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn contains_bond() {
        let mol = chain();
        assert!(mol.contains_bond(bond(1)));
        assert!(mol.contains_bond(bond(2)));
        assert!(!mol.contains_bond(bond(99)));
    }

    #[test]
    fn bond_between_present() {
        let mol = chain();
        assert_eq!(mol.bond_between(site(1), site(2)), Some(bond(1)));
        assert_eq!(mol.bond_between(site(2), site(3)), Some(bond(2)));
    }

    #[test]
    fn bond_between_symmetric() {
        let mol = chain();
        assert_eq!(
            mol.bond_between(site(1), site(2)),
            mol.bond_between(site(2), site(1))
        );
        assert_eq!(
            mol.bond_between(site(2), site(3)),
            mol.bond_between(site(3), site(2))
        );
    }

    #[test]
    fn bond_between_absent() {
        let mol = chain();
        assert_eq!(mol.bond_between(site(1), site(3)), None);
    }

    #[test]
    fn bonds_of() {
        let mol = chain();
        assert_eq!(
            mol.bonds_of(site(1)).collect::<Vec<_>>(),
            vec![(bond(1), site(2))]
        );
        assert_eq!(
            mol.bonds_of(site(2)).collect::<Vec<_>>(),
            vec![(bond(1), site(1)), (bond(2), site(3))]
        );
        assert_eq!(
            mol.bonds_of(site(3)).collect::<Vec<_>>(),
            vec![(bond(2), site(2))]
        );
    }

    #[test]
    fn bonds_of_isolated_site() {
        let mol = Bare {
            sites: vec![site(1), site(2)],
            bonds: vec![],
            endpoints: vec![],
        };
        assert_eq!(mol.bonds_of(site(1)).count(), 0);
    }

    #[test]
    fn neighbors() {
        let mol = chain();
        assert_eq!(mol.neighbors(site(1)).collect::<Vec<_>>(), vec![site(2)]);
        assert_eq!(
            mol.neighbors(site(2)).collect::<Vec<_>>(),
            vec![site(1), site(3)]
        );
        assert_eq!(mol.neighbors(site(3)).collect::<Vec<_>>(), vec![site(2)]);
    }

    #[test]
    fn degree() {
        let mol = chain();
        assert_eq!(mol.degree(site(1)), 1);
        assert_eq!(mol.degree(site(2)), 2);
        assert_eq!(mol.degree(site(3)), 1);
    }

    #[test]
    fn degree_isolated_site() {
        let mol = Bare {
            sites: vec![site(1)],
            bonds: vec![],
            endpoints: vec![],
        };
        assert_eq!(mol.degree(site(1)), 0);
    }

    #[test]
    fn override_matches_default() {
        use std::collections::BTreeMap;

        let sites = vec![site(1), site(2), site(3)];
        let bonds_vec = vec![bond(1), bond(2)];
        let endpoints_vec = vec![(site(1), site(2)), (site(2), site(3))];

        let bare = Bare {
            sites: sites.clone(),
            bonds: bonds_vec.clone(),
            endpoints: endpoints_vec.clone(),
        };
        let adj = Adjacency {
            sites,
            bonds: bonds_vec,
            endpoints: endpoints_vec,
        };

        assert_eq!(bare.bond_count(), adj.bond_count());

        for b in [bond(1), bond(2), bond(99)] {
            assert_eq!(bare.contains_bond(b), adj.contains_bond(b));
        }

        let all_sites = [site(1), site(2), site(3)];
        for &a in &all_sites {
            for &b in &all_sites {
                assert_eq!(bare.bond_between(a, b), adj.bond_between(a, b));
            }
        }

        let bare_adj: BTreeMap<SiteId, Vec<(BondId, SiteId)>> = all_sites
            .iter()
            .map(|&s| {
                let mut inc: Vec<_> = bare.bonds_of(s).collect();
                inc.sort();
                (s, inc)
            })
            .collect();
        let adj_adj: BTreeMap<SiteId, Vec<(BondId, SiteId)>> = all_sites
            .iter()
            .map(|&s| {
                let mut inc: Vec<_> = adj.bonds_of(s).collect();
                inc.sort();
                (s, inc)
            })
            .collect();
        assert_eq!(bare_adj, adj_adj);
    }
}

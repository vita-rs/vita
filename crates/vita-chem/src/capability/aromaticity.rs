use vita_core::SiteId;

use crate::{BondId, HasBonds};

/// Per-bond aromaticity: whether each bond lies in an aromatic π system.
///
/// Access is by keyed lookup:
/// [`is_aromatic`](HasAromaticity::is_aromatic) maps a [`BondId`] to whether
/// that bond is aromatic.
/// [`is_aromatic_site`](HasAromaticity::is_aromatic_site) reports the same for
/// a [`SiteId`], derived from the bonds incident to it.
///
/// # Contract
///
/// [`is_aromatic`](HasAromaticity::is_aromatic) is total over
/// [`bonds`](HasBonds::bonds): every bond is either aromatic or not.
pub trait HasAromaticity: HasBonds {
    /// Returns whether `bond` is aromatic.
    ///
    /// # Panics
    ///
    /// Panics if `bond` is not in [`bonds`](HasBonds::bonds).
    fn is_aromatic(&self, bond: BondId) -> bool;

    /// Returns whether `site` lies in an aromatic system.
    ///
    /// A site is aromatic exactly when one of its bonds is aromatic; in
    /// biphenyl the two sites joined by the central bond are thus aromatic
    /// while that bond is not. The default implementation scans the bonds
    /// incident to `site`; override it when site aromaticity can be determined
    /// directly.
    #[inline]
    fn is_aromatic_site(&self, site: SiteId) -> bool {
        self.bonds_of(site).any(|(bond, _)| self.is_aromatic(bond))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vita_core::HasSites;

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
        aromatic: Vec<bool>,
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

        fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.endpoints[i]
        }
    }
    impl HasAromaticity for Bare {
        fn is_aromatic(&self, bond: BondId) -> bool {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.aromatic[i]
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        aromatic: Vec<bool>,
        aromatic_sites: Vec<bool>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasBonds for Columnar {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.endpoints[i]
        }
    }
    impl HasAromaticity for Columnar {
        fn is_aromatic(&self, bond: BondId) -> bool {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.aromatic[i]
        }

        fn is_aromatic_site(&self, site: SiteId) -> bool {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.aromatic_sites[i]
        }
    }

    fn biphenyl() -> Bare {
        let mut aromatic = vec![true; 12];
        aromatic.push(false);
        Bare {
            sites: (1..=12).map(site).collect(),
            bonds: (1..=13).map(bond).collect(),
            endpoints: vec![
                (site(1), site(2)),
                (site(2), site(3)),
                (site(3), site(4)),
                (site(4), site(5)),
                (site(5), site(6)),
                (site(6), site(1)),
                (site(7), site(8)),
                (site(8), site(9)),
                (site(9), site(10)),
                (site(10), site(11)),
                (site(11), site(12)),
                (site(12), site(7)),
                (site(1), site(7)),
            ],
            aromatic,
        }
    }

    fn ethane() -> Bare {
        Bare {
            sites: vec![site(1), site(2)],
            bonds: vec![bond(1)],
            endpoints: vec![(site(1), site(2))],
            aromatic: vec![false],
        }
    }

    #[test]
    fn is_aromatic() {
        let mol = biphenyl();
        assert!(mol.is_aromatic(bond(1)));
        assert!(!mol.is_aromatic(bond(13)));
        assert!(!ethane().is_aromatic(bond(1)));
    }

    #[test]
    fn is_aromatic_site() {
        assert!(biphenyl().is_aromatic_site(site(2)));
        assert!(!ethane().is_aromatic_site(site(1)));
    }

    #[test]
    fn plain_bond_between_aromatic_sites() {
        let mol = biphenyl();
        assert!(!mol.is_aromatic(bond(13)));
        assert!(mol.is_aromatic_site(site(1)));
        assert!(mol.is_aromatic_site(site(7)));
    }

    #[test]
    fn override_matches_default() {
        let sites = vec![site(1), site(2), site(3)];
        let bonds = vec![bond(1), bond(2)];
        let endpoints = vec![(site(1), site(2)), (site(2), site(3))];
        let aromatic = vec![true, false];

        let bare = Bare {
            sites: sites.clone(),
            bonds: bonds.clone(),
            endpoints: endpoints.clone(),
            aromatic: aromatic.clone(),
        };
        let col = Columnar {
            sites: sites.clone(),
            bonds,
            endpoints,
            aromatic,
            aromatic_sites: vec![true, true, false],
        };

        for s in sites {
            assert_eq!(bare.is_aromatic_site(s), col.is_aromatic_site(s));
        }
    }
}

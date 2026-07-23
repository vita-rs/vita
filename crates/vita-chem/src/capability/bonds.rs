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

    /// Returns an iterator over the site identifiers of every neighbor of
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

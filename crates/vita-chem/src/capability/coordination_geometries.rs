use vita_core::{HasSites, SiteId};

use crate::CoordinationGeometry;

/// Per-site coordination: the [`CoordinationGeometry`] of each site's substituents.
///
/// Access is by keyed lookup:
/// [`coordination_geometry`](HasCoordinationGeometries::coordination_geometry) maps
/// a [`SiteId`] to the geometry its substituents take, if they take one.
/// [`coordination_geometries`](HasCoordinationGeometries::coordination_geometries)
/// yields one answer per site in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`coordination_geometry`](HasCoordinationGeometries::coordination_geometry) is
/// total over [`sites`](HasSites::sites): every site has exactly one answer, `None`
/// where its substituents take no named arrangement.
/// [`coordination_geometries`](HasCoordinationGeometries::coordination_geometries)
/// yields values in the same order as [`sites`](HasSites::sites).
pub trait HasCoordinationGeometries: HasSites {
    /// Returns the geometry `site`'s substituents take, or `None` if they take no
    /// named arrangement.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn coordination_geometry(&self, site: SiteId) -> Option<CoordinationGeometry>;

    /// Yields one answer per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up
    /// [`coordination_geometry`](HasCoordinationGeometries::coordination_geometry)
    /// per site; override it when the geometries can be produced directly.
    #[inline]
    fn coordination_geometries(&self) -> impl Iterator<Item = Option<CoordinationGeometry>> + '_ {
        self.sites()
            .map(move |site| self.coordination_geometry(site))
    }
}

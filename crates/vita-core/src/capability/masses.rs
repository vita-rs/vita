use crate::units::mass::{Mass, MassUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site mass: the [`Mass`] of the particle at each site.
///
/// Access is by keyed lookup: [`mass`](HasMasses::mass) maps a [`SiteId`] to its mass,
/// in any requested [unit](MassUnit). [`masses`](HasMasses::masses) yields one mass per
/// site in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`mass`](HasMasses::mass) is total over [`sites`](HasSites::sites): every site has
/// exactly one mass.
/// [`masses`](HasMasses::masses) yields values in the same order as
/// [`sites`](HasSites::sites).
pub trait HasMasses<V: Scalar>: HasSites {
    /// Returns the mass of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn mass<U: MassUnit>(&self, site: SiteId) -> Mass<V, U>;

    /// Yields one mass per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up [`mass`](HasMasses::mass) per site; override
    /// it when the masses can be produced directly.
    #[inline]
    fn masses<U: MassUnit>(&self) -> impl Iterator<Item = Mass<V, U>> + '_ {
        self.sites().map(move |site| self.mass::<U>(site))
    }
}

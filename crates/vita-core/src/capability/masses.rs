use crate::units::mass::{Mass, MassUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site mass: the [`Mass`] of the particle at each site.
///
/// Access is by lookup: [`mass`](HasMasses::mass) maps a [`SiteId`] to its mass, in any
/// requested [unit](MassUnit). [`masses`](HasMasses::masses) iterates every `(site, mass)`
/// pair.
///
/// # Contract
///
/// [`mass`](HasMasses::mass) is total over [`sites`](HasSites::sites): every site has
/// exactly one mass.
pub trait HasMasses<V: Scalar>: HasSites {
    /// Returns the mass of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn mass<U: MassUnit>(&self, site: SiteId) -> Mass<V, U>;

    /// Returns an iterator over every `(site, mass)` pair, each mass in unit `U`.
    ///
    /// Each mass is yielded with its [`SiteId`]. The default implementation looks up
    /// [`mass`](HasMasses::mass) per site; override it when the pairs can be produced
    /// directly.
    #[inline]
    fn masses<U: MassUnit>(&self) -> impl Iterator<Item = (SiteId, Mass<V, U>)> + '_ {
        self.sites().map(move |site| (site, self.mass::<U>(site)))
    }
}

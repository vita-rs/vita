use crate::tensor::Point3;
use crate::units::length::{Length, LengthUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site position: the [`Point3`] locating each site in space.
///
/// Access is by keyed lookup: [`position`](HasPositions::position) maps a [`SiteId`] to
/// its position, in any requested [unit](LengthUnit). [`positions`](HasPositions::positions)
/// yields one position per site in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`position`](HasPositions::position) is total over [`sites`](HasSites::sites): every
/// site has exactly one position.
/// [`positions`](HasPositions::positions) yields values in the same order as
/// [`sites`](HasSites::sites).
pub trait HasPositions<V: Scalar>: HasSites {
    /// Returns the position of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<V, U>>;

    /// Yields one position per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up [`position`](HasPositions::position) per site;
    /// override it when the positions can be produced directly.
    #[inline]
    fn positions<U: LengthUnit>(&self) -> impl Iterator<Item = Point3<Length<V, U>>> + '_ {
        self.sites().map(move |site| self.position::<U>(site))
    }
}

use crate::tensor::Point3;
use crate::units::length::{Length, LengthUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site position: the [`Point3`] locating each site in space.
///
/// Access is by lookup: [`position`](HasPositions::position) maps a [`SiteId`] to its
/// position, in any requested [unit](LengthUnit). [`positions`](HasPositions::positions)
/// iterates every `(site, position)` pair.
///
/// # Contract
///
/// [`position`](HasPositions::position) is total over [`sites`](HasSites::sites): every
/// site has exactly one position.
pub trait HasPositions<V: Scalar>: HasSites {
    /// Returns the position of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<V, U>>;

    /// Returns an iterator over every `(site, position)` pair, each position in unit `U`.
    ///
    /// Each position is yielded with its [`SiteId`]. The default implementation looks up
    /// [`position`](HasPositions::position) per site; override it when the pairs can be
    /// produced directly.
    #[inline]
    fn positions<U: LengthUnit>(
        &self,
    ) -> impl Iterator<Item = (SiteId, Point3<Length<V, U>>)> + '_ {
        self.sites()
            .map(move |site| (site, self.position::<U>(site)))
    }
}

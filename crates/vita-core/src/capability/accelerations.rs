use crate::tensor::Vector3;
use crate::units::acceleration::{Acceleration, AccelerationUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site acceleration: the [`Vector3`] acceleration of each site.
///
/// Access is by keyed lookup: [`acceleration`](HasAccelerations::acceleration) maps a
/// [`SiteId`] to its acceleration, in any requested [unit](AccelerationUnit).
/// [`accelerations`](HasAccelerations::accelerations) yields one acceleration per site
/// in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`acceleration`](HasAccelerations::acceleration) is total over [`sites`](HasSites::sites):
/// every site has exactly one acceleration.
/// [`accelerations`](HasAccelerations::accelerations) yields values in the same order as
/// [`sites`](HasSites::sites).
pub trait HasAccelerations<V: Scalar>: HasSites {
    /// Returns the acceleration of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn acceleration<U: AccelerationUnit>(&self, site: SiteId) -> Vector3<Acceleration<V, U>>;

    /// Yields one acceleration per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up
    /// [`acceleration`](HasAccelerations::acceleration) per site; override it when the
    /// accelerations can be produced directly.
    #[inline]
    fn accelerations<U: AccelerationUnit>(
        &self,
    ) -> impl Iterator<Item = Vector3<Acceleration<V, U>>> + '_ {
        self.sites().map(move |site| self.acceleration::<U>(site))
    }
}

use crate::tensor::Vector3;
use crate::units::velocity::{Velocity, VelocityUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site velocity: the [`Vector3`] velocity of each site.
///
/// Access is by keyed lookup: [`velocity`](HasVelocities::velocity) maps a [`SiteId`]
/// to its velocity, in any requested [unit](VelocityUnit).
/// [`velocities`](HasVelocities::velocities) yields one velocity per site in
/// [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`velocity`](HasVelocities::velocity) is total over [`sites`](HasSites::sites): every
/// site has exactly one velocity.
/// [`velocities`](HasVelocities::velocities) yields values in the same order as
/// [`sites`](HasSites::sites).
pub trait HasVelocities<V: Scalar>: HasSites {
    /// Returns the velocity of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn velocity<U: VelocityUnit>(&self, site: SiteId) -> Vector3<Velocity<V, U>>;

    /// Yields one velocity per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up [`velocity`](HasVelocities::velocity) per
    /// site; override it when the velocities can be produced directly.
    #[inline]
    fn velocities<U: VelocityUnit>(&self) -> impl Iterator<Item = Vector3<Velocity<V, U>>> + '_ {
        self.sites().map(move |site| self.velocity::<U>(site))
    }
}

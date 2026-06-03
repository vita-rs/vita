use crate::tensor::Vector3;
use crate::units::velocity::{Velocity, VelocityUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site velocity: the [`Vector3`] velocity of each site.
///
/// Access is by lookup: [`velocity`](HasVelocities::velocity) maps a [`SiteId`] to its
/// velocity, in any requested [unit](VelocityUnit). [`velocities`](HasVelocities::velocities)
/// iterates every `(site, velocity)` pair.
///
/// # Contract
///
/// [`velocity`](HasVelocities::velocity) is total over [`sites`](HasSites::sites): every
/// site has exactly one velocity.
pub trait HasVelocities<V: Scalar>: HasSites {
    /// Returns the velocity of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn velocity<U: VelocityUnit>(&self, site: SiteId) -> Vector3<Velocity<V, U>>;

    /// Returns an iterator over every `(site, velocity)` pair, each velocity in unit `U`.
    ///
    /// Each velocity is yielded with its [`SiteId`]. The default implementation looks up
    /// [`velocity`](HasVelocities::velocity) per site; override it when the pairs can be
    /// produced directly.
    #[inline]
    fn velocities<U: VelocityUnit>(
        &self,
    ) -> impl Iterator<Item = (SiteId, Vector3<Velocity<V, U>>)> + '_ {
        self.sites()
            .map(move |site| (site, self.velocity::<U>(site)))
    }
}

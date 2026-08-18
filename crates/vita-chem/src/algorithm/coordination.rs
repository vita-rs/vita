//! The arrangement a site's substituents take about it, read two ways.
//!
//! [`perceive`] predicts it from the Lewis structure within VSEPR's model — the
//! electron domains settle where they separate furthest, and the substituents take
//! whichever of that parent arrangement's vertices the lone pairs leave. [`departures`]
//! instead measures it, reporting how far a site's observed directions lie from every
//! idealized arrangement its coordination number admits. One answers within a model and
//! names a geometry; the other reads coordinates and names none, because the nearest
//! arrangement is a fact about distances while *near enough* is a judgement the caller
//! owns.

mod departures;
mod perceive;

pub use departures::{Departures, departures};
pub use perceive::perceive;

use vita_core::{HasSites, SiteId};

use crate::CoordinationGeometry;
use crate::algorithm::utils::SortedMap;
use crate::capability::HasCoordinationGeometries;
use crate::capability::delegation::forward_capabilities;

/// The [`CoordinationGeometry`] of each of a molecule's sites.
///
/// A site whose substituents take no named arrangement — fewer than two of them, more
/// than the vocabulary reaches, or a count the reading could not settle — is absent.
///
/// Obtain via [`perceive`] or [`Departures::geometries`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoordinationGeometries {
    geometries: SortedMap<SiteId, CoordinationGeometry>,
}

impl CoordinationGeometries {
    /// Number of sites with a geometry.
    pub fn len(&self) -> usize {
        self.geometries.len()
    }

    /// Returns `true` if no site has one.
    pub fn is_empty(&self) -> bool {
        self.geometries.is_empty()
    }

    /// The geometry `site`'s substituents take, or `None` if they take none.
    pub fn get(&self, site: SiteId) -> Option<CoordinationGeometry> {
        self.geometries.get(&site).copied()
    }

    /// Iterates `(site, geometry)` pairs in ascending site order.
    pub fn iter(&self) -> impl Iterator<Item = (SiteId, CoordinationGeometry)> + '_ {
        self.geometries
            .iter()
            .map(|(&site, &geometry)| (site, geometry))
    }

    /// Binds these geometries to `mol`, yielding a view that implements
    /// [`HasCoordinationGeometries`].
    ///
    /// The view borrows both, so `mol` stays immutable while it is held — the geometries
    /// cannot silently fall out of step with the molecule they describe. Use it to feed a
    /// molecule's computed geometries to anything that reads the
    /// [`HasCoordinationGeometries`] capability.
    pub fn bind<'a, M: HasSites>(&'a self, mol: &'a M) -> WithCoordinationGeometries<'a, M> {
        WithCoordinationGeometries {
            mol,
            geometries: self,
        }
    }

    /// The `geometries` gathered into a site-ordered set.
    fn from_pairs(geometries: impl IntoIterator<Item = (SiteId, CoordinationGeometry)>) -> Self {
        CoordinationGeometries {
            geometries: SortedMap::from_pairs(geometries),
        }
    }
}

/// A molecule viewed together with a set of [`CoordinationGeometries`].
///
/// Answers the coordination geometries from that set and forwards every other core and
/// chem capability to the molecule, so a computed result reads as the
/// [`HasCoordinationGeometries`] capability its consumers expect — at no cost beyond the
/// two references it holds.
///
/// Obtain via [`CoordinationGeometries::bind`].
pub struct WithCoordinationGeometries<'a, M> {
    mol: &'a M,
    geometries: &'a CoordinationGeometries,
}

impl<M> Copy for WithCoordinationGeometries<'_, M> {}

impl<M> Clone for WithCoordinationGeometries<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

forward_capabilities!(
    WithCoordinationGeometries,
    mol,
    HasAccelerations,
    HasElements,
    HasIsotopes,
    HasLattice,
    HasMasses,
    HasNetCharge,
    HasPositions,
    HasSites,
    HasVelocities,
    HasAromaticity,
    HasBondOrders,
    HasBonds,
    HasFormalCharges,
    HasPartialCharges,
    HasRadicalElectrons,
    HasStereoConfigurations,
);

impl<M: HasSites> HasCoordinationGeometries for WithCoordinationGeometries<'_, M> {
    fn coordination_geometry(&self, site: SiteId) -> Option<CoordinationGeometry> {
        assert!(
            self.mol.contains_site(site),
            "site is not in the molecule the geometries are bound to"
        );
        self.geometries.get(site)
    }

    fn coordination_geometries(&self) -> impl Iterator<Item = Option<CoordinationGeometry>> + '_ {
        self.mol.sites().map(|site| self.geometries.get(site))
    }
}

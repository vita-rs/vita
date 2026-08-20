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

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasElements};

    use crate::CoordinationGeometry::{Angular, Tetrahedral};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    struct Mol;

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            (1..=3).map(s)
        }
    }

    impl HasElements for Mol {
        fn element(&self, _site: SiteId) -> Element {
            Element::from_symbol("C").unwrap()
        }
    }

    fn geometries() -> CoordinationGeometries {
        CoordinationGeometries::from_pairs([(s(3), Angular), (s(1), Tetrahedral)])
    }

    #[test]
    fn a_set_of_no_geometries_is_empty() {
        let none = CoordinationGeometries::from_pairs([]);
        assert_eq!(none.len(), 0);
        assert!(none.is_empty());
    }

    #[test]
    fn a_site_answers_the_geometry_it_was_given() {
        assert_eq!(geometries().get(s(1)), Some(Tetrahedral));
        assert_eq!(geometries().get(s(3)), Some(Angular));
    }

    #[test]
    fn a_site_that_was_given_none_answers_none() {
        assert_eq!(geometries().get(s(2)), None);
    }

    #[test]
    fn the_geometries_count_the_sites_that_have_one() {
        assert_eq!(geometries().len(), 2);
        assert!(!geometries().is_empty());
    }

    #[test]
    fn the_geometries_come_out_in_ascending_site_order() {
        assert_eq!(
            geometries().iter().collect::<Vec<_>>(),
            vec![(s(1), Tetrahedral), (s(3), Angular)]
        );
    }

    #[test]
    fn a_bound_view_answers_the_capability_from_the_geometries() {
        let geometries = geometries();
        let view = geometries.bind(&Mol);
        assert_eq!(view.coordination_geometry(s(1)), Some(Tetrahedral));
        assert_eq!(view.coordination_geometry(s(2)), None);
    }

    #[test]
    fn a_bound_view_yields_one_answer_per_site() {
        let geometries = geometries();
        let view = geometries.bind(&Mol);
        assert_eq!(
            view.coordination_geometries().collect::<Vec<_>>(),
            vec![Some(Tetrahedral), None, Some(Angular)]
        );
    }

    #[test]
    fn a_bound_view_forwards_the_skeleton() {
        let geometries = geometries();
        let view = geometries.bind(&Mol);
        assert_eq!(view.site_count(), Mol.site_count());
        assert_eq!(view.element(s(1)), Mol.element(s(1)));
    }

    #[test]
    #[should_panic(expected = "site is not in the molecule")]
    fn a_bound_view_refuses_a_site_the_molecule_does_not_hold() {
        let geometries = geometries();
        geometries.bind(&Mol).coordination_geometry(s(9));
    }
}

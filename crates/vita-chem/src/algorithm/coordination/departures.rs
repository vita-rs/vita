use vita_core::tensor::{Matrix3, Vector3};
use vita_core::units::angle::{Angle, AngleUnit, Radian};
use vita_core::units::length::Angstrom;
use vita_core::{HasPositions, Quantity, Scalar, SiteId};

use super::CoordinationGeometries;
use crate::algorithm::utils::{SortedMap, next_permutation};
use crate::{CoordinationGeometry, HasBonds};

/// How far each of a molecule's sites departs from every arrangement its coordination
/// number admits.
///
/// A departure is an angle: zero where the substituents sit exactly on an idealized
/// arrangement's directions, and otherwise the angle whose chord is the root-mean-square
/// chord between each substituent and the slot it fills, taken over the best rotation and
/// the best assignment of substituents to slots. Only directions enter, so a site is
/// measured by the angles it holds its substituents at and not by how far away it holds
/// them.
///
/// Every candidate is measured, not just the nearest, because the standard names the
/// nearest idealized geometry but nowhere says how near is near enough: a site 3° from
/// one arrangement and 40° from the rest is a different report from one lying 20° from
/// two of them, and only the whole ranking tells them apart.
///
/// Obtain via [`departures`].
#[derive(Clone, Debug, PartialEq)]
pub struct Departures<V> {
    departures: SortedMap<SiteId, Vec<(CoordinationGeometry, V)>>,
}

impl<V: Scalar> Departures<V> {
    /// Number of sites measured.
    pub fn len(&self) -> usize {
        self.departures.len()
    }

    /// Returns `true` if no site was measured.
    pub fn is_empty(&self) -> bool {
        self.departures.is_empty()
    }

    /// The angle by which `site` departs from `geometry`, or `None` if the site was not
    /// measured or the geometry does not hold its substituents.
    pub fn get<U: AngleUnit>(
        &self,
        site: SiteId,
        geometry: CoordinationGeometry,
    ) -> Option<Angle<V, U>> {
        self.departures
            .get(&site)?
            .iter()
            .find(|&&(candidate, _)| candidate == geometry)
            .map(|&(_, departure)| Angle::<V, Radian>::new(departure).to::<U>())
    }

    /// Iterates the arrangements `site` was measured against with its departure from
    /// each, nearest first and ties in ascending [`CoordinationGeometry`] order.
    pub fn of_site<U: AngleUnit>(
        &self,
        site: SiteId,
    ) -> impl Iterator<Item = (CoordinationGeometry, Angle<V, U>)> + '_ {
        self.departures
            .get(&site)
            .into_iter()
            .flatten()
            .map(|&(geometry, departure)| (geometry, Angle::<V, Radian>::new(departure).to::<U>()))
    }

    /// The arrangement `site` departs least from — the first of
    /// [`of_site`](Departures::of_site) — or `None` if it was not measured.
    ///
    /// Nearest is a fact about distances; whether it is near enough to call the site that
    /// shape is not, and this does not decide it — read the departure alongside.
    pub fn nearest(&self, site: SiteId) -> Option<CoordinationGeometry> {
        self.departures
            .get(&site)?
            .first()
            .map(|&(geometry, _)| geometry)
    }

    /// The nearest arrangement of every measured site, as the geometries a molecule can
    /// be read through.
    pub fn geometries(&self) -> CoordinationGeometries {
        CoordinationGeometries::from_pairs(
            self.departures
                .iter()
                .filter_map(|(&site, ranked)| Some((site, ranked.first()?.0))),
        )
    }
}

/// Measures how far each of a molecule's sites departs from every arrangement its
/// coordination number admits.
///
/// Each site is read as the unit directions to its substituents, and each candidate
/// arrangement as the unit directions of its own idealized slots. The two are brought
/// together over every assignment of substituents to slots and, for each, the rotation
/// that fits them best; the least residual left over is the site's departure from that
/// arrangement. A site bearing fewer than two substituents, or more than the vocabulary
/// reaches, is left unmeasured — one direction fixes no arrangement — as is one whose
/// substituent sits on top of it.
///
/// # Complexity
///
/// O(V · n · n! + V · log V) time and O(V) space, over the molecule's `V` sites, for a
/// largest coordination number `n`, assuming [`neighbors`](HasBonds::neighbors) runs in
/// O(degree); every assignment of every candidate is tried and each is scored
/// over `n` directions, `n` is capped at six, and the log factor orders the sites.
pub fn departures<M, V>(mol: &M) -> Departures<V>
where
    M: HasBonds + HasPositions<V>,
    V: Scalar,
{
    let position = |site: SiteId| mol.position::<Angstrom>(site).map(|length| length.value());
    Departures {
        departures: SortedMap::from_pairs(mol.sites().filter_map(|site| {
            let center = position(site);
            let observed: Vec<Vector3<V>> = mol
                .neighbors(site)
                .map(|neighbor| (position(neighbor) - center).try_normalize())
                .collect::<Option<_>>()?;
            Some((site, ranked(&observed)?))
        })),
    }
}

/// Every arrangement holding as many slots as there are `observed` directions, paired
/// with the departure from it and ordered nearest first, or `None` if the vocabulary
/// holds no arrangement of that size.
fn ranked<V: Scalar>(observed: &[Vector3<V>]) -> Option<Vec<(CoordinationGeometry, V)>> {
    let mut ranked: Vec<(CoordinationGeometry, V)> = CoordinationGeometry::ALL
        .into_iter()
        .filter(|geometry| geometry.slot_count() == observed.len())
        .map(|geometry| (geometry, departure(observed, geometry)))
        .collect();
    ranked.sort_by(|(_, near), (_, far)| {
        near.partial_cmp(far)
            .expect("a departure is a finite angle")
    });
    (!ranked.is_empty()).then_some(ranked)
}

/// The angle by which the `observed` directions depart from `geometry`: the angle whose
/// chord is the root-mean-square chord left by the best assignment of directions to slots
/// under the rotation that fits it best.
fn departure<V: Scalar>(observed: &[Vector3<V>], geometry: CoordinationGeometry) -> V {
    let reference: Vec<Vector3<V>> = geometry
        .directions()
        .iter()
        .map(|&[x, y, z]| {
            Vector3::from_array([V::from_f64(x), V::from_f64(y), V::from_f64(z)]).normalize()
        })
        .collect();

    let mut assignment: Vec<u8> = (0..observed.len() as u8).collect();
    let mut least = residual(observed, &reference, &assignment);
    while next_permutation(&mut assignment) {
        let candidate = residual(observed, &reference, &assignment);
        if candidate < least {
            least = candidate;
        }
    }

    let mean_square = least / V::from_f64(observed.len() as f64);
    let half_chord = if mean_square > V::ZERO {
        mean_square.sqrt() / V::from_f64(2.0)
    } else {
        V::ZERO
    };
    let sine = if half_chord < V::ONE {
        half_chord
    } else {
        V::ONE
    };
    sine.asin() * V::from_f64(2.0)
}

/// The squared chord the `observed` directions leave against the `reference` slots they
/// fill under `assignment`, minimized over rotations.
///
/// The rotation drops out of the arithmetic: the greatest a rotation can bring the two
/// sets into register is the sum of the cross-covariance's singular values, the least of
/// them entering negated where the two sets are related by a reflection rather than a
/// rotation.
fn residual<V: Scalar>(observed: &[Vector3<V>], reference: &[Vector3<V>], assignment: &[u8]) -> V {
    let mut covariance = Matrix3::ZERO;
    for (slot, &filled_by) in assignment.iter().enumerate() {
        covariance += Matrix3::outer_product(observed[filled_by as usize], reference[slot]);
    }
    let values = covariance.singular_values();
    let least = if covariance.determinant() < V::ZERO {
        -values.z
    } else {
        values.z
    };
    (V::from_f64(observed.len() as f64) - values.x - values.y - least) * V::from_f64(2.0)
}

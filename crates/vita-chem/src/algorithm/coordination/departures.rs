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

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;
    use vita_core::tensor::Point3;
    use vita_core::units::angle::Degree;
    use vita_core::units::length::{Length, LengthUnit};

    use crate::BondId;

    const TOL: f64 = 1e-7;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        coords: Vec<[f64; 3]>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasPositions<f64> for Mol {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<f64, U>> {
            let [x, y, z] = self.coords[self.sites.iter().position(|&i| i == site).unwrap()];
            Point3::new(Length::new(x), Length::new(y), Length::new(z))
        }
    }

    impl HasBonds for Mol {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
            self.endpoints[self.bonds.iter().position(|&i| i == bond).unwrap()]
        }
    }

    fn centered(placements: &[[f64; 3]]) -> Mol {
        let mut coords = vec![[0.0, 0.0, 0.0]];
        coords.extend_from_slice(placements);
        Mol {
            sites: (1..=placements.len() as u32 + 1).map(s).collect(),
            coords,
            bonds: (1..=placements.len() as u32).map(b).collect(),
            endpoints: (1..=placements.len() as u32)
                .map(|ligand| (s(1), s(ligand + 1)))
                .collect(),
        }
    }

    fn ideal(geometry: CoordinationGeometry) -> Mol {
        centered(geometry.directions())
    }

    fn departure(site: &Mol, geometry: CoordinationGeometry) -> f64 {
        departures(site)
            .get::<Radian>(s(1), geometry)
            .unwrap()
            .value()
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL
    }

    fn mapped(placements: &[[f64; 3]], f: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<[f64; 3]> {
        placements.iter().map(|&point| f(point)).collect()
    }

    #[test]
    fn a_site_no_arrangement_fits_goes_unmeasured() {
        assert!(departures::<_, f64>(&centered(&[[1.0, 0.0, 0.0]])).is_empty());
        assert!(departures::<_, f64>(&centered(&[])).is_empty());
    }

    #[test]
    fn a_substituent_sitting_on_its_site_leaves_it_unmeasured() {
        let coincident = centered(&[[0.0, 0.0, 0.0], [0.0, 0.0, -1.0]]);
        assert_eq!(
            departures(&coincident).get::<Radian>(s(1), CoordinationGeometry::Linear),
            None
        );
    }

    #[test]
    fn a_site_on_an_arrangements_own_directions_departs_from_it_by_nothing() {
        for geometry in CoordinationGeometry::ALL {
            assert!(
                close(departure(&ideal(geometry), geometry), 0.0),
                "{geometry:?}"
            );
        }
    }

    #[test]
    fn bending_a_straight_site_departs_by_half_the_bend() {
        let bend = 0.4;
        let (sine, cosine) = bend.sin_cos();
        let bent = centered(&[[0.0, 0.0, 1.0], [sine, 0.0, -cosine]]);
        assert!(close(
            departure(&bent, CoordinationGeometry::Linear),
            bend / 2.0
        ));
    }

    #[test]
    fn a_departure_is_given_in_the_requested_unit() {
        let bent = centered(&[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]);
        let in_radians = departure(&bent, CoordinationGeometry::Linear);
        let in_degrees = departures(&bent)
            .get::<Degree>(s(1), CoordinationGeometry::Linear)
            .unwrap();
        assert!(close(in_degrees.value(), in_radians.to_degrees()));
    }

    #[test]
    fn every_other_candidate_of_the_same_size_lies_further() {
        for geometry in CoordinationGeometry::ALL {
            let site = ideal(geometry);
            for (candidate, apart) in departures(&site).of_site::<Radian>(s(1)) {
                assert_eq!(
                    candidate == geometry,
                    close(apart.value(), 0.0),
                    "{geometry:?} against {candidate:?}"
                );
            }
        }
    }

    #[test]
    fn a_bent_site_is_nearer_angular_than_linear() {
        let water = centered(&[[0.0, 0.0, 1.0], [0.968, 0.0, -0.25]]);
        assert_eq!(
            departures(&water).nearest(s(1)),
            Some(CoordinationGeometry::Angular)
        );
    }

    #[test]
    fn only_arrangements_of_the_right_size_are_measured() {
        let measured: Vec<CoordinationGeometry> =
            departures(&ideal(CoordinationGeometry::Octahedral))
                .of_site::<Radian>(s(1))
                .map(|(geometry, _)| geometry)
                .collect();
        assert_eq!(
            measured.len(),
            CoordinationGeometry::ALL
                .into_iter()
                .filter(|candidate| candidate.slot_count() == 6)
                .count()
        );
        assert!(measured.iter().all(|candidate| candidate.slot_count() == 6));
    }

    #[test]
    fn an_arrangement_of_another_size_answers_nothing() {
        let measured = departures(&ideal(CoordinationGeometry::Tetrahedral));
        assert_eq!(
            measured.get::<Radian>(s(1), CoordinationGeometry::Octahedral),
            None
        );
    }

    #[test]
    fn an_unmeasured_site_answers_nothing() {
        let measured = departures(&ideal(CoordinationGeometry::Tetrahedral));
        assert_eq!(measured.nearest(s(2)), None);
        assert_eq!(
            measured.get::<Radian>(s(2), CoordinationGeometry::Linear),
            None
        );
        assert_eq!(measured.of_site::<Radian>(s(2)).count(), 0);
    }

    #[test]
    fn an_arrangement_is_the_nearest_to_its_own_directions() {
        for geometry in CoordinationGeometry::ALL {
            assert_eq!(
                departures(&ideal(geometry)).nearest(s(1)),
                Some(geometry),
                "{geometry:?}"
            );
        }
    }

    #[test]
    fn the_candidates_come_out_nearest_first() {
        let apart: Vec<f64> = departures(&ideal(CoordinationGeometry::Tetrahedral))
            .of_site::<Radian>(s(1))
            .map(|(_, departure)| departure.value())
            .collect();
        assert_eq!(apart.len(), 4);
        assert!(apart.is_sorted_by(|near, far| near <= far));
    }

    #[test]
    fn the_geometries_are_the_nearest_arrangement_of_every_measured_site() {
        let measured = departures(&ideal(CoordinationGeometry::Tetrahedral));
        let geometries = measured.geometries();
        assert_eq!(geometries.len(), measured.len());
        assert_eq!(geometries.get(s(1)), measured.nearest(s(1)));
    }

    #[test]
    fn how_far_a_substituent_lies_does_not_move_the_departure() {
        for geometry in CoordinationGeometry::ALL {
            let stretched = centered(&mapped(geometry.directions(), |point| {
                point.map(|value| value * 3.0)
            }));
            let mut uneven = geometry.directions().to_vec();
            uneven[0] = uneven[0].map(|value| value * 7.0);
            assert!(close(departure(&stretched, geometry), 0.0), "{geometry:?}");
            assert!(
                close(departure(&centered(&uneven), geometry), 0.0),
                "{geometry:?}"
            );
        }
    }

    #[test]
    fn turning_a_site_does_not_move_its_departures() {
        let (sine, cosine) = 0.7f64.sin_cos();
        for geometry in CoordinationGeometry::ALL {
            let turned = centered(&mapped(geometry.directions(), |[x, y, z]| {
                [x * cosine - y * sine, x * sine + y * cosine, z]
            }));
            assert!(close(departure(&turned, geometry), 0.0), "{geometry:?}");
        }
    }

    #[test]
    fn reflecting_a_site_does_not_move_its_departures() {
        for geometry in CoordinationGeometry::ALL {
            let mirrored = centered(&mapped(geometry.directions(), |[x, y, z]| [-x, y, z]));
            assert!(close(departure(&mirrored, geometry), 0.0), "{geometry:?}");
        }
    }

    #[test]
    fn listing_the_substituents_in_another_order_does_not_move_the_departure() {
        for geometry in CoordinationGeometry::ALL {
            let mut reordered = geometry.directions().to_vec();
            reordered.reverse();
            assert!(
                close(departure(&centered(&reordered), geometry), 0.0),
                "{geometry:?}"
            );
        }
    }
}

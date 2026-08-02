use core::ops::{Add, AddAssign, Index, IndexMut, Sub, SubAssign};

use super::Vector3;
use crate::{Quantity, Scalar};

/// A point with three coordinates `x`, `y`, and `z`.
///
/// A point is a location, not a displacement: it translates by a [`Vector3`]
/// and subtracts from another point, but it does not scale or negate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Point3<T> {
    /// The first coordinate.
    pub x: T,
    /// The second coordinate.
    pub y: T,
    /// The third coordinate.
    pub z: T,
}

impl<T> Point3<T> {
    /// Constructs a point from its three coordinates.
    #[inline]
    pub const fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    /// Constructs a point from an array `[x, y, z]`.
    #[inline]
    pub fn from_array(array: [T; 3]) -> Self {
        let [x, y, z] = array;
        Self { x, y, z }
    }

    /// Returns the coordinates as an array `[x, y, z]`.
    #[inline]
    pub fn to_array(self) -> [T; 3] {
        [self.x, self.y, self.z]
    }

    /// Returns a copy of `self` with the `x` coordinate replaced.
    #[inline]
    pub fn with_x(self, x: T) -> Self {
        Self {
            x,
            y: self.y,
            z: self.z,
        }
    }

    /// Returns a copy of `self` with the `y` coordinate replaced.
    #[inline]
    pub fn with_y(self, y: T) -> Self {
        Self {
            x: self.x,
            y,
            z: self.z,
        }
    }

    /// Returns a copy of `self` with the `z` coordinate replaced.
    #[inline]
    pub fn with_z(self, z: T) -> Self {
        Self {
            x: self.x,
            y: self.y,
            z,
        }
    }

    /// Reinterprets a displacement from the origin as a point.
    #[inline]
    pub fn from_vector(vector: Vector3<T>) -> Self {
        Self {
            x: vector.x,
            y: vector.y,
            z: vector.z,
        }
    }

    /// Returns the position vector of `self` relative to the origin.
    #[inline]
    pub fn to_vector(self) -> Vector3<T> {
        Vector3::new(self.x, self.y, self.z)
    }

    /// Applies `f` to every coordinate, returning the resulting point.
    #[inline]
    pub fn map<U, F: FnMut(T) -> U>(self, mut f: F) -> Point3<U> {
        Point3 {
            x: f(self.x),
            y: f(self.y),
            z: f(self.z),
        }
    }

    /// Combines `self` and `rhs` coordinate-wise through `f`.
    #[inline]
    pub fn zip_map<U, R, F: FnMut(T, U) -> R>(self, rhs: Point3<U>, mut f: F) -> Point3<R> {
        Point3 {
            x: f(self.x, rhs.x),
            y: f(self.y, rhs.y),
            z: f(self.z, rhs.z),
        }
    }
}

impl<T: Copy> Point3<T> {
    /// Constructs a point with all three coordinates set to `value`.
    #[inline]
    pub const fn splat(value: T) -> Self {
        Self {
            x: value,
            y: value,
            z: value,
        }
    }

    /// Constructs a point from the first three elements of a slice.
    ///
    /// # Panics
    ///
    /// Panics if `slice` has fewer than three elements.
    #[inline]
    pub fn from_slice(slice: &[T]) -> Self {
        Self {
            x: slice[0],
            y: slice[1],
            z: slice[2],
        }
    }
}

impl<T> Index<usize> for Point3<T> {
    type Output = T;

    /// Returns the coordinate at `index`, where `0`, `1`, and `2` map to `x`,
    /// `y`, and `z`.
    ///
    /// # Panics
    ///
    /// Panics if `index` is greater than `2`.
    #[inline]
    fn index(&self, index: usize) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("index out of bounds: Point3 has 3 coordinates but the index is {index}"),
        }
    }
}

impl<T> IndexMut<usize> for Point3<T> {
    /// Returns the coordinate at `index`, where `0`, `1`, and `2` map to `x`,
    /// `y`, and `z`.
    ///
    /// # Panics
    ///
    /// Panics if `index` is greater than `2`.
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("index out of bounds: Point3 has 3 coordinates but the index is {index}"),
        }
    }
}

impl<T: Add<Output = T>> Add<Vector3<T>> for Point3<T> {
    type Output = Self;
    /// Translates the point by the displacement `rhs`.
    #[inline]
    fn add(self, rhs: Vector3<T>) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<T: AddAssign> AddAssign<Vector3<T>> for Point3<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Vector3<T>) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T: Sub<Output = T>> Sub<Vector3<T>> for Point3<T> {
    type Output = Self;
    /// Translates the point by the negated displacement `rhs`.
    #[inline]
    fn sub(self, rhs: Vector3<T>) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: SubAssign> SubAssign<Vector3<T>> for Point3<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Vector3<T>) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl<T: Sub<Output = T>> Sub for Point3<T> {
    type Output = Vector3<T>;
    /// Returns the displacement from `rhs` to `self`.
    #[inline]
    fn sub(self, rhs: Self) -> Vector3<T> {
        Vector3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<Q: Quantity> Point3<Q> {
    /// The origin, `(0, 0, 0)`.
    pub const ORIGIN: Self = Self::new(Q::ZERO, Q::ZERO, Q::ZERO);

    /// Returns the Euclidean distance between `self` and `rhs`.
    #[inline]
    pub fn distance(self, rhs: Self) -> Q {
        (self - rhs).norm()
    }

    /// Returns the centroid (arithmetic mean) of `points`, or the
    /// [`ORIGIN`][Self::ORIGIN] when `points` is empty.
    #[inline]
    pub fn centroid(points: &[Self]) -> Self {
        if points.is_empty() {
            return Self::ORIGIN;
        }
        let mut sum = Vector3::<Q>::ZERO;
        for point in points {
            sum += point.to_vector();
        }
        Self::from_vector(sum / Q::Value::from_f64(points.len() as f64))
    }

    /// Linearly interpolates from `self` toward `rhs` by the dimensionless
    /// factor `t`.
    ///
    /// `t == 0` yields `self`, `t == 1` yields `rhs`.
    #[inline]
    pub fn lerp(self, rhs: Self, t: Q::Value) -> Self {
        self + (rhs - self) * t
    }

    /// Returns the coordinate-wise absolute value.
    #[inline]
    pub fn abs(self) -> Self {
        self.map(Quantity::abs)
    }

    /// Returns the coordinate-wise minimum of `self` and `other`, ignoring
    /// NaN.
    ///
    /// Where one of the two is NaN, the other is taken.
    #[inline]
    pub fn min(self, other: Self) -> Self {
        self.zip_map(other, Quantity::min)
    }

    /// Returns the coordinate-wise maximum of `self` and `other`, ignoring
    /// NaN.
    ///
    /// Where one of the two is NaN, the other is taken.
    #[inline]
    pub fn max(self, other: Self) -> Self {
        self.zip_map(other, Quantity::max)
    }

    /// Restricts every coordinate to the interval `[lo, hi]`, confining the
    /// point to an axis-aligned bounding box.
    ///
    /// # Panics
    ///
    /// Panics if any coordinate of `min` exceeds the corresponding coordinate
    /// of `max`.
    #[inline]
    pub fn clamp(self, lo: Self, hi: Self) -> Self {
        Self::new(
            self.x.clamp(lo.x, hi.x),
            self.y.clamp(lo.y, hi.y),
            self.z.clamp(lo.z, hi.z),
        )
    }

    /// Returns the midpoint of the segment between `self` and `other`.
    #[inline]
    pub fn midpoint(self, other: Self) -> Self {
        self.zip_map(other, Quantity::midpoint)
    }

    /// Returns the smallest of the three coordinates.
    #[inline]
    pub fn min_element(self) -> Q {
        self.x.min(self.y).min(self.z)
    }

    /// Returns the largest of the three coordinates.
    #[inline]
    pub fn max_element(self) -> Q {
        self.x.max(self.y).max(self.z)
    }

    /// Returns the coordinate-wise sign, each `1.0`, `-1.0`, or NaN.
    #[inline]
    pub fn signum(self) -> Point3<Q::Value> {
        self.map(Quantity::signum)
    }

    /// Returns a point with the magnitudes of `self` and the coordinate-wise
    /// signs of `sign`.
    #[inline]
    pub fn copysign(self, sign: Self) -> Self {
        self.zip_map(sign, Quantity::copysign)
    }

    /// Returns the coordinate-wise floor (largest integer not exceeding each
    /// coordinate).
    #[inline]
    pub fn floor(self) -> Self {
        self.map(Quantity::floor)
    }

    /// Returns the coordinate-wise ceiling (smallest integer not less than each
    /// coordinate).
    #[inline]
    pub fn ceil(self) -> Self {
        self.map(Quantity::ceil)
    }

    /// Returns the coordinate-wise nearest integer, rounding halves away from
    /// zero.
    #[inline]
    pub fn round(self) -> Self {
        self.map(Quantity::round)
    }

    /// Returns the coordinate-wise nearest integer, rounding halves to even.
    #[inline]
    pub fn round_ties_even(self) -> Self {
        self.map(Quantity::round_ties_even)
    }

    /// Returns the coordinate-wise truncation toward zero.
    #[inline]
    pub fn trunc(self) -> Self {
        self.map(Quantity::trunc)
    }

    /// Returns the coordinate-wise fractional part, i.e. the position of
    /// `self` within the unit-grid cell it occupies.
    #[inline]
    pub fn fract(self) -> Self {
        self.map(Quantity::fract)
    }

    /// Returns the coordinate-wise Euclidean quotient against `rhs`.
    #[inline]
    pub fn div_euclid(self, rhs: Self) -> Point3<Q::Value> {
        self.zip_map(rhs, Quantity::div_euclid)
    }

    /// Returns the coordinate-wise least nonnegative remainder against `rhs`.
    #[inline]
    pub fn rem_euclid(self, rhs: Self) -> Self {
        self.zip_map(rhs, Quantity::rem_euclid)
    }

    /// Returns the coordinate-wise fused multiply-add `self * a + b`, each computed with a
    /// single rounding error.
    #[inline]
    pub fn mul_add(self, a: Q::Value, b: Self) -> Self {
        self.zip_map(b, |factor, addend| factor.mul_add(a, addend))
    }

    /// Returns the coordinate-wise hypotenuse of `self` and `other`, each computed without
    /// unnecessary overflow or underflow.
    #[inline]
    pub fn hypot(self, other: Self) -> Self {
        self.zip_map(other, Quantity::hypot)
    }

    /// Returns `true` if any coordinate is NaN.
    #[inline]
    pub fn is_nan(self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }

    /// Returns `true` if any coordinate is positive or negative infinity.
    #[inline]
    pub fn is_infinite(self) -> bool {
        self.x.is_infinite() || self.y.is_infinite() || self.z.is_infinite()
    }

    /// Returns `true` if every coordinate is finite.
    #[inline]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Returns `true` if every coordinate is neither zero, subnormal, infinite,
    /// nor NaN.
    #[inline]
    pub fn is_normal(self) -> bool {
        self.x.is_normal() && self.y.is_normal() && self.z.is_normal()
    }

    /// Returns `true` if any coordinate is subnormal.
    #[inline]
    pub fn is_subnormal(self) -> bool {
        self.x.is_subnormal() || self.y.is_subnormal() || self.z.is_subnormal()
    }
}

impl<V: Scalar> Point3<V> {
    /// Returns the squared Euclidean distance between `self` and `rhs`.
    ///
    /// Cheaper than [`distance`][Self::distance] and sufficient whenever only
    /// relative distances are compared.
    #[inline]
    pub fn distance_squared(self, rhs: Self) -> V {
        (self - rhs).norm_squared()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::units::length::{Angstrom, Length};

    const TOL: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL
    }

    fn points_close(a: Point3<f64>, b: Point3<f64>) -> bool {
        close(a.x, b.x) && close(a.y, b.y) && close(a.z, b.z)
    }

    fn length(value: f64) -> Length<f64, Angstrom> {
        Length::new(value)
    }

    fn point() -> Point3<f64> {
        Point3::new(1.0, 2.0, 3.0)
    }

    #[test]
    fn default_is_the_origin() {
        assert_eq!(Point3::<f64>::default(), Point3::ORIGIN);
    }

    #[test]
    fn the_origin_has_zero_coordinates() {
        assert_eq!(Point3::<f64>::ORIGIN, Point3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn centroid_of_no_points_is_the_origin() {
        assert_eq!(Point3::<f64>::centroid(&[]), Point3::ORIGIN);
    }

    #[test]
    fn centroid_of_one_point_is_that_point() {
        assert_eq!(Point3::centroid(&[point()]), point());
    }

    #[test]
    fn distance_from_a_point_to_itself_is_zero() {
        assert_eq!(point().distance(point()), 0.0);
    }

    #[test]
    fn distance_squared_from_a_point_to_itself_is_zero() {
        assert_eq!(point().distance_squared(point()), 0.0);
    }

    #[test]
    fn subtracting_a_point_from_itself_yields_the_zero_vector() {
        assert_eq!(point() - point(), Vector3::ZERO);
    }

    #[test]
    fn lerp_at_zero_yields_the_starting_point() {
        assert_eq!(point().lerp(Point3::new(4.0, 8.0, 12.0), 0.0), point());
    }

    #[test]
    fn lerp_at_one_yields_the_ending_point() {
        let target = Point3::new(4.0, 8.0, 12.0);
        assert_eq!(point().lerp(target, 1.0), target);
    }

    #[test]
    fn new_sets_the_three_coordinates() {
        let p = Point3::new(1.0, 2.0, 3.0);
        assert_eq!((p.x, p.y, p.z), (1.0, 2.0, 3.0));
    }

    #[test]
    fn from_array_takes_coordinates_in_order() {
        assert_eq!(Point3::from_array([1.0, 2.0, 3.0]), point());
    }

    #[test]
    fn to_array_yields_coordinates_in_order() {
        assert_eq!(point().to_array(), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn with_x_replaces_only_the_first_coordinate() {
        assert_eq!(point().with_x(9.0), Point3::new(9.0, 2.0, 3.0));
    }

    #[test]
    fn with_y_replaces_only_the_second_coordinate() {
        assert_eq!(point().with_y(9.0), Point3::new(1.0, 9.0, 3.0));
    }

    #[test]
    fn with_z_replaces_only_the_third_coordinate() {
        assert_eq!(point().with_z(9.0), Point3::new(1.0, 2.0, 9.0));
    }

    #[test]
    fn from_vector_reinterprets_a_displacement_as_a_point() {
        assert_eq!(Point3::from_vector(Vector3::new(1.0, 2.0, 3.0)), point());
    }

    #[test]
    fn to_vector_yields_the_position_relative_to_the_origin() {
        assert_eq!(point().to_vector(), Vector3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn map_applies_the_function_to_every_coordinate() {
        assert_eq!(point().map(|c| c as i32), Point3::new(1, 2, 3));
    }

    #[test]
    fn zip_map_combines_the_two_points_coordinate_wise() {
        let other = Point3::new(4.0, 5.0, 6.0);
        let combined = point().zip_map(other, |a, b| (a + b) as i32);
        assert_eq!(combined, Point3::new(5, 7, 9));
    }

    #[test]
    fn splat_repeats_one_value_across_all_coordinates() {
        assert_eq!(Point3::splat(5.0), Point3::new(5.0, 5.0, 5.0));
    }

    #[test]
    fn from_slice_takes_the_first_three_elements() {
        assert_eq!(Point3::from_slice(&[1.0, 2.0, 3.0, 4.0]), point());
    }

    #[test]
    fn indexing_yields_the_coordinate_at_that_position() {
        let p = point();
        assert_eq!((p[0], p[1], p[2]), (1.0, 2.0, 3.0));
    }

    #[test]
    fn index_mut_replaces_the_coordinate_at_that_position() {
        let mut p = point();
        p[1] = 9.0;
        assert_eq!(p, Point3::new(1.0, 9.0, 3.0));
    }

    #[test]
    fn adding_a_vector_translates_the_point() {
        assert_eq!(
            point() + Vector3::new(1.0, 1.0, 1.0),
            Point3::new(2.0, 3.0, 4.0)
        );
    }

    #[test]
    fn add_assign_translates_the_point_in_place() {
        let mut p = point();
        p += Vector3::new(1.0, 1.0, 1.0);
        assert_eq!(p, Point3::new(2.0, 3.0, 4.0));
    }

    #[test]
    fn subtracting_a_vector_translates_the_point_backward() {
        assert_eq!(
            point() - Vector3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 2.0)
        );
    }

    #[test]
    fn sub_assign_translates_the_point_in_place() {
        let mut p = point();
        p -= Vector3::new(1.0, 1.0, 1.0);
        assert_eq!(p, Point3::new(0.0, 1.0, 2.0));
    }

    #[test]
    fn subtracting_two_points_yields_the_displacement_between_them() {
        let from = Point3::new(1.0, 1.0, 1.0);
        assert_eq!(point() - from, Vector3::new(0.0, 1.0, 2.0));
    }

    #[test]
    fn distance_is_the_euclidean_distance_between_the_points() {
        assert_eq!(Point3::ORIGIN.distance(Point3::new(3.0, 4.0, 0.0)), 5.0);
    }

    #[test]
    fn distance_squared_is_the_squared_euclidean_distance() {
        assert_eq!(
            Point3::ORIGIN.distance_squared(Point3::new(3.0, 4.0, 0.0)),
            25.0
        );
    }

    #[test]
    fn centroid_is_the_arithmetic_mean_of_the_points() {
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 6.0, 9.0),
            Point3::new(6.0, 12.0, 18.0),
        ];
        assert_eq!(Point3::centroid(&points), Point3::new(3.0, 6.0, 9.0));
    }

    #[test]
    fn lerp_interpolates_between_the_points() {
        let target = Point3::new(4.0, 8.0, 12.0);
        assert_eq!(
            Point3::ORIGIN.lerp(target, 0.25),
            Point3::new(1.0, 2.0, 3.0)
        );
    }

    #[test]
    fn abs_takes_the_magnitude_of_every_coordinate() {
        let p = Point3::new(-1.0, 2.0, -3.0);
        assert_eq!(p.abs(), point());
    }

    #[test]
    fn min_takes_the_smaller_of_each_pair_of_coordinates() {
        let other = Point3::new(4.0, 2.0, 3.0);
        assert_eq!(Point3::new(1.0, 5.0, 3.0).min(other), point());
    }

    #[test]
    fn max_takes_the_larger_of_each_pair_of_coordinates() {
        let other = Point3::new(4.0, 2.0, 3.0);
        assert_eq!(
            Point3::new(1.0, 5.0, 3.0).max(other),
            Point3::new(4.0, 5.0, 3.0)
        );
    }

    #[test]
    fn clamp_raises_a_coordinate_below_the_interval_to_the_lower_bound() {
        let p = Point3::new(-1.0, 2.0, 3.0);
        assert_eq!(
            p.clamp(Point3::ORIGIN, Point3::splat(9.0)),
            Point3::new(0.0, 2.0, 3.0)
        );
    }

    #[test]
    fn clamp_leaves_a_coordinate_inside_the_interval_unchanged() {
        assert_eq!(point().clamp(Point3::ORIGIN, Point3::splat(9.0)), point());
    }

    #[test]
    fn clamp_lowers_a_coordinate_above_the_interval_to_the_upper_bound() {
        let p = Point3::new(1.0, 2.0, 99.0);
        assert_eq!(
            p.clamp(Point3::ORIGIN, Point3::splat(9.0)),
            Point3::new(1.0, 2.0, 9.0)
        );
    }

    #[test]
    fn midpoint_is_halfway_between_the_points() {
        let target = Point3::new(2.0, 4.0, 6.0);
        assert_eq!(Point3::ORIGIN.midpoint(target), Point3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn min_element_is_the_smallest_coordinate() {
        assert_eq!(Point3::new(3.0, -1.0, 2.0).min_element(), -1.0);
    }

    #[test]
    fn max_element_is_the_largest_coordinate() {
        assert_eq!(Point3::new(3.0, -1.0, 2.0).max_element(), 3.0);
    }

    #[test]
    fn signum_takes_the_sign_of_every_coordinate() {
        let p = Point3::new(-2.0, 0.5, -3.0);
        assert_eq!(p.signum(), Point3::new(-1.0, 1.0, -1.0));
    }

    #[test]
    fn copysign_keeps_the_magnitudes_and_takes_the_signs_of_its_argument() {
        let p = Point3::new(3.0, -4.0, 5.0);
        let signs = Point3::new(-1.0, 1.0, -1.0);
        assert_eq!(p.copysign(signs), Point3::new(-3.0, 4.0, -5.0));
    }

    #[test]
    fn floor_rounds_every_coordinate_toward_negative_infinity() {
        let p = Point3::new(-2.5, 2.5, 3.0);
        assert_eq!(p.floor(), Point3::new(-3.0, 2.0, 3.0));
    }

    #[test]
    fn ceil_rounds_every_coordinate_toward_positive_infinity() {
        let p = Point3::new(-2.5, 2.5, 3.0);
        assert_eq!(p.ceil(), Point3::new(-2.0, 3.0, 3.0));
    }

    #[test]
    fn round_sends_a_half_away_from_zero() {
        let p = Point3::new(2.5, 3.5, -2.5);
        assert_eq!(p.round(), Point3::new(3.0, 4.0, -3.0));
    }

    #[test]
    fn round_ties_even_sends_a_half_to_the_even_integer() {
        let p = Point3::new(2.5, 3.5, -2.5);
        assert_eq!(p.round_ties_even(), Point3::new(2.0, 4.0, -2.0));
    }

    #[test]
    fn trunc_drops_the_fractional_part_of_every_coordinate() {
        let p = Point3::new(-2.75, 2.75, 3.0);
        assert_eq!(p.trunc(), Point3::new(-2.0, 2.0, 3.0));
    }

    #[test]
    fn fract_keeps_the_fractional_part_of_every_coordinate() {
        let p = Point3::new(-2.75, 2.75, 3.0);
        assert_eq!(p.fract(), Point3::new(-0.75, 0.75, 0.0));
    }

    #[test]
    fn div_euclid_is_the_euclidean_quotient_of_every_coordinate() {
        let p = Point3::new(7.0, -7.0, 0.0);
        assert_eq!(
            p.div_euclid(Point3::splat(3.0)),
            Point3::new(2.0, -3.0, 0.0)
        );
    }

    #[test]
    fn rem_euclid_is_nonnegative_for_a_negative_coordinate() {
        let p = Point3::new(-1.0, 7.0, 5.0);
        assert_eq!(p.rem_euclid(Point3::splat(3.0)), Point3::new(2.0, 1.0, 2.0));
    }

    #[test]
    fn mul_add_scales_then_offsets_every_coordinate() {
        let offset = Point3::new(10.0, 20.0, 30.0);
        assert_eq!(point().mul_add(2.0, offset), Point3::new(12.0, 24.0, 36.0));
    }

    #[test]
    fn hypot_combines_the_coordinates_pairwise() {
        let legs = Point3::new(4.0, 0.0, 5.0);
        let combined = Point3::new(3.0, 4.0, 12.0).hypot(legs);
        assert!(points_close(combined, Point3::new(5.0, 4.0, 13.0)));
    }

    #[test]
    fn is_nan_holds_when_a_coordinate_is_not_a_number() {
        assert!(Point3::new(1.0, f64::NAN, 3.0).is_nan());
    }

    #[test]
    fn is_infinite_holds_when_a_coordinate_is_infinite() {
        assert!(Point3::new(1.0, f64::INFINITY, 3.0).is_infinite());
    }

    #[test]
    fn is_finite_holds_when_every_coordinate_is_finite() {
        assert!(point().is_finite());
    }

    #[test]
    fn is_normal_holds_when_every_coordinate_is_normal() {
        assert!(point().is_normal());
    }

    #[test]
    fn is_subnormal_holds_when_a_coordinate_is_subnormal() {
        assert!(Point3::new(1.0, f64::MIN_POSITIVE / 2.0, 3.0).is_subnormal());
    }

    #[test]
    fn signum_of_a_length_point_is_dimensionless() {
        let p = Point3::new(length(-2.0), length(0.5), length(-3.0));
        assert_eq!(p.signum(), Point3::new(-1.0, 1.0, -1.0));
    }

    #[test]
    fn div_euclid_of_length_points_is_dimensionless() {
        let p = Point3::new(length(7.0), length(-7.0), length(0.0));
        let divisor = Point3::splat(length(3.0));
        assert_eq!(p.div_euclid(divisor), Point3::new(2.0, -3.0, 0.0));
    }

    #[test]
    fn mul_add_of_length_points_takes_a_dimensionless_factor() {
        let p = Point3::new(length(1.0), length(2.0), length(3.0));
        let offset = Point3::new(length(10.0), length(20.0), length(30.0));
        let expected = Point3::new(length(12.0), length(24.0), length(36.0));
        assert_eq!(p.mul_add(2.0, offset), expected);
    }

    #[test]
    fn hypot_of_length_points_stays_a_length() {
        let p = Point3::new(length(3.0), length(4.0), length(12.0));
        let legs = Point3::new(length(4.0), length(0.0), length(5.0));
        let combined = p.hypot(legs).map(Quantity::value);
        assert!(points_close(combined, Point3::new(5.0, 4.0, 13.0)));
    }

    #[test]
    fn lerp_of_length_points_takes_a_dimensionless_factor() {
        let target = Point3::new(length(2.0), length(4.0), length(6.0));
        let expected = Point3::new(length(1.0), length(2.0), length(3.0));
        assert_eq!(Point3::ORIGIN.lerp(target, 0.5), expected);
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn from_slice_shorter_than_three_panics() {
        let _ = Point3::from_slice(&[1.0_f64, 2.0]);
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn indexing_out_of_bounds_panics() {
        let _ = point()[3];
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn index_mut_out_of_bounds_panics() {
        point()[3] = 0.0;
    }

    #[test]
    #[should_panic(expected = "min > max")]
    fn clamp_with_an_inverted_interval_panics() {
        let _ = point().clamp(Point3::splat(9.0), Point3::ORIGIN);
    }

    #[test]
    fn min_ignores_a_not_a_number_coordinate() {
        let nan = Point3::new(f64::NAN, f64::NAN, f64::NAN);
        assert_eq!(nan.min(point()), point());
    }

    #[test]
    fn max_ignores_a_not_a_number_coordinate() {
        let nan = Point3::new(f64::NAN, f64::NAN, f64::NAN);
        assert_eq!(nan.max(point()), point());
    }

    #[test]
    fn is_nan_does_not_hold_when_no_coordinate_is_not_a_number() {
        assert!(!point().is_nan());
    }

    #[test]
    fn is_infinite_does_not_hold_when_every_coordinate_is_finite() {
        assert!(!point().is_infinite());
    }

    #[test]
    fn is_finite_does_not_hold_when_a_coordinate_is_infinite() {
        assert!(!Point3::new(1.0, f64::INFINITY, 3.0).is_finite());
    }

    #[test]
    fn is_normal_does_not_hold_when_a_coordinate_is_zero() {
        assert!(!Point3::new(1.0, 0.0, 3.0).is_normal());
    }

    #[test]
    fn is_subnormal_does_not_hold_when_every_coordinate_is_normal() {
        assert!(!point().is_subnormal());
    }

    #[test]
    fn from_slice_of_exactly_three_elements_yields_the_point() {
        assert_eq!(Point3::from_slice(&[1.0, 2.0, 3.0]), point());
    }

    #[test]
    fn floor_leaves_an_integer_coordinate_unchanged() {
        assert_eq!(point().floor(), point());
    }

    #[test]
    fn ceil_leaves_an_integer_coordinate_unchanged() {
        assert_eq!(point().ceil(), point());
    }

    #[test]
    fn fract_of_an_integer_coordinate_is_zero() {
        assert_eq!(point().fract(), Point3::ORIGIN);
    }

    #[test]
    fn signum_of_a_negative_zero_coordinate_is_minus_one() {
        assert_eq!(Point3::new(-0.0, 1.0, 1.0).signum().x, -1.0);
    }

    #[test]
    fn lerp_beyond_one_extrapolates_past_the_ending_point() {
        let target = Point3::new(4.0, 8.0, 12.0);
        assert_eq!(
            Point3::ORIGIN.lerp(target, 2.0),
            Point3::new(8.0, 16.0, 24.0)
        );
    }

    #[test]
    fn array_roundtrip_preserves_the_point() {
        assert_eq!(Point3::from_array(point().to_array()), point());
    }

    #[test]
    fn vector_roundtrip_preserves_the_point() {
        assert_eq!(Point3::from_vector(point().to_vector()), point());
    }

    #[test]
    fn translating_by_a_vector_and_back_returns_the_point() {
        let displacement = Vector3::new(4.0, -5.0, 6.0);
        assert_eq!(point() + displacement - displacement, point());
    }

    #[test]
    fn the_displacement_between_points_translates_one_to_the_other() {
        let other = Point3::new(4.0, -5.0, 6.0);
        assert_eq!(point() + (other - point()), other);
    }

    #[test]
    fn distance_is_symmetric() {
        let other = Point3::new(4.0, -5.0, 6.0);
        assert_eq!(point().distance(other), other.distance(point()));
    }

    #[test]
    fn distance_squared_is_the_square_of_the_distance() {
        let other = Point3::new(4.0, 6.0, 3.0);
        let distance = point().distance(other);
        assert_eq!(point().distance_squared(other), distance * distance);
    }

    #[test]
    fn midpoint_is_equidistant_from_both_points() {
        let other = Point3::new(4.0, 8.0, 12.0);
        let middle = point().midpoint(other);
        assert_eq!(point().distance(middle), other.distance(middle));
    }

    #[test]
    fn midpoint_agrees_with_interpolating_halfway() {
        let other = Point3::new(4.0, 8.0, 12.0);
        assert_eq!(point().midpoint(other), point().lerp(other, 0.5));
    }

    #[test]
    fn centroid_is_independent_of_the_order_of_the_points() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(2.0, 4.0, 6.0);
        let c = Point3::new(4.0, 8.0, 12.0);
        assert_eq!(Point3::centroid(&[a, b, c]), Point3::centroid(&[c, a, b]));
    }

    #[test]
    fn truncating_and_the_fractional_part_reconstruct_the_point() {
        let p = Point3::new(-2.75, 2.75, 3.0);
        assert_eq!(p.trunc() + p.fract().to_vector(), p);
    }

    #[test]
    fn min_of_a_point_with_itself_is_that_point() {
        assert_eq!(point().min(point()), point());
    }

    #[test]
    fn max_of_a_point_with_itself_is_that_point() {
        assert_eq!(point().max(point()), point());
    }

    #[test]
    fn clamp_agrees_with_taking_the_maximum_then_the_minimum() {
        let (lo, hi) = (Point3::splat(1.5), Point3::splat(2.5));
        assert_eq!(point().clamp(lo, hi), point().max(lo).min(hi));
    }

    #[test]
    fn points_are_equal_exactly_when_all_coordinates_match() {
        assert_eq!(point(), Point3::new(1.0, 2.0, 3.0));
        assert_ne!(point(), Point3::new(1.0, 2.0, 4.0));
    }
}

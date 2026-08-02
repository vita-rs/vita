use core::iter::Sum;
use core::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

use crate::{Quantity, Scalar};

/// A vector of three components `x`, `y`, and `z`.
///
/// A vector is a displacement: it adds to another vector, scales by a number,
/// and negates. Its length and direction are intrinsic, its position is not.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vector3<T> {
    /// The first component.
    pub x: T,
    /// The second component.
    pub y: T,
    /// The third component.
    pub z: T,
}

impl<T> Vector3<T> {
    /// Constructs a vector from its three components.
    #[inline]
    pub const fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    /// Constructs a vector from an array `[x, y, z]`.
    #[inline]
    pub fn from_array(array: [T; 3]) -> Self {
        let [x, y, z] = array;
        Self { x, y, z }
    }

    /// Returns the components as an array `[x, y, z]`.
    #[inline]
    pub fn to_array(self) -> [T; 3] {
        [self.x, self.y, self.z]
    }

    /// Returns a copy of `self` with the `x` component replaced.
    #[inline]
    pub fn with_x(self, x: T) -> Self {
        Self {
            x,
            y: self.y,
            z: self.z,
        }
    }

    /// Returns a copy of `self` with the `y` component replaced.
    #[inline]
    pub fn with_y(self, y: T) -> Self {
        Self {
            x: self.x,
            y,
            z: self.z,
        }
    }

    /// Returns a copy of `self` with the `z` component replaced.
    #[inline]
    pub fn with_z(self, z: T) -> Self {
        Self {
            x: self.x,
            y: self.y,
            z,
        }
    }

    /// Applies `f` to every component, returning the resulting vector.
    #[inline]
    pub fn map<U, F: FnMut(T) -> U>(self, mut f: F) -> Vector3<U> {
        Vector3 {
            x: f(self.x),
            y: f(self.y),
            z: f(self.z),
        }
    }

    /// Combines `self` and `rhs` component-wise through `f`.
    #[inline]
    pub fn zip_map<U, R, F: FnMut(T, U) -> R>(self, rhs: Vector3<U>, mut f: F) -> Vector3<R> {
        Vector3 {
            x: f(self.x, rhs.x),
            y: f(self.y, rhs.y),
            z: f(self.z, rhs.z),
        }
    }
}

impl<T: Copy> Vector3<T> {
    /// Constructs a vector with all three components set to `value`.
    #[inline]
    pub const fn splat(value: T) -> Self {
        Self {
            x: value,
            y: value,
            z: value,
        }
    }

    /// Constructs a vector from the first three elements of a slice.
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

impl<T> Index<usize> for Vector3<T> {
    type Output = T;

    /// Returns the component at `index`, where `0`, `1`, and `2` map to `x`,
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
            _ => panic!("index out of bounds: Vector3 has 3 components but the index is {index}"),
        }
    }
}

impl<T> IndexMut<usize> for Vector3<T> {
    /// Returns the component at `index`, where `0`, `1`, and `2` map to `x`,
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
            _ => panic!("index out of bounds: Vector3 has 3 components but the index is {index}"),
        }
    }
}

impl<T: Neg<Output = T>> Neg for Vector3<T> {
    type Output = Self;
    /// Returns the component-wise negation of `self`.
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl<T: Add<Output = T>> Add for Vector3<T> {
    type Output = Self;
    /// Returns the component-wise sum of `self` and `rhs`.
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<T: AddAssign> AddAssign for Vector3<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T: Sub<Output = T>> Sub for Vector3<T> {
    type Output = Self;
    /// Returns the component-wise difference of `self` and `rhs`.
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: SubAssign> SubAssign for Vector3<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl<T, S: Quantity> Mul<S> for Vector3<T>
where
    T: Mul<S>,
{
    type Output = Vector3<T::Output>;
    /// Scales every component by `rhs`, whose dimension multiplies into the
    /// component's.
    #[inline]
    fn mul(self, rhs: S) -> Self::Output {
        Vector3::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<T, S: Quantity> MulAssign<S> for Vector3<T>
where
    T: MulAssign<S>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: S) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl<T> Mul<Vector3<T>> for f32
where
    f32: Mul<T>,
{
    type Output = Vector3<<f32 as Mul<T>>::Output>;
    /// Scales every component of `rhs` by `self`.
    #[inline]
    fn mul(self, rhs: Vector3<T>) -> Self::Output {
        rhs.map(|component| self * component)
    }
}

impl<T> Mul<Vector3<T>> for f64
where
    f64: Mul<T>,
{
    type Output = Vector3<<f64 as Mul<T>>::Output>;
    /// Scales every component of `rhs` by `self`.
    #[inline]
    fn mul(self, rhs: Vector3<T>) -> Self::Output {
        rhs.map(|component| self * component)
    }
}

impl<T, S: Quantity> Div<S> for Vector3<T>
where
    T: Div<S>,
{
    type Output = Vector3<T::Output>;
    /// Divides every component by `rhs`, whose dimension divides out of the
    /// component's.
    #[inline]
    fn div(self, rhs: S) -> Self::Output {
        Vector3::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl<T, S: Quantity> DivAssign<S> for Vector3<T>
where
    T: DivAssign<S>,
{
    #[inline]
    fn div_assign(&mut self, rhs: S) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl<T: Add<Output = T> + Default> Sum for Vector3<T> {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |acc, v| acc + v)
    }
}

impl<'a, T: Add<Output = T> + Default + Copy> Sum<&'a Vector3<T>> for Vector3<T> {
    #[inline]
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.copied().fold(Self::default(), |acc, v| acc + v)
    }
}

impl<Q: Quantity> Vector3<Q> {
    /// The zero vector, `(0, 0, 0)`.
    pub const ZERO: Self = Self::new(Q::ZERO, Q::ZERO, Q::ZERO);

    /// Returns the Euclidean norm (length) of the vector.
    ///
    /// The norm is representable whenever the components are, so it is computed by the
    /// direct route only while the sum of squares stays normal, and by a pairwise
    /// [`hypot`](Quantity::hypot) otherwise — squaring alone would halve the usable range.
    #[inline]
    pub fn norm(self) -> Q {
        let values = self.values();
        let square = values.norm_squared();
        if square.is_normal() {
            return Q::from_value(square.sqrt());
        }
        Q::from_value(values.x.hypot(values.y).hypot(values.z))
    }

    /// Returns the unit vector along `self`, whose components carry no
    /// dimension.
    ///
    /// Yields a non-finite vector when the norm is zero or non-finite; use
    /// [`try_normalize`][Self::try_normalize] or
    /// [`normalize_or_zero`][Self::normalize_or_zero] to handle those cases.
    #[inline]
    pub fn normalize(self) -> Vector3<Q::Value> {
        self.values() / self.norm().value()
    }

    /// Returns the unit vector along `self`, or `None` if the result would not
    /// be finite (e.g. for the zero vector).
    #[inline]
    pub fn try_normalize(self) -> Option<Vector3<Q::Value>> {
        let norm = self.norm().value();
        if norm.is_finite() && norm > Q::Value::ZERO {
            Some(self.values() / norm)
        } else {
            None
        }
    }

    /// Returns the unit vector along `self`, or the zero vector if the result
    /// would not be finite.
    #[inline]
    pub fn normalize_or_zero(self) -> Vector3<Q::Value> {
        self.try_normalize().unwrap_or(Vector3::ZERO)
    }

    /// Returns the unsigned angle between `self` and `rhs`, in radians within
    /// `[0, π]`.
    ///
    /// The computation is numerically stable across the whole range, including
    /// near-parallel and near-antiparallel inputs. An angle depends only on the two
    /// directions, so once the products carrying it leave the representable range — which
    /// the angle itself never does — the directions are taken first instead.
    #[inline]
    pub fn angle_between(self, rhs: Self) -> Q::Value {
        let (left, right) = (self.values(), rhs.values());
        let (perpendicular, parallel) = (left.cross(right).norm(), left.dot(right));
        if perpendicular.max(parallel.abs()).is_normal() {
            return perpendicular.atan2(parallel);
        }
        let (left, right) = (self.normalize(), rhs.normalize());
        left.cross(right).norm().atan2(left.dot(right))
    }

    /// Returns the vector projection of `self` onto `onto`.
    ///
    /// The projection depends on `onto` only through its direction, so once the dot
    /// products carrying it leave the representable range — which the projection itself
    /// never does — that direction is taken first instead.
    #[inline]
    pub fn project_onto(self, onto: Self) -> Self {
        let values = onto.values();
        let (square, along) = (values.norm_squared(), self.values().dot(values));
        if square.is_normal() && along.is_finite() {
            return onto * (along / square);
        }
        let direction = onto.normalize();
        let length = self.x * direction.x + self.y * direction.y + self.z * direction.z;
        direction.map(|component| length * component)
    }

    /// Returns the component of `self` orthogonal to `from`.
    #[inline]
    pub fn reject_from(self, from: Self) -> Self {
        self - self.project_onto(from)
    }

    /// Reflects `self` across the plane through the origin with dimensionless
    /// unit normal `normal`.
    #[inline]
    pub fn reflect(self, normal: Vector3<Q::Value>) -> Self {
        let values = self.values();
        let twice = Q::Value::ONE + Q::Value::ONE;
        self - (normal * (values.dot(normal) * twice)).map(Q::from_value)
    }

    /// Linearly interpolates from `self` toward `rhs` by the dimensionless
    /// factor `t`.
    ///
    /// `t == 0` yields `self`, `t == 1` yields `rhs`.
    #[inline]
    pub fn lerp(self, rhs: Self, t: Q::Value) -> Self {
        let complement = Q::Value::ONE - t;
        self.zip_map(rhs, |start, end| start * complement + end * t)
    }

    /// Returns the component-wise absolute value.
    #[inline]
    pub fn abs(self) -> Self {
        self.map(Quantity::abs)
    }

    /// Returns the component-wise minimum of `self` and `other`, ignoring
    /// NaN.
    ///
    /// Where one of the two is NaN, the other is taken.
    #[inline]
    pub fn min(self, other: Self) -> Self {
        self.zip_map(other, Quantity::min)
    }

    /// Returns the component-wise maximum of `self` and `other`, ignoring
    /// NaN.
    ///
    /// Where one of the two is NaN, the other is taken.
    #[inline]
    pub fn max(self, other: Self) -> Self {
        self.zip_map(other, Quantity::max)
    }

    /// Restricts every component to the interval `[lo, hi]`.
    ///
    /// # Panics
    ///
    /// Panics if any component of `lo` exceeds the corresponding component of
    /// `hi`, or if either is NaN.
    #[inline]
    pub fn clamp(self, lo: Self, hi: Self) -> Self {
        Self::new(
            self.x.clamp(lo.x, hi.x),
            self.y.clamp(lo.y, hi.y),
            self.z.clamp(lo.z, hi.z),
        )
    }

    /// Returns the component-wise midpoint of `self` and `other`.
    #[inline]
    pub fn midpoint(self, other: Self) -> Self {
        self.zip_map(other, Quantity::midpoint)
    }

    /// Returns the smallest of the three components.
    #[inline]
    pub fn min_element(self) -> Q {
        self.x.min(self.y).min(self.z)
    }

    /// Returns the largest of the three components.
    #[inline]
    pub fn max_element(self) -> Q {
        self.x.max(self.y).max(self.z)
    }

    /// Returns the sum of the three components, `x + y + z`.
    #[inline]
    pub fn element_sum(self) -> Q {
        self.x + self.y + self.z
    }

    /// Returns the component-wise sign, each `1.0`, `-1.0`, or NaN.
    #[inline]
    pub fn signum(self) -> Vector3<Q::Value> {
        self.map(Quantity::signum)
    }

    /// Returns a vector with the magnitudes of `self` and the component-wise
    /// signs of `sign`.
    #[inline]
    pub fn copysign(self, sign: Self) -> Self {
        self.zip_map(sign, Quantity::copysign)
    }

    /// Returns the component-wise floor.
    #[inline]
    pub fn floor(self) -> Self {
        self.map(Quantity::floor)
    }

    /// Returns the component-wise ceiling.
    #[inline]
    pub fn ceil(self) -> Self {
        self.map(Quantity::ceil)
    }

    /// Returns the component-wise nearest integer, rounding halves away from
    /// zero.
    #[inline]
    pub fn round(self) -> Self {
        self.map(Quantity::round)
    }

    /// Returns the component-wise nearest integer, rounding halves to even.
    #[inline]
    pub fn round_ties_even(self) -> Self {
        self.map(Quantity::round_ties_even)
    }

    /// Returns the component-wise truncation toward zero.
    #[inline]
    pub fn trunc(self) -> Self {
        self.map(Quantity::trunc)
    }

    /// Returns the component-wise fractional part.
    #[inline]
    pub fn fract(self) -> Self {
        self.map(Quantity::fract)
    }

    /// Returns the component-wise Euclidean quotient against `rhs`.
    #[inline]
    pub fn div_euclid(self, rhs: Self) -> Vector3<Q::Value> {
        self.zip_map(rhs, Quantity::div_euclid)
    }

    /// Returns the component-wise least nonnegative remainder against `rhs`.
    #[inline]
    pub fn rem_euclid(self, rhs: Self) -> Self {
        self.zip_map(rhs, Quantity::rem_euclid)
    }

    /// Returns the component-wise fused multiply-add `self * a + b`, each computed with a
    /// single rounding error.
    #[inline]
    pub fn mul_add(self, a: Q::Value, b: Self) -> Self {
        self.zip_map(b, |factor, addend| factor.mul_add(a, addend))
    }

    /// Returns the component-wise hypotenuse of `self` and `other`, each computed without
    /// unnecessary overflow or underflow.
    #[inline]
    pub fn hypot(self, other: Self) -> Self {
        self.zip_map(other, Quantity::hypot)
    }

    /// Returns `true` if any component is NaN.
    #[inline]
    pub fn is_nan(self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }

    /// Returns `true` if any component is positive or negative infinity.
    #[inline]
    pub fn is_infinite(self) -> bool {
        self.x.is_infinite() || self.y.is_infinite() || self.z.is_infinite()
    }

    /// Returns `true` if every component is finite.
    #[inline]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Returns `true` if every component is neither zero, subnormal, infinite,
    /// nor NaN.
    #[inline]
    pub fn is_normal(self) -> bool {
        self.x.is_normal() && self.y.is_normal() && self.z.is_normal()
    }

    /// Returns `true` if any component is subnormal.
    #[inline]
    pub fn is_subnormal(self) -> bool {
        self.x.is_subnormal() || self.y.is_subnormal() || self.z.is_subnormal()
    }

    /// Returns the component-wise underlying numbers, read in `Q`'s own unit.
    #[inline]
    fn values(self) -> Vector3<Q::Value> {
        self.map(Quantity::value)
    }
}

impl<V: Scalar> Vector3<V> {
    /// The vector with every component set to one, `(1, 1, 1)`.
    pub const ONE: Self = Self::new(V::ONE, V::ONE, V::ONE);

    /// The unit vector along the `x` axis, `(1, 0, 0)`.
    pub const X: Self = Self::new(V::ONE, V::ZERO, V::ZERO);

    /// The unit vector along the `y` axis, `(0, 1, 0)`.
    pub const Y: Self = Self::new(V::ZERO, V::ONE, V::ZERO);

    /// The unit vector along the `z` axis, `(0, 0, 1)`.
    pub const Z: Self = Self::new(V::ZERO, V::ZERO, V::ONE);

    /// Returns the dot product `self · rhs`.
    #[inline]
    pub fn dot(self, rhs: Self) -> V {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    /// Returns the cross product `self × rhs`.
    #[inline]
    pub fn cross(self, rhs: Self) -> Self {
        Self::new(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x,
        )
    }

    /// Returns the squared Euclidean norm `self · self`.
    ///
    /// Cheaper than [`norm`][Self::norm] and sufficient whenever only relative
    /// magnitudes are compared.
    #[inline]
    pub fn norm_squared(self) -> V {
        self.dot(self)
    }

    /// Returns the component-wise reciprocal `(1/x, 1/y, 1/z)`.
    #[inline]
    pub fn recip(self) -> Self {
        self.map(Scalar::recip)
    }

    /// Returns the product of the three components, `x * y * z`.
    #[inline]
    pub fn element_product(self) -> V {
        self.x * self.y * self.z
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use core::f64::consts::{FRAC_PI_2, PI};

    use crate::units::length::{Angstrom, Length};

    const TOL: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL
    }

    fn vectors_close(a: Vector3<f64>, b: Vector3<f64>) -> bool {
        close(a.x, b.x) && close(a.y, b.y) && close(a.z, b.z)
    }

    fn length(value: f64) -> Length<f64, Angstrom> {
        Length::new(value)
    }

    fn vector() -> Vector3<f64> {
        Vector3::new(1.0, 2.0, 3.0)
    }

    #[test]
    fn default_is_the_zero_vector() {
        assert_eq!(Vector3::<f64>::default(), Vector3::ZERO);
    }

    #[test]
    fn the_zero_vector_has_zero_components() {
        assert_eq!(Vector3::<f64>::ZERO, Vector3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn norm_of_the_zero_vector_is_zero() {
        assert_eq!(Vector3::<f64>::ZERO.norm(), 0.0);
    }

    #[test]
    fn try_normalize_of_the_zero_vector_is_none() {
        assert!(Vector3::<f64>::ZERO.try_normalize().is_none());
    }

    #[test]
    fn normalize_or_zero_of_the_zero_vector_is_the_zero_vector() {
        assert_eq!(Vector3::<f64>::ZERO.normalize_or_zero(), Vector3::ZERO);
    }

    #[test]
    fn element_sum_of_the_zero_vector_is_zero() {
        assert_eq!(Vector3::<f64>::ZERO.element_sum(), 0.0);
    }

    #[test]
    fn summing_no_vectors_yields_the_zero_vector() {
        let none: [Vector3<f64>; 0] = [];
        assert_eq!(none.into_iter().sum::<Vector3<f64>>(), Vector3::ZERO);
    }

    #[test]
    fn angle_between_a_vector_and_itself_is_zero() {
        assert!(close(vector().angle_between(vector()), 0.0));
    }

    #[test]
    fn lerp_at_zero_yields_the_starting_vector() {
        assert_eq!(vector().lerp(Vector3::new(4.0, 8.0, 12.0), 0.0), vector());
    }

    #[test]
    fn lerp_at_one_yields_the_ending_vector() {
        let target = Vector3::new(4.0, 8.0, 12.0);
        assert_eq!(vector().lerp(target, 1.0), target);
    }

    #[test]
    fn new_sets_the_three_components() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!((v.x, v.y, v.z), (1.0, 2.0, 3.0));
    }

    #[test]
    fn from_array_takes_components_in_order() {
        assert_eq!(Vector3::from_array([1.0, 2.0, 3.0]), vector());
    }

    #[test]
    fn to_array_yields_components_in_order() {
        assert_eq!(vector().to_array(), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn with_x_replaces_only_the_first_component() {
        assert_eq!(vector().with_x(9.0), Vector3::new(9.0, 2.0, 3.0));
    }

    #[test]
    fn with_y_replaces_only_the_second_component() {
        assert_eq!(vector().with_y(9.0), Vector3::new(1.0, 9.0, 3.0));
    }

    #[test]
    fn with_z_replaces_only_the_third_component() {
        assert_eq!(vector().with_z(9.0), Vector3::new(1.0, 2.0, 9.0));
    }

    #[test]
    fn map_applies_the_function_to_every_component() {
        assert_eq!(vector().map(|c| c as i32), Vector3::new(1, 2, 3));
    }

    #[test]
    fn zip_map_combines_the_two_vectors_component_wise() {
        let other = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(
            vector().zip_map(other, |a, b| (a + b) as i32),
            Vector3::new(5, 7, 9)
        );
    }

    #[test]
    fn splat_repeats_one_value_across_all_components() {
        assert_eq!(Vector3::splat(5.0), Vector3::new(5.0, 5.0, 5.0));
    }

    #[test]
    fn from_slice_takes_the_first_three_elements() {
        assert_eq!(Vector3::from_slice(&[1.0, 2.0, 3.0, 4.0]), vector());
    }

    #[test]
    fn indexing_yields_the_component_at_that_position() {
        let v = vector();
        assert_eq!((v[0], v[1], v[2]), (1.0, 2.0, 3.0));
    }

    #[test]
    fn index_mut_replaces_the_component_at_that_position() {
        let mut v = vector();
        v[1] = 9.0;
        assert_eq!(v, Vector3::new(1.0, 9.0, 3.0));
    }

    #[test]
    fn negation_flips_the_sign_of_every_component() {
        assert_eq!(-vector(), Vector3::new(-1.0, -2.0, -3.0));
    }

    #[test]
    fn addition_sums_the_components() {
        assert_eq!(vector() + Vector3::splat(1.0), Vector3::new(2.0, 3.0, 4.0));
    }

    #[test]
    fn add_assign_sums_the_components_in_place() {
        let mut v = vector();
        v += Vector3::splat(1.0);
        assert_eq!(v, Vector3::new(2.0, 3.0, 4.0));
    }

    #[test]
    fn subtraction_differences_the_components() {
        assert_eq!(vector() - Vector3::splat(1.0), Vector3::new(0.0, 1.0, 2.0));
    }

    #[test]
    fn sub_assign_differences_the_components_in_place() {
        let mut v = vector();
        v -= Vector3::splat(1.0);
        assert_eq!(v, Vector3::new(0.0, 1.0, 2.0));
    }

    #[test]
    fn multiplying_by_a_number_scales_every_component() {
        assert_eq!(vector() * 2.0, Vector3::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn mul_assign_scales_every_component_in_place() {
        let mut v = vector();
        v *= 2.0;
        assert_eq!(v, Vector3::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn a_double_on_the_left_scales_every_component() {
        assert_eq!(2.0 * vector(), Vector3::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn a_single_on_the_left_scales_every_component() {
        let v = Vector3::new(1.0_f32, 2.0, 3.0);
        assert_eq!(2.0_f32 * v, Vector3::new(2.0_f32, 4.0, 6.0));
    }

    #[test]
    fn dividing_by_a_number_scales_every_component_down() {
        assert_eq!(Vector3::new(2.0, 4.0, 6.0) / 2.0, vector());
    }

    #[test]
    fn div_assign_scales_every_component_down_in_place() {
        let mut v = Vector3::new(2.0, 4.0, 6.0);
        v /= 2.0;
        assert_eq!(v, vector());
    }

    #[test]
    fn summing_owned_vectors_adds_them() {
        assert_eq!(
            [vector(), vector()].into_iter().sum::<Vector3<f64>>(),
            vector() * 2.0
        );
    }

    #[test]
    fn summing_borrowed_vectors_adds_them() {
        assert_eq!(
            [vector(), vector()].iter().sum::<Vector3<f64>>(),
            vector() * 2.0
        );
    }

    #[test]
    fn norm_is_the_euclidean_length() {
        assert_eq!(Vector3::new(3.0, 4.0, 0.0).norm(), 5.0);
    }

    #[test]
    fn normalize_rescales_to_unit_length() {
        assert_eq!(
            Vector3::new(0.0, 4.0, 0.0).normalize(),
            Vector3::new(0.0, 1.0, 0.0)
        );
    }

    #[test]
    fn try_normalize_of_a_nonzero_vector_is_some() {
        let unit = Vector3::new(0.0, 4.0, 0.0).try_normalize();
        assert_eq!(unit.unwrap(), Vector3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn normalize_or_zero_rescales_a_nonzero_vector() {
        let unit = Vector3::new(0.0, 4.0, 0.0).normalize_or_zero();
        assert_eq!(unit, Vector3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn angle_between_perpendicular_vectors_is_a_quarter_turn() {
        assert!(close(
            Vector3::<f64>::X.angle_between(Vector3::Y),
            FRAC_PI_2
        ));
    }

    #[test]
    fn project_onto_yields_the_component_along_the_target() {
        let v = Vector3::new(2.0, 3.0, 0.0);
        assert_eq!(v.project_onto(Vector3::X), Vector3::new(2.0, 0.0, 0.0));
    }

    #[test]
    fn reject_from_yields_the_component_orthogonal_to_the_target() {
        let v = Vector3::new(2.0, 3.0, 0.0);
        assert_eq!(v.reject_from(Vector3::X), Vector3::new(0.0, 3.0, 0.0));
    }

    #[test]
    fn reflect_mirrors_the_vector_across_the_plane() {
        let v = Vector3::new(1.0, -2.0, 3.0);
        assert_eq!(v.reflect(Vector3::Y), Vector3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn lerp_interpolates_between_the_vectors() {
        let target = Vector3::new(4.0, 8.0, 12.0);
        assert_eq!(Vector3::ZERO.lerp(target, 0.25), vector());
    }

    #[test]
    fn abs_takes_the_magnitude_of_every_component() {
        assert_eq!(Vector3::new(-1.0, 2.0, -3.0).abs(), vector());
    }

    #[test]
    fn min_takes_the_smaller_of_each_pair_of_components() {
        let other = Vector3::new(4.0, 2.0, 3.0);
        assert_eq!(Vector3::new(1.0, 5.0, 3.0).min(other), vector());
    }

    #[test]
    fn max_takes_the_larger_of_each_pair_of_components() {
        let other = Vector3::new(4.0, 2.0, 3.0);
        assert_eq!(
            Vector3::new(1.0, 5.0, 3.0).max(other),
            Vector3::new(4.0, 5.0, 3.0)
        );
    }

    #[test]
    fn clamp_raises_a_component_below_the_interval_to_the_lower_bound() {
        let v = Vector3::new(-1.0, 2.0, 3.0);
        assert_eq!(
            v.clamp(Vector3::ZERO, Vector3::splat(9.0)),
            Vector3::new(0.0, 2.0, 3.0)
        );
    }

    #[test]
    fn clamp_leaves_a_component_inside_the_interval_unchanged() {
        assert_eq!(vector().clamp(Vector3::ZERO, Vector3::splat(9.0)), vector());
    }

    #[test]
    fn clamp_lowers_a_component_above_the_interval_to_the_upper_bound() {
        let v = Vector3::new(1.0, 2.0, 99.0);
        assert_eq!(
            v.clamp(Vector3::ZERO, Vector3::splat(9.0)),
            Vector3::new(1.0, 2.0, 9.0)
        );
    }

    #[test]
    fn midpoint_is_halfway_between_the_vectors() {
        assert_eq!(
            Vector3::ZERO.midpoint(Vector3::new(2.0, 4.0, 6.0)),
            vector()
        );
    }

    #[test]
    fn min_element_is_the_smallest_component() {
        assert_eq!(Vector3::new(3.0, -1.0, 2.0).min_element(), -1.0);
    }

    #[test]
    fn max_element_is_the_largest_component() {
        assert_eq!(Vector3::new(3.0, -1.0, 2.0).max_element(), 3.0);
    }

    #[test]
    fn element_sum_adds_the_three_components() {
        assert_eq!(vector().element_sum(), 6.0);
    }

    #[test]
    fn signum_takes_the_sign_of_every_component() {
        assert_eq!(
            Vector3::new(-2.0, 0.5, -3.0).signum(),
            Vector3::new(-1.0, 1.0, -1.0)
        );
    }

    #[test]
    fn copysign_keeps_the_magnitudes_and_takes_the_signs_of_its_argument() {
        let v = Vector3::new(3.0, -4.0, 5.0);
        let signs = Vector3::new(-1.0, 1.0, -1.0);
        assert_eq!(v.copysign(signs), Vector3::new(-3.0, 4.0, -5.0));
    }

    #[test]
    fn floor_rounds_every_component_toward_negative_infinity() {
        assert_eq!(
            Vector3::new(-2.5, 2.5, 3.0).floor(),
            Vector3::new(-3.0, 2.0, 3.0)
        );
    }

    #[test]
    fn ceil_rounds_every_component_toward_positive_infinity() {
        assert_eq!(
            Vector3::new(-2.5, 2.5, 3.0).ceil(),
            Vector3::new(-2.0, 3.0, 3.0)
        );
    }

    #[test]
    fn round_sends_a_half_away_from_zero() {
        assert_eq!(
            Vector3::new(2.5, 3.5, -2.5).round(),
            Vector3::new(3.0, 4.0, -3.0)
        );
    }

    #[test]
    fn round_ties_even_sends_a_half_to_the_even_integer() {
        let v = Vector3::new(2.5, 3.5, -2.5);
        assert_eq!(v.round_ties_even(), Vector3::new(2.0, 4.0, -2.0));
    }

    #[test]
    fn trunc_drops_the_fractional_part_of_every_component() {
        assert_eq!(
            Vector3::new(-2.75, 2.75, 3.0).trunc(),
            Vector3::new(-2.0, 2.0, 3.0)
        );
    }

    #[test]
    fn fract_keeps_the_fractional_part_of_every_component() {
        let v = Vector3::new(-2.75, 2.75, 3.0);
        assert_eq!(v.fract(), Vector3::new(-0.75, 0.75, 0.0));
    }

    #[test]
    fn div_euclid_is_the_euclidean_quotient_of_every_component() {
        let v = Vector3::new(7.0, -7.0, 0.0);
        assert_eq!(
            v.div_euclid(Vector3::splat(3.0)),
            Vector3::new(2.0, -3.0, 0.0)
        );
    }

    #[test]
    fn rem_euclid_is_nonnegative_for_a_negative_component() {
        let v = Vector3::new(-1.0, 7.0, 5.0);
        assert_eq!(
            v.rem_euclid(Vector3::splat(3.0)),
            Vector3::new(2.0, 1.0, 2.0)
        );
    }

    #[test]
    fn mul_add_scales_then_offsets_every_component() {
        let offset = Vector3::new(10.0, 20.0, 30.0);
        assert_eq!(
            vector().mul_add(2.0, offset),
            Vector3::new(12.0, 24.0, 36.0)
        );
    }

    #[test]
    fn hypot_combines_the_components_pairwise() {
        let legs = Vector3::new(4.0, 0.0, 5.0);
        let combined = Vector3::new(3.0, 4.0, 12.0).hypot(legs);
        assert!(vectors_close(combined, Vector3::new(5.0, 4.0, 13.0)));
    }

    #[test]
    fn is_nan_holds_when_a_component_is_not_a_number() {
        assert!(Vector3::new(1.0, f64::NAN, 3.0).is_nan());
    }

    #[test]
    fn is_infinite_holds_when_a_component_is_infinite() {
        assert!(Vector3::new(1.0, f64::INFINITY, 3.0).is_infinite());
    }

    #[test]
    fn is_finite_holds_when_every_component_is_finite() {
        assert!(vector().is_finite());
    }

    #[test]
    fn is_normal_holds_when_every_component_is_normal() {
        assert!(vector().is_normal());
    }

    #[test]
    fn is_subnormal_holds_when_a_component_is_subnormal() {
        assert!(Vector3::new(1.0, f64::MIN_POSITIVE / 2.0, 3.0).is_subnormal());
    }

    #[test]
    fn the_one_vector_has_every_component_set_to_one() {
        assert_eq!(Vector3::<f64>::ONE, Vector3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn the_x_axis_vector_points_along_the_first_component() {
        assert_eq!(Vector3::<f64>::X, Vector3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn the_y_axis_vector_points_along_the_second_component() {
        assert_eq!(Vector3::<f64>::Y, Vector3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn the_z_axis_vector_points_along_the_third_component() {
        assert_eq!(Vector3::<f64>::Z, Vector3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn dot_multiplies_and_sums_the_components() {
        assert_eq!(vector().dot(Vector3::new(4.0, 5.0, 6.0)), 32.0);
    }

    #[test]
    fn cross_returns_the_right_handed_product() {
        assert_eq!(Vector3::<f64>::X.cross(Vector3::Y), Vector3::Z);
    }

    #[test]
    fn norm_squared_is_the_sum_of_the_squared_components() {
        assert_eq!(vector().norm_squared(), 14.0);
    }

    #[test]
    fn recip_inverts_every_component() {
        let v = Vector3::new(2.0, 4.0, 8.0);
        assert_eq!(v.recip(), Vector3::new(0.5, 0.25, 0.125));
    }

    #[test]
    fn element_product_multiplies_the_three_components() {
        assert_eq!(Vector3::new(2.0, 3.0, 4.0).element_product(), 24.0);
    }

    #[test]
    fn normalize_of_a_length_vector_is_dimensionless() {
        let v = Vector3::new(length(0.0), length(4.0), length(0.0));
        assert_eq!(v.normalize(), Vector3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn try_normalize_of_a_length_vector_is_dimensionless() {
        let v = Vector3::new(length(0.0), length(4.0), length(0.0));
        assert_eq!(v.try_normalize().unwrap(), Vector3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn normalize_or_zero_of_a_length_vector_is_dimensionless() {
        let v = Vector3::new(length(0.0), length(4.0), length(0.0));
        assert_eq!(v.normalize_or_zero(), Vector3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn angle_between_length_vectors_is_dimensionless() {
        let a = Vector3::new(length(1.0), length(0.0), length(0.0));
        let b = Vector3::new(length(0.0), length(2.0), length(0.0));
        assert!(close(a.angle_between(b), FRAC_PI_2));
    }

    #[test]
    fn reflect_takes_a_dimensionless_normal() {
        let v = Vector3::new(length(1.0), length(-2.0), length(3.0));
        let expected = Vector3::new(length(1.0), length(2.0), length(3.0));
        assert_eq!(v.reflect(Vector3::Y), expected);
    }

    #[test]
    fn lerp_of_length_vectors_takes_a_dimensionless_factor() {
        let target = Vector3::new(length(2.0), length(4.0), length(6.0));
        let expected = Vector3::new(length(1.0), length(2.0), length(3.0));
        assert_eq!(Vector3::ZERO.lerp(target, 0.5), expected);
    }

    #[test]
    fn multiplying_a_dimensionless_vector_by_a_length_carries_the_dimension() {
        let expected = Vector3::new(length(2.0), length(4.0), length(6.0));
        assert_eq!(vector() * length(2.0), expected);
    }

    #[test]
    fn a_double_on_the_left_scales_a_length_vector() {
        let v = Vector3::new(length(1.0), length(2.0), length(3.0));
        let expected = Vector3::new(length(2.0), length(4.0), length(6.0));
        assert_eq!(2.0 * v, expected);
    }

    #[test]
    fn dividing_a_length_vector_by_a_length_leaves_it_dimensionless() {
        let v = Vector3::new(length(2.0), length(4.0), length(6.0));
        assert_eq!(v / length(2.0), vector());
    }

    #[test]
    fn signum_of_a_length_vector_is_dimensionless() {
        let v = Vector3::new(length(-2.0), length(0.5), length(-3.0));
        assert_eq!(v.signum(), Vector3::new(-1.0, 1.0, -1.0));
    }

    #[test]
    fn div_euclid_of_length_vectors_is_dimensionless() {
        let v = Vector3::new(length(7.0), length(-7.0), length(0.0));
        assert_eq!(
            v.div_euclid(Vector3::splat(length(3.0))),
            Vector3::new(2.0, -3.0, 0.0)
        );
    }

    #[test]
    fn mul_add_of_length_vectors_takes_a_dimensionless_factor() {
        let v = Vector3::new(length(1.0), length(2.0), length(3.0));
        let offset = Vector3::new(length(10.0), length(20.0), length(30.0));
        let expected = Vector3::new(length(12.0), length(24.0), length(36.0));
        assert_eq!(v.mul_add(2.0, offset), expected);
    }

    #[test]
    fn hypot_of_length_vectors_stays_a_length() {
        let v = Vector3::new(length(3.0), length(4.0), length(12.0));
        let legs = Vector3::new(length(4.0), length(0.0), length(5.0));
        let combined = v.hypot(legs).map(Quantity::value);
        assert!(vectors_close(combined, Vector3::new(5.0, 4.0, 13.0)));
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn from_slice_shorter_than_three_panics() {
        let _ = Vector3::from_slice(&[1.0_f64, 2.0]);
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn indexing_out_of_bounds_panics() {
        let _ = vector()[3];
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn index_mut_out_of_bounds_panics() {
        vector()[3] = 0.0;
    }

    #[test]
    #[should_panic(expected = "min > max")]
    fn clamp_with_an_inverted_interval_panics() {
        let _ = vector().clamp(Vector3::splat(9.0), Vector3::ZERO);
    }

    #[test]
    fn try_normalize_of_a_non_finite_vector_is_none() {
        let v = Vector3::new(f64::INFINITY, 0.0, 0.0);
        assert!(v.try_normalize().is_none());
    }

    #[test]
    fn min_ignores_a_not_a_number_component() {
        assert_eq!(Vector3::splat(f64::NAN).min(vector()), vector());
    }

    #[test]
    fn max_ignores_a_not_a_number_component() {
        assert_eq!(Vector3::splat(f64::NAN).max(vector()), vector());
    }

    #[test]
    fn is_nan_does_not_hold_when_no_component_is_not_a_number() {
        assert!(!vector().is_nan());
    }

    #[test]
    fn is_infinite_does_not_hold_when_every_component_is_finite() {
        assert!(!vector().is_infinite());
    }

    #[test]
    fn is_finite_does_not_hold_when_a_component_is_infinite() {
        assert!(!Vector3::new(1.0, f64::INFINITY, 3.0).is_finite());
    }

    #[test]
    fn is_normal_does_not_hold_when_a_component_is_zero() {
        assert!(!Vector3::new(1.0, 0.0, 3.0).is_normal());
    }

    #[test]
    fn is_subnormal_does_not_hold_when_every_component_is_normal() {
        assert!(!vector().is_subnormal());
    }

    #[test]
    fn from_slice_of_exactly_three_elements_yields_the_vector() {
        assert_eq!(Vector3::from_slice(&[1.0, 2.0, 3.0]), vector());
    }

    #[test]
    fn floor_leaves_an_integer_component_unchanged() {
        assert_eq!(vector().floor(), vector());
    }

    #[test]
    fn ceil_leaves_an_integer_component_unchanged() {
        assert_eq!(vector().ceil(), vector());
    }

    #[test]
    fn fract_of_an_integer_component_is_zero() {
        assert_eq!(vector().fract(), Vector3::ZERO);
    }

    #[test]
    fn signum_of_a_negative_zero_component_is_minus_one() {
        assert_eq!(Vector3::new(-0.0, 1.0, 1.0).signum().x, -1.0);
    }

    #[test]
    fn angle_between_antiparallel_vectors_is_a_half_turn() {
        assert!(close(Vector3::<f64>::X.angle_between(-Vector3::X), PI));
    }

    #[test]
    fn lerp_beyond_one_extrapolates_past_the_ending_vector() {
        let target = Vector3::new(4.0, 8.0, 12.0);
        assert_eq!(
            Vector3::ZERO.lerp(target, 2.0),
            Vector3::new(8.0, 16.0, 24.0)
        );
    }

    #[test]
    fn array_roundtrip_preserves_the_vector() {
        assert_eq!(Vector3::from_array(vector().to_array()), vector());
    }

    #[test]
    fn addition_is_commutative() {
        let other = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(vector() + other, other + vector());
    }

    #[test]
    fn the_zero_vector_is_the_additive_identity() {
        assert_eq!(vector() + Vector3::ZERO, vector());
    }

    #[test]
    fn negation_is_the_additive_inverse() {
        assert_eq!(vector() + -vector(), Vector3::ZERO);
    }

    #[test]
    fn scaling_on_either_side_gives_the_same_vector() {
        assert_eq!(2.0 * vector(), vector() * 2.0);
    }

    #[test]
    fn dot_product_is_symmetric() {
        let other = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(vector().dot(other), other.dot(vector()));
    }

    #[test]
    fn cross_product_is_anticommutative() {
        let other = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(vector().cross(other), -other.cross(vector()));
    }

    #[test]
    fn cross_product_is_orthogonal_to_both_factors() {
        let other = Vector3::new(4.0, 5.0, 6.0);
        let normal = vector().cross(other);
        assert!(close(normal.dot(vector()), 0.0) && close(normal.dot(other), 0.0));
    }

    #[test]
    fn norm_squared_is_the_square_of_the_norm() {
        let norm = vector().norm();
        assert!(close(vector().norm_squared(), norm * norm));
    }

    #[test]
    fn norm_holds_above_the_squared_range() {
        let scale = 2.0_f64.powi(600);
        let v = Vector3::new(3.0 * scale, 4.0 * scale, 0.0);
        assert!(close(v.norm() / scale, 5.0));
    }

    #[test]
    fn norm_holds_below_the_squared_range() {
        let scale = 2.0_f64.powi(-600);
        let v = Vector3::new(3.0 * scale, 4.0 * scale, 0.0);
        assert!(close(v.norm() / scale, 5.0));
    }

    #[test]
    fn project_onto_holds_above_the_squared_range() {
        let scale = 2.0_f64.powi(600);
        let v = Vector3::new(3.0 * scale, 4.0 * scale, 0.0);
        let projected = v.project_onto(Vector3::new(scale, 0.0, 0.0));
        assert!(close(projected.x / scale, 3.0) && projected.y == 0.0);
    }

    #[test]
    fn reject_from_holds_above_the_squared_range() {
        let scale = 2.0_f64.powi(600);
        let v = Vector3::new(3.0 * scale, 4.0 * scale, 0.0);
        let rejected = v.reject_from(Vector3::new(scale, 0.0, 0.0));
        assert!(close(rejected.y / scale, 4.0) && rejected.x == 0.0);
    }

    #[test]
    fn norm_of_a_vector_with_an_infinite_component_is_infinite() {
        assert_eq!(Vector3::new(f64::INFINITY, 1.0, 2.0).norm(), f64::INFINITY);
    }

    #[test]
    fn norm_of_a_vector_with_a_not_a_number_component_is_not_a_number() {
        assert!(Vector3::new(f64::NAN, 1.0, 2.0).norm().is_nan());
    }

    #[test]
    fn angle_between_is_unchanged_by_a_scale_beyond_the_squared_range() {
        let scale = 2.0_f64.powi(600);
        let far =
            Vector3::new(scale, 2.0 * scale, 0.0).angle_between(Vector3::new(scale, 0.0, 0.0));
        let near = Vector3::new(1.0, 2.0, 0.0).angle_between(Vector3::X);
        assert!(close(far, near));
    }

    #[test]
    fn angle_between_is_unchanged_by_a_scale_below_the_squared_range() {
        let scale = 2.0_f64.powi(-600);
        let near =
            Vector3::new(scale, 2.0 * scale, 0.0).angle_between(Vector3::new(scale, 0.0, 0.0));
        assert!(close(
            near,
            Vector3::new(1.0, 2.0, 0.0).angle_between(Vector3::X)
        ));
    }

    #[test]
    fn lerp_holds_when_the_displacement_between_the_ends_overflows() {
        let start = Vector3::splat(-f64::MAX);
        assert_eq!(start.lerp(Vector3::splat(f64::MAX), 0.5), Vector3::ZERO);
    }

    #[test]
    fn a_vector_with_a_zero_component_is_not_normal() {
        assert!(!Vector3::<f64>::X.is_normal());
    }

    #[test]
    fn normalizing_yields_a_vector_of_unit_norm() {
        assert!(close(Vector3::new(3.0, 4.0, 12.0).normalize().norm(), 1.0));
    }

    #[test]
    fn projection_and_rejection_reconstruct_the_vector() {
        let onto = Vector3::new(1.0, 1.0, 0.0);
        let reconstructed = vector().project_onto(onto) + vector().reject_from(onto);
        assert!(vectors_close(reconstructed, vector()));
    }

    #[test]
    fn reflecting_twice_returns_the_vector() {
        let v = Vector3::new(1.0, -2.0, 3.0);
        assert_eq!(v.reflect(Vector3::Y).reflect(Vector3::Y), v);
    }

    #[test]
    fn truncating_and_the_fractional_part_reconstruct_the_vector() {
        let v = Vector3::new(-2.75, 2.75, 3.0);
        assert_eq!(v.trunc() + v.fract(), v);
    }

    #[test]
    fn min_of_a_vector_with_itself_is_that_vector() {
        assert_eq!(vector().min(vector()), vector());
    }

    #[test]
    fn max_of_a_vector_with_itself_is_that_vector() {
        assert_eq!(vector().max(vector()), vector());
    }

    #[test]
    fn clamp_agrees_with_taking_the_maximum_then_the_minimum() {
        let (lo, hi) = (Vector3::splat(1.5), Vector3::splat(2.5));
        assert_eq!(vector().clamp(lo, hi), vector().max(lo).min(hi));
    }

    #[test]
    fn midpoint_agrees_with_interpolating_halfway() {
        let other = Vector3::new(4.0, 8.0, 12.0);
        assert_eq!(vector().midpoint(other), vector().lerp(other, 0.5));
    }

    #[test]
    fn recip_is_its_own_inverse() {
        let v = Vector3::new(2.0, 4.0, 8.0);
        assert_eq!(v.recip().recip(), v);
    }

    #[test]
    fn vectors_are_equal_exactly_when_all_components_match() {
        assert_eq!(vector(), Vector3::new(1.0, 2.0, 3.0));
        assert_ne!(vector(), Vector3::new(1.0, 2.0, 4.0));
    }
}

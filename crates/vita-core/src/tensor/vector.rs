use core::iter::Sum;
use core::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

use crate::Scalar;

/// A vector of three components `x`, `y`, and `z`.
pub struct Vector3<T> {
    /// The first component.
    pub x: T,
    /// The second component.
    pub y: T,
    /// The third component.
    pub z: T,
}

impl<T: ::core::marker::Copy> ::core::marker::Copy for Vector3<T> {}

impl<T: ::core::clone::Clone> ::core::clone::Clone for Vector3<T> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.clone(),
            z: self.z.clone(),
        }
    }
}

impl<T: ::core::fmt::Debug> ::core::fmt::Debug for Vector3<T> {
    fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
        f.debug_struct("Vector3")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("z", &self.z)
            .finish()
    }
}

impl<T: ::core::cmp::PartialEq> ::core::cmp::PartialEq for Vector3<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
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

impl<T: Default> Default for Vector3<T> {
    /// Returns the zero vector.
    #[inline]
    fn default() -> Self {
        Self {
            x: T::default(),
            y: T::default(),
            z: T::default(),
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

impl<T, S: Scalar> Mul<S> for Vector3<T>
where
    T: Mul<S, Output = T>,
{
    type Output = Self;
    /// Scales every component by the scalar `rhs`, preserving the element's unit.
    #[inline]
    fn mul(self, rhs: S) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<T, S: Scalar> MulAssign<S> for Vector3<T>
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

impl<T, S: Scalar> Div<S> for Vector3<T>
where
    T: Div<S, Output = T>,
{
    type Output = Self;
    /// Divides every component by the scalar `rhs`, preserving the element's unit.
    #[inline]
    fn div(self, rhs: S) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl<T, S: Scalar> DivAssign<S> for Vector3<T>
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

impl<T: Add<Output = T> + Default + Copy> Sum for Vector3<T> {
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

impl<T: Add<Output = T>> Vector3<T> {
    /// Returns the sum of the three components, `x + y + z`.
    #[inline]
    pub fn element_sum(self) -> T {
        self.x + self.y + self.z
    }
}

impl<T: Add<Output = T> + Sub<Output = T> + Copy> Vector3<T> {
    /// Linearly interpolates from `self` toward `rhs` by the factor `t`.
    ///
    /// `t == 0` yields `self`, `t == 1` yields `rhs`.
    #[inline]
    pub fn lerp<S: Scalar>(self, rhs: Self, t: S) -> Self
    where
        T: Mul<S, Output = T>,
    {
        self + (rhs - self) * t
    }
}

impl<V: Scalar> Vector3<V> {
    /// The zero vector, `(0, 0, 0)`.
    pub const ZERO: Self = Self::new(V::ZERO, V::ZERO, V::ZERO);

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

    /// Returns the Euclidean norm (length) of the vector.
    #[inline]
    pub fn norm(self) -> V {
        self.norm_squared().sqrt()
    }

    /// Returns `self` rescaled to unit length.
    ///
    /// Yields a non-finite vector when the norm is zero or non-finite; use
    /// [`try_normalize`][Self::try_normalize] or
    /// [`normalize_or_zero`][Self::normalize_or_zero] to handle those cases.
    #[inline]
    pub fn normalize(self) -> Self {
        self / self.norm()
    }

    /// Returns `self` rescaled to unit length, or `None` if the result would
    /// not be finite (e.g. for the zero vector).
    #[inline]
    pub fn try_normalize(self) -> Option<Self> {
        let norm = self.norm();
        if norm.is_finite() && norm > V::ZERO {
            Some(self / norm)
        } else {
            None
        }
    }

    /// Returns `self` rescaled to unit length, or the zero vector if the
    /// result would not be finite.
    #[inline]
    pub fn normalize_or_zero(self) -> Self {
        self.try_normalize().unwrap_or(Self::ZERO)
    }

    /// Returns `true` if the vector is of unit length within a small
    /// tolerance (`2.0e-4`).
    #[inline]
    pub fn is_normalized(self) -> bool {
        (self.norm_squared() - V::ONE).abs() <= V::from_f64(2.0e-4)
    }

    /// Returns the unsigned angle between `self` and `rhs`, in radians within
    /// `[0, π]`.
    ///
    /// The computation is numerically stable across the whole range, including
    /// near-parallel and near-antiparallel inputs.
    #[inline]
    pub fn angle_between(self, rhs: Self) -> V {
        self.cross(rhs).norm().atan2(self.dot(rhs))
    }

    /// Returns the vector projection of `self` onto `onto`.
    #[inline]
    pub fn project_onto(self, onto: Self) -> Self {
        onto * (self.dot(onto) / onto.norm_squared())
    }

    /// Returns the component of `self` orthogonal to `from`.
    #[inline]
    pub fn reject_from(self, from: Self) -> Self {
        self - self.project_onto(from)
    }

    /// Reflects `self` across the plane through the origin with unit normal
    /// `normal`.
    #[inline]
    pub fn reflect(self, normal: Self) -> Self {
        self - normal * (self.dot(normal) * (V::ONE + V::ONE))
    }

    /// Returns the component-wise reciprocal `(1/x, 1/y, 1/z)`.
    #[inline]
    pub fn recip(self) -> Self {
        Self::new(self.x.recip(), self.y.recip(), self.z.recip())
    }

    /// Returns the product of the three components, `x * y * z`.
    #[inline]
    pub fn element_product(self) -> V {
        self.x * self.y * self.z
    }

    /// Returns the component-wise absolute value.
    #[inline]
    pub fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// Returns the component-wise minimum of `self` and `rhs`.
    #[inline]
    pub fn min(self, rhs: Self) -> Self {
        Self::new(self.x.min(rhs.x), self.y.min(rhs.y), self.z.min(rhs.z))
    }

    /// Returns the component-wise maximum of `self` and `rhs`.
    #[inline]
    pub fn max(self, rhs: Self) -> Self {
        Self::new(self.x.max(rhs.x), self.y.max(rhs.y), self.z.max(rhs.z))
    }

    /// Restricts every component to the interval `[min, max]`.
    ///
    /// # Panics
    ///
    /// Panics if any component of `min` exceeds the corresponding component of
    /// `max`.
    #[inline]
    pub fn clamp(self, min: Self, max: Self) -> Self {
        Self::new(
            self.x.clamp(min.x, max.x),
            self.y.clamp(min.y, max.y),
            self.z.clamp(min.z, max.z),
        )
    }

    /// Returns the smallest of the three components.
    #[inline]
    pub fn min_element(self) -> V {
        self.x.min(self.y).min(self.z)
    }

    /// Returns the largest of the three components.
    #[inline]
    pub fn max_element(self) -> V {
        self.x.max(self.y).max(self.z)
    }

    /// Returns the component-wise floor.
    #[inline]
    pub fn floor(self) -> Self {
        Self::new(self.x.floor(), self.y.floor(), self.z.floor())
    }

    /// Returns the component-wise ceiling.
    #[inline]
    pub fn ceil(self) -> Self {
        Self::new(self.x.ceil(), self.y.ceil(), self.z.ceil())
    }

    /// Returns the component-wise nearest integer, rounding halves away from
    /// zero.
    #[inline]
    pub fn round(self) -> Self {
        Self::new(self.x.round(), self.y.round(), self.z.round())
    }

    /// Returns the component-wise nearest integer, rounding halves to even.
    #[inline]
    pub fn round_ties_even(self) -> Self {
        Self::new(
            self.x.round_ties_even(),
            self.y.round_ties_even(),
            self.z.round_ties_even(),
        )
    }

    /// Returns the component-wise truncation toward zero.
    #[inline]
    pub fn trunc(self) -> Self {
        Self::new(self.x.trunc(), self.y.trunc(), self.z.trunc())
    }

    /// Returns the component-wise fractional part.
    #[inline]
    pub fn fract(self) -> Self {
        Self::new(self.x.fract(), self.y.fract(), self.z.fract())
    }

    /// Returns a vector with the magnitudes of `self` and the component-wise
    /// signs of `sign`.
    #[inline]
    pub fn copysign(self, sign: Self) -> Self {
        Self::new(
            self.x.copysign(sign.x),
            self.y.copysign(sign.y),
            self.z.copysign(sign.z),
        )
    }

    /// Returns the component-wise sign, each `1`, `-1`, or `NaN`.
    #[inline]
    pub fn signum(self) -> Self {
        Self::new(self.x.signum(), self.y.signum(), self.z.signum())
    }

    /// Returns the component-wise least nonnegative remainder against `rhs`.
    #[inline]
    pub fn rem_euclid(self, rhs: Self) -> Self {
        Self::new(
            self.x.rem_euclid(rhs.x),
            self.y.rem_euclid(rhs.y),
            self.z.rem_euclid(rhs.z),
        )
    }

    /// Returns `true` if every component is finite.
    #[inline]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Returns `true` if any component is positive or negative infinity.
    #[inline]
    pub fn is_infinite(self) -> bool {
        self.x.is_infinite() || self.y.is_infinite() || self.z.is_infinite()
    }

    /// Returns `true` if any component is `NaN`.
    #[inline]
    pub fn is_nan(self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use core::f64::consts::{FRAC_PI_2, PI};

    const TOL: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL
    }

    #[test]
    fn new_sets_the_three_components() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(v.x, 1.0);
        assert_eq!(v.y, 2.0);
        assert_eq!(v.z, 3.0);
    }

    #[test]
    fn splat_repeats_one_value_across_all_components() {
        assert_eq!(Vector3::splat(5.0), Vector3::new(5.0, 5.0, 5.0));
    }

    #[test]
    fn from_array_takes_components_in_order() {
        assert_eq!(
            Vector3::from_array([1.0, 2.0, 3.0]),
            Vector3::new(1.0, 2.0, 3.0)
        );
    }

    #[test]
    fn to_array_returns_components_in_order() {
        assert_eq!(Vector3::new(1.0, 2.0, 3.0).to_array(), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn from_slice_reads_the_first_three_elements() {
        assert_eq!(
            Vector3::from_slice(&[1.0, 2.0, 3.0, 4.0]),
            Vector3::new(1.0, 2.0, 3.0),
        );
    }

    #[test]
    #[should_panic]
    fn from_slice_shorter_than_three_panics() {
        let _ = Vector3::from_slice(&[1.0_f64, 2.0]);
    }

    #[test]
    fn with_x_replaces_only_the_x_component() {
        assert_eq!(
            Vector3::new(1.0, 2.0, 3.0).with_x(9.0),
            Vector3::new(9.0, 2.0, 3.0),
        );
    }

    #[test]
    fn with_y_replaces_only_the_y_component() {
        assert_eq!(
            Vector3::new(1.0, 2.0, 3.0).with_y(9.0),
            Vector3::new(1.0, 9.0, 3.0),
        );
    }

    #[test]
    fn with_z_replaces_only_the_z_component() {
        assert_eq!(
            Vector3::new(1.0, 2.0, 3.0).with_z(9.0),
            Vector3::new(1.0, 2.0, 9.0),
        );
    }

    #[test]
    fn index_reads_the_component_at_each_position() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);
    }

    #[test]
    fn index_mut_writes_the_component_at_each_position() {
        let mut v = Vector3::new(0.0, 0.0, 0.0);
        v[0] = 1.0;
        v[1] = 2.0;
        v[2] = 3.0;
        assert_eq!(v, Vector3::new(1.0, 2.0, 3.0));
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn index_out_of_bounds_panics() {
        let _ = Vector3::new(1.0_f64, 2.0, 3.0)[3];
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn index_mut_out_of_bounds_panics() {
        let mut v = Vector3::new(1.0_f64, 2.0, 3.0);
        v[3] = 0.0;
    }

    #[test]
    fn map_applies_the_function_to_every_component() {
        assert_eq!(
            Vector3::new(1.0, 2.0, 3.0).map(|c| c * 2.0),
            Vector3::new(2.0, 4.0, 6.0),
        );
    }

    #[test]
    fn zip_map_combines_the_two_vectors_component_wise() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a.zip_map(b, |x, y| x * y), Vector3::new(4.0, 10.0, 18.0));
    }

    #[test]
    fn neg_negates_every_component() {
        assert_eq!(-Vector3::new(1.0, -2.0, 3.0), Vector3::new(-1.0, 2.0, -3.0));
    }

    #[test]
    fn add_sums_components_pairwise() {
        assert_eq!(
            Vector3::new(1.0, 2.0, 3.0) + Vector3::new(4.0, 5.0, 6.0),
            Vector3::new(5.0, 7.0, 9.0),
        );
    }

    #[test]
    fn sub_subtracts_components_pairwise() {
        assert_eq!(
            Vector3::new(5.0, 7.0, 9.0) - Vector3::new(4.0, 5.0, 6.0),
            Vector3::new(1.0, 2.0, 3.0),
        );
    }

    #[test]
    fn add_assign_adds_in_place() {
        let mut v = Vector3::new(1.0, 2.0, 3.0);
        v += Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(v, Vector3::new(5.0, 7.0, 9.0));
    }

    #[test]
    fn sub_assign_subtracts_in_place() {
        let mut v = Vector3::new(5.0, 7.0, 9.0);
        v -= Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(v, Vector3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn mul_scales_every_component() {
        assert_eq!(
            Vector3::new(1.0, 2.0, 3.0) * 2.0,
            Vector3::new(2.0, 4.0, 6.0)
        );
    }

    #[test]
    fn div_divides_every_component() {
        assert_eq!(
            Vector3::new(2.0, 4.0, 6.0) / 2.0,
            Vector3::new(1.0, 2.0, 3.0)
        );
    }

    #[test]
    fn mul_assign_scales_in_place() {
        let mut v = Vector3::new(1.0, 2.0, 3.0);
        v *= 2.0;
        assert_eq!(v, Vector3::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn div_assign_divides_in_place() {
        let mut v = Vector3::new(2.0, 4.0, 6.0);
        v /= 2.0;
        assert_eq!(v, Vector3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn element_sum_adds_the_three_components() {
        assert_eq!(Vector3::new(1.0, 2.0, 3.0).element_sum(), 6.0);
    }

    #[test]
    fn element_product_multiplies_the_three_components() {
        assert_eq!(Vector3::new(2.0, 3.0, 4.0).element_product(), 24.0);
    }

    #[test]
    fn sum_of_no_vectors_is_the_zero_vector() {
        let total: Vector3<f64> = core::iter::empty::<Vector3<f64>>().sum();
        assert_eq!(total, Vector3::ZERO);
    }

    #[test]
    fn sum_folds_owned_vectors() {
        let total: Vector3<f64> = [Vector3::X, Vector3::Y, Vector3::Z].into_iter().sum();
        assert_eq!(total, Vector3::ONE);
    }

    #[test]
    fn sum_folds_borrowed_vectors() {
        let vectors = [Vector3::<f64>::X, Vector3::Y];
        let total: Vector3<f64> = vectors.iter().sum();
        assert_eq!(total, Vector3::new(1.0, 1.0, 0.0));
    }

    #[test]
    fn lerp_at_one_half_returns_the_midpoint() {
        let a = Vector3::new(0.0, 0.0, 0.0);
        let b = Vector3::new(10.0, 20.0, 30.0);
        assert_eq!(a.lerp(b, 0.5), Vector3::new(5.0, 10.0, 15.0));
    }

    #[test]
    fn lerp_at_the_endpoints_returns_each_bound() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a.lerp(b, 0.0), a);
        assert_eq!(a.lerp(b, 1.0), b);
    }

    #[test]
    fn default_is_the_zero_vector() {
        assert_eq!(Vector3::<f64>::default(), Vector3::ZERO);
    }

    #[test]
    fn zero_has_all_components_zero() {
        assert_eq!(Vector3::<f64>::ZERO, Vector3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn one_has_all_components_one() {
        assert_eq!(Vector3::<f64>::ONE, Vector3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn the_basis_constants_are_the_unit_axes() {
        assert_eq!(Vector3::<f64>::X, Vector3::new(1.0, 0.0, 0.0));
        assert_eq!(Vector3::<f64>::Y, Vector3::new(0.0, 1.0, 0.0));
        assert_eq!(Vector3::<f64>::Z, Vector3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn dot_multiplies_and_sums_the_components() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a.dot(b), 32.0);
    }

    #[test]
    fn cross_returns_the_right_hand_product() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a.cross(b), Vector3::new(-3.0, 6.0, -3.0));
    }

    #[test]
    fn norm_squared_is_the_sum_of_squares() {
        assert_eq!(Vector3::new(2.0, 3.0, 6.0).norm_squared(), 49.0);
    }

    #[test]
    fn norm_is_the_euclidean_length() {
        assert_eq!(Vector3::new(2.0, 3.0, 6.0).norm(), 7.0);
    }

    #[test]
    fn normalize_rescales_to_unit_length_in_the_same_direction() {
        assert_eq!(
            Vector3::new(2.0, 3.0, 6.0).normalize(),
            Vector3::new(2.0 / 7.0, 3.0 / 7.0, 6.0 / 7.0),
        );
    }

    #[test]
    fn normalize_of_the_zero_vector_is_not_finite() {
        assert!(!Vector3::<f64>::ZERO.normalize().is_finite());
    }

    #[test]
    fn try_normalize_of_a_nonzero_vector_is_some() {
        assert_eq!(
            Vector3::new(2.0, 3.0, 6.0).try_normalize(),
            Some(Vector3::new(2.0 / 7.0, 3.0 / 7.0, 6.0 / 7.0)),
        );
    }

    #[test]
    fn try_normalize_of_the_zero_vector_is_none() {
        assert_eq!(Vector3::<f64>::ZERO.try_normalize(), None);
    }

    #[test]
    fn normalize_or_zero_rescales_a_nonzero_vector() {
        assert_eq!(
            Vector3::new(2.0, 3.0, 6.0).normalize_or_zero(),
            Vector3::new(2.0 / 7.0, 3.0 / 7.0, 6.0 / 7.0),
        );
    }

    #[test]
    fn normalize_or_zero_of_the_zero_vector_is_zero() {
        assert_eq!(Vector3::<f64>::ZERO.normalize_or_zero(), Vector3::ZERO);
    }

    #[test]
    fn is_normalized_is_true_for_a_unit_vector() {
        assert!(Vector3::<f64>::X.is_normalized());
    }

    #[test]
    fn is_normalized_is_false_for_a_non_unit_vector() {
        assert!(!Vector3::new(2.0, 0.0, 0.0).is_normalized());
    }

    #[test]
    fn angle_between_perpendicular_vectors_is_a_right_angle() {
        let x = Vector3::<f64>::X;
        assert!(close(x.angle_between(Vector3::Y), FRAC_PI_2));
    }

    #[test]
    fn angle_between_parallel_vectors_is_zero() {
        let x = Vector3::<f64>::X;
        assert!(close(x.angle_between(x), 0.0));
    }

    #[test]
    fn angle_between_antiparallel_vectors_is_pi() {
        let x = Vector3::<f64>::X;
        assert!(close(x.angle_between(-x), PI));
    }

    #[test]
    fn project_onto_returns_the_parallel_component() {
        let v = Vector3::new(2.0, 3.0, 0.0);
        let axis = Vector3::new(1.0, 0.0, 0.0);
        assert_eq!(v.project_onto(axis), Vector3::new(2.0, 0.0, 0.0));
    }

    #[test]
    fn reject_from_returns_the_orthogonal_component() {
        let v = Vector3::new(2.0, 3.0, 0.0);
        let axis = Vector3::new(1.0, 0.0, 0.0);
        assert_eq!(v.reject_from(axis), Vector3::new(0.0, 3.0, 0.0));
    }

    #[test]
    fn reflect_mirrors_across_the_plane_normal() {
        let v = Vector3::new(1.0, -1.0, 0.0);
        let normal = Vector3::new(0.0, 1.0, 0.0);
        assert_eq!(v.reflect(normal), Vector3::new(1.0, 1.0, 0.0));
    }

    #[test]
    fn recip_inverts_every_component() {
        assert_eq!(
            Vector3::new(2.0, 4.0, 8.0).recip(),
            Vector3::new(0.5, 0.25, 0.125),
        );
    }

    #[test]
    fn abs_takes_the_magnitude_of_every_component() {
        assert_eq!(
            Vector3::new(-1.0, 2.0, -3.0).abs(),
            Vector3::new(1.0, 2.0, 3.0),
        );
    }

    #[test]
    fn min_selects_the_smaller_of_each_component() {
        let a = Vector3::new(1.0, 5.0, 2.0);
        let b = Vector3::new(4.0, 2.0, 3.0);
        assert_eq!(a.min(b), Vector3::new(1.0, 2.0, 2.0));
    }

    #[test]
    fn max_selects_the_larger_of_each_component() {
        let a = Vector3::new(1.0, 5.0, 2.0);
        let b = Vector3::new(4.0, 2.0, 3.0);
        assert_eq!(a.max(b), Vector3::new(4.0, 5.0, 3.0));
    }

    #[test]
    fn clamp_restricts_each_component_to_the_range() {
        let v = Vector3::new(-1.0, 5.0, 2.0);
        let lo = Vector3::new(0.0, 0.0, 0.0);
        let hi = Vector3::new(3.0, 3.0, 3.0);
        assert_eq!(v.clamp(lo, hi), Vector3::new(0.0, 3.0, 2.0));
    }

    #[test]
    #[should_panic]
    fn clamp_with_min_above_max_panics() {
        let _ = Vector3::new(1.0_f64, 2.0, 3.0).clamp(Vector3::splat(5.0), Vector3::splat(0.0));
    }

    #[test]
    fn min_element_returns_the_smallest_component() {
        assert_eq!(Vector3::new(3.0, 1.0, 2.0).min_element(), 1.0);
    }

    #[test]
    fn max_element_returns_the_largest_component() {
        assert_eq!(Vector3::new(3.0, 1.0, 2.0).max_element(), 3.0);
    }

    #[test]
    fn floor_rounds_each_component_down() {
        assert_eq!(
            Vector3::new(1.5, -1.5, 2.9).floor(),
            Vector3::new(1.0, -2.0, 2.0),
        );
    }

    #[test]
    fn ceil_rounds_each_component_up() {
        assert_eq!(
            Vector3::new(1.1, -1.9, 2.5).ceil(),
            Vector3::new(2.0, -1.0, 3.0),
        );
    }

    #[test]
    fn round_rounds_halves_away_from_zero() {
        assert_eq!(
            Vector3::new(0.5, 2.5, -0.5).round(),
            Vector3::new(1.0, 3.0, -1.0),
        );
    }

    #[test]
    fn round_ties_even_rounds_halves_to_even() {
        assert_eq!(
            Vector3::new(0.5, 2.5, -1.5).round_ties_even(),
            Vector3::new(0.0, 2.0, -2.0),
        );
    }

    #[test]
    fn trunc_discards_each_fractional_part() {
        assert_eq!(
            Vector3::new(1.7, -1.7, 2.2).trunc(),
            Vector3::new(1.0, -1.0, 2.0),
        );
    }

    #[test]
    fn fract_keeps_each_fractional_part() {
        assert_eq!(
            Vector3::new(1.25, -1.25, 2.5).fract(),
            Vector3::new(0.25, -0.25, 0.5),
        );
    }

    #[test]
    fn copysign_takes_each_sign_from_the_argument() {
        let magnitudes = Vector3::new(1.0, 2.0, 3.0);
        let signs = Vector3::new(-1.0, 1.0, -5.0);
        assert_eq!(magnitudes.copysign(signs), Vector3::new(-1.0, 2.0, -3.0));
    }

    #[test]
    fn signum_returns_the_sign_of_each_component() {
        assert_eq!(
            Vector3::new(-3.0, 4.0, -5.0).signum(),
            Vector3::new(-1.0, 1.0, -1.0),
        );
    }

    #[test]
    fn rem_euclid_returns_the_nonnegative_remainder() {
        let v = Vector3::new(7.0, -1.0, 8.0);
        let m = Vector3::new(3.0, 3.0, 3.0);
        assert_eq!(v.rem_euclid(m), Vector3::new(1.0, 2.0, 2.0));
    }

    #[test]
    fn is_finite_is_true_when_all_components_are_finite() {
        assert!(Vector3::new(1.0, 2.0, 3.0).is_finite());
    }

    #[test]
    fn is_finite_is_false_when_a_component_is_infinite() {
        assert!(!Vector3::new(1.0, f64::INFINITY, 3.0).is_finite());
    }

    #[test]
    fn is_infinite_is_true_when_a_component_is_infinite() {
        assert!(Vector3::new(1.0, f64::INFINITY, 3.0).is_infinite());
    }

    #[test]
    fn is_infinite_is_false_when_all_components_are_finite() {
        assert!(!Vector3::new(1.0, 2.0, 3.0).is_infinite());
    }

    #[test]
    fn is_nan_is_true_when_a_component_is_nan() {
        assert!(Vector3::new(1.0, f64::NAN, 3.0).is_nan());
    }

    #[test]
    fn is_nan_is_false_when_no_component_is_nan() {
        assert!(!Vector3::new(1.0, 2.0, 3.0).is_nan());
    }

    #[test]
    fn array_roundtrip_preserves_the_vector() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(Vector3::from_array(v.to_array()), v);
    }

    #[test]
    fn addition_is_commutative() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, b + a);
    }

    #[test]
    fn the_zero_vector_is_the_additive_identity() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(v + Vector3::ZERO, v);
    }

    #[test]
    fn negation_is_the_additive_inverse() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(v + (-v), Vector3::ZERO);
    }

    #[test]
    fn dot_product_is_symmetric() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a.dot(b), b.dot(a));
    }

    #[test]
    fn cross_product_is_anticommutative() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a.cross(b), -(b.cross(a)));
    }

    #[test]
    fn cross_product_is_orthogonal_to_both_factors() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        let c = a.cross(b);
        assert_eq!(c.dot(a), 0.0);
        assert_eq!(c.dot(b), 0.0);
    }

    #[test]
    fn projection_and_rejection_reconstruct_the_vector() {
        let v = Vector3::new(2.0, 3.0, 0.0);
        let axis = Vector3::new(1.0, 0.0, 0.0);
        assert_eq!(v.project_onto(axis) + v.reject_from(axis), v);
    }

    #[test]
    fn equality_holds_only_when_all_components_match() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(v, Vector3::new(1.0, 2.0, 3.0));
        assert_ne!(v, Vector3::new(9.0, 2.0, 3.0));
        assert_ne!(v, Vector3::new(1.0, 9.0, 3.0));
        assert_ne!(v, Vector3::new(1.0, 2.0, 9.0));
    }

    #[test]
    fn the_operations_are_generic_over_f32() {
        let v = Vector3::new(2.0_f32, 3.0, 6.0);
        assert_eq!(v.norm(), 7.0);
        assert_eq!(v.dot(v), 49.0);
    }
}

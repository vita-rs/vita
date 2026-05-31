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

impl<T, S: Copy> Mul<S> for Vector3<T>
where
    T: Mul<S, Output = T>,
{
    type Output = Self;
    /// Scales every component by `rhs`, preserving the element's unit.
    #[inline]
    fn mul(self, rhs: S) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<T, S: Copy> MulAssign<S> for Vector3<T>
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

impl<T, S: Copy> Div<S> for Vector3<T>
where
    T: Div<S, Output = T>,
{
    type Output = Self;
    /// Divides every component by `rhs`, preserving the element's unit.
    #[inline]
    fn div(self, rhs: S) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl<T, S: Copy> DivAssign<S> for Vector3<T>
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
    pub fn lerp<S: Copy>(self, rhs: Self, t: S) -> Self
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

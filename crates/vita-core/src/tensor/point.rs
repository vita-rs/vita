use core::ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign};

use crate::Scalar;
use crate::tensor::Vector3;

/// A point with three coordinates `x`, `y`, and `z`.
pub struct Point3<T> {
    /// The first coordinate.
    pub x: T,
    /// The second coordinate.
    pub y: T,
    /// The third coordinate.
    pub z: T,
}

impl<T: ::core::marker::Copy> ::core::marker::Copy for Point3<T> {}

impl<T: ::core::clone::Clone> ::core::clone::Clone for Point3<T> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.clone(),
            z: self.z.clone(),
        }
    }
}

impl<T: ::core::fmt::Debug> ::core::fmt::Debug for Point3<T> {
    fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
        f.debug_struct("Point3")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("z", &self.z)
            .finish()
    }
}

impl<T: ::core::cmp::PartialEq> ::core::cmp::PartialEq for Point3<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
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

impl<T: Default> Default for Point3<T> {
    /// Returns the origin.
    #[inline]
    fn default() -> Self {
        Self {
            x: T::default(),
            y: T::default(),
            z: T::default(),
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

impl<T: Add<Output = T> + Sub<Output = T> + Copy> Point3<T> {
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

impl<V: Scalar> Point3<V> {
    /// The origin, `(0, 0, 0)`.
    pub const ORIGIN: Self = Self::new(V::ZERO, V::ZERO, V::ZERO);

    /// Returns the squared Euclidean distance between `self` and `rhs`.
    ///
    /// Cheaper than [`distance`][Self::distance] and sufficient whenever only
    /// relative distances are compared.
    #[inline]
    pub fn distance_squared(self, rhs: Self) -> V {
        (self - rhs).norm_squared()
    }

    /// Returns the Euclidean distance between `self` and `rhs`.
    #[inline]
    pub fn distance(self, rhs: Self) -> V {
        (self - rhs).norm()
    }

    /// Returns the midpoint of the segment between `self` and `rhs`.
    #[inline]
    pub fn midpoint(self, rhs: Self) -> Self {
        self + (rhs - self) * V::from_f64(0.5)
    }

    /// Returns the centroid (arithmetic mean) of `points`, or the
    /// [`ORIGIN`][Self::ORIGIN] when `points` is empty.
    #[inline]
    pub fn centroid(points: &[Self]) -> Self {
        if points.is_empty() {
            return Self::ORIGIN;
        }
        let mut sum = Vector3::<V>::ZERO;
        for p in points {
            sum += p.to_vector();
        }
        Self::from_vector(sum / V::from_f64(points.len() as f64))
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

    /// Restricts every coordinate to the interval `[min, max]`, confining the
    /// point to an axis-aligned bounding box.
    ///
    /// # Panics
    ///
    /// Panics if any coordinate of `min` exceeds the corresponding coordinate
    /// of `max`.
    #[inline]
    pub fn clamp(self, min: Self, max: Self) -> Self {
        Self::new(
            self.x.clamp(min.x, max.x),
            self.y.clamp(min.y, max.y),
            self.z.clamp(min.z, max.z),
        )
    }

    /// Returns the component-wise floor (largest integer not exceeding each
    /// coordinate).
    #[inline]
    pub fn floor(self) -> Self {
        Self::new(self.x.floor(), self.y.floor(), self.z.floor())
    }

    /// Returns the component-wise ceiling (smallest integer not less than each
    /// coordinate).
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

    /// Returns the component-wise fractional part, i.e. the position of
    /// `self` within the unit-grid cell it occupies.
    #[inline]
    pub fn fract(self) -> Self {
        Self::new(self.x.fract(), self.y.fract(), self.z.fract())
    }

    /// Returns `true` if every coordinate is finite.
    #[inline]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Returns `true` if any coordinate is positive or negative infinity.
    #[inline]
    pub fn is_infinite(self) -> bool {
        self.x.is_infinite() || self.y.is_infinite() || self.z.is_infinite()
    }

    /// Returns `true` if any coordinate is `NaN`.
    #[inline]
    pub fn is_nan(self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }
}

use core::fmt;
use core::iter::Sum;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::Vector3;
use crate::{Quantity, Scalar};

/// A 3×3 matrix stored as three column vectors.
///
/// A matrix is a linear map: its columns are the images of the basis vectors,
/// so it applies to a [`Vector3`] and composes with another matrix. It does not
/// apply to a [`Point3`](super::Point3), which would need an origin.
#[derive(Clone, Copy, Default, PartialEq)]
pub struct Matrix3<T> {
    cols: [Vector3<T>; 3],
}

impl<T: fmt::Debug> fmt::Debug for Matrix3<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Matrix3")
            .field("x_col", &self.cols[0])
            .field("y_col", &self.cols[1])
            .field("z_col", &self.cols[2])
            .finish()
    }
}

impl<T> Matrix3<T> {
    /// Constructs a matrix from its three column vectors.
    #[inline]
    pub const fn from_cols(x_col: Vector3<T>, y_col: Vector3<T>, z_col: Vector3<T>) -> Self {
        Self {
            cols: [x_col, y_col, z_col],
        }
    }

    /// Applies `f` to every element in column-major order, returning the
    /// resulting matrix.
    #[inline]
    pub fn map<U, F: FnMut(T) -> U>(self, mut f: F) -> Matrix3<U> {
        let [c0, c1, c2] = self.cols;
        Matrix3::from_cols(c0.map(&mut f), c1.map(&mut f), c2.map(&mut f))
    }

    /// Combines `self` and `rhs` element-wise in column-major order through
    /// `f`.
    #[inline]
    pub fn zip_map<U, R, F: FnMut(T, U) -> R>(self, rhs: Matrix3<U>, mut f: F) -> Matrix3<R> {
        let [c0, c1, c2] = self.cols;
        let [r0, r1, r2] = rhs.cols;
        Matrix3::from_cols(
            c0.zip_map(r0, &mut f),
            c1.zip_map(r1, &mut f),
            c2.zip_map(r2, &mut f),
        )
    }
}

impl<T: Copy> Matrix3<T> {
    /// Constructs a matrix from its three row vectors.
    #[inline]
    pub fn from_rows(x_row: Vector3<T>, y_row: Vector3<T>, z_row: Vector3<T>) -> Self {
        Self::from_cols(
            Vector3::new(x_row.x, y_row.x, z_row.x),
            Vector3::new(x_row.y, y_row.y, z_row.y),
            Vector3::new(x_row.z, y_row.z, z_row.z),
        )
    }

    /// Constructs a matrix from a column-major array of nine elements.
    #[inline]
    pub fn from_cols_array(m: &[T; 9]) -> Self {
        Self::from_cols(
            Vector3::new(m[0], m[1], m[2]),
            Vector3::new(m[3], m[4], m[5]),
            Vector3::new(m[6], m[7], m[8]),
        )
    }

    /// Returns the elements as a column-major array of nine elements.
    #[inline]
    pub fn to_cols_array(self) -> [T; 9] {
        [
            self.cols[0].x,
            self.cols[0].y,
            self.cols[0].z,
            self.cols[1].x,
            self.cols[1].y,
            self.cols[1].z,
            self.cols[2].x,
            self.cols[2].y,
            self.cols[2].z,
        ]
    }

    /// Returns the column at `index`.
    ///
    /// # Panics
    ///
    /// Panics if `index` is greater than `2`.
    #[inline]
    pub fn col(self, index: usize) -> Vector3<T> {
        self.cols[index]
    }

    /// Returns the row at `index`.
    ///
    /// # Panics
    ///
    /// Panics if `index` is greater than `2`.
    #[inline]
    pub fn row(self, index: usize) -> Vector3<T> {
        Vector3::new(
            self.cols[0][index],
            self.cols[1][index],
            self.cols[2][index],
        )
    }

    /// Returns the main diagonal as a vector.
    #[inline]
    pub fn diagonal(self) -> Vector3<T> {
        Vector3::new(self.cols[0].x, self.cols[1].y, self.cols[2].z)
    }

    /// Returns the transpose, exchanging rows and columns.
    #[inline]
    pub fn transpose(self) -> Self {
        Self::from_cols(self.row(0), self.row(1), self.row(2))
    }
}

impl<T: Copy + Default> Matrix3<T> {
    /// Constructs the diagonal matrix carrying `diagonal`, whose off-diagonal
    /// elements are the default (zero) value of `T`.
    #[inline]
    pub fn from_diagonal(diagonal: Vector3<T>) -> Self {
        let zero = T::default();
        Self::from_cols(
            Vector3::new(diagonal.x, zero, zero),
            Vector3::new(zero, diagonal.y, zero),
            Vector3::new(zero, zero, diagonal.z),
        )
    }
}

impl<T: Neg<Output = T>> Neg for Matrix3<T> {
    type Output = Self;
    /// Returns the element-wise negation of `self`.
    #[inline]
    fn neg(self) -> Self {
        let [c0, c1, c2] = self.cols;
        Self::from_cols(-c0, -c1, -c2)
    }
}

impl<T: Add<Output = T>> Add for Matrix3<T> {
    type Output = Self;
    /// Returns the element-wise sum of `self` and `rhs`.
    #[inline]
    fn add(self, rhs: Self) -> Self {
        let [a0, a1, a2] = self.cols;
        let [b0, b1, b2] = rhs.cols;
        Self::from_cols(a0 + b0, a1 + b1, a2 + b2)
    }
}

impl<T: AddAssign> AddAssign for Matrix3<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        let [b0, b1, b2] = rhs.cols;
        self.cols[0] += b0;
        self.cols[1] += b1;
        self.cols[2] += b2;
    }
}

impl<T: Sub<Output = T>> Sub for Matrix3<T> {
    type Output = Self;
    /// Returns the element-wise difference of `self` and `rhs`.
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let [a0, a1, a2] = self.cols;
        let [b0, b1, b2] = rhs.cols;
        Self::from_cols(a0 - b0, a1 - b1, a2 - b2)
    }
}

impl<T: SubAssign> SubAssign for Matrix3<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        let [b0, b1, b2] = rhs.cols;
        self.cols[0] -= b0;
        self.cols[1] -= b1;
        self.cols[2] -= b2;
    }
}

impl<T, S: Quantity> Mul<S> for Matrix3<T>
where
    T: Mul<S>,
{
    type Output = Matrix3<T::Output>;
    /// Scales every element by `rhs`, whose dimension multiplies into the
    /// element's.
    #[inline]
    fn mul(self, rhs: S) -> Self::Output {
        self.map(|element| element * rhs)
    }
}

impl<T, S: Quantity> MulAssign<S> for Matrix3<T>
where
    T: MulAssign<S>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: S) {
        self.cols[0] *= rhs;
        self.cols[1] *= rhs;
        self.cols[2] *= rhs;
    }
}

impl<T> Mul<Matrix3<T>> for f32
where
    f32: Mul<T>,
{
    type Output = Matrix3<<f32 as Mul<T>>::Output>;
    /// Scales every element of `rhs` by `self`.
    #[inline]
    fn mul(self, rhs: Matrix3<T>) -> Self::Output {
        rhs.map(|element| self * element)
    }
}

impl<T> Mul<Matrix3<T>> for f64
where
    f64: Mul<T>,
{
    type Output = Matrix3<<f64 as Mul<T>>::Output>;
    /// Scales every element of `rhs` by `self`.
    #[inline]
    fn mul(self, rhs: Matrix3<T>) -> Self::Output {
        rhs.map(|element| self * element)
    }
}

impl<T, S: Quantity> Div<S> for Matrix3<T>
where
    T: Div<S>,
{
    type Output = Matrix3<T::Output>;
    /// Divides every element by `rhs`, whose dimension divides out of the
    /// element's.
    #[inline]
    fn div(self, rhs: S) -> Self::Output {
        self.map(|element| element / rhs)
    }
}

impl<T, S: Quantity> DivAssign<S> for Matrix3<T>
where
    T: DivAssign<S>,
{
    #[inline]
    fn div_assign(&mut self, rhs: S) {
        self.cols[0] /= rhs;
        self.cols[1] /= rhs;
        self.cols[2] /= rhs;
    }
}

impl<T, U: Quantity> Mul<Vector3<U>> for Matrix3<T>
where
    T: Mul<U> + Copy,
    T::Output: Add<Output = T::Output>,
{
    type Output = Vector3<T::Output>;
    /// Applies the matrix to the vector, returning `self * rhs`.
    #[inline]
    fn mul(self, rhs: Vector3<U>) -> Self::Output {
        self.cols[0] * rhs.x + self.cols[1] * rhs.y + self.cols[2] * rhs.z
    }
}

impl<T, U: Quantity> Mul<Matrix3<U>> for Matrix3<T>
where
    T: Mul<U> + Copy,
    T::Output: Add<Output = T::Output>,
{
    type Output = Matrix3<T::Output>;
    /// Returns the matrix product `self * rhs`, the composition of the two
    /// maps.
    #[inline]
    fn mul(self, rhs: Matrix3<U>) -> Self::Output {
        Matrix3::from_cols(self * rhs.cols[0], self * rhs.cols[1], self * rhs.cols[2])
    }
}

impl<T, U: Quantity> MulAssign<Matrix3<U>> for Matrix3<T>
where
    T: Mul<U, Output = T> + Add<Output = T> + Copy,
{
    #[inline]
    fn mul_assign(&mut self, rhs: Matrix3<U>) {
        *self = *self * rhs;
    }
}

impl<T: Add<Output = T> + Default> Sum for Matrix3<T> {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |acc, m| acc + m)
    }
}

impl<'a, T: Add<Output = T> + Default + Copy> Sum<&'a Matrix3<T>> for Matrix3<T> {
    #[inline]
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.copied().fold(Self::default(), |acc, m| acc + m)
    }
}

impl<Q: Quantity> Matrix3<Q> {
    /// The zero matrix.
    pub const ZERO: Self = Self::from_cols(Vector3::ZERO, Vector3::ZERO, Vector3::ZERO);

    /// Returns the trace, the sum of the diagonal elements.
    #[inline]
    pub fn trace(self) -> Q {
        self.cols[0].x + self.cols[1].y + self.cols[2].z
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

    /// Returns the element-wise absolute value.
    #[inline]
    pub fn abs(self) -> Self {
        self.map(Quantity::abs)
    }

    /// Returns the element-wise minimum of `self` and `other`, ignoring
    /// NaN.
    ///
    /// Where one of the two is NaN, the other is taken.
    #[inline]
    pub fn min(self, other: Self) -> Self {
        self.zip_map(other, Quantity::min)
    }

    /// Returns the element-wise maximum of `self` and `other`, ignoring
    /// NaN.
    ///
    /// Where one of the two is NaN, the other is taken.
    #[inline]
    pub fn max(self, other: Self) -> Self {
        self.zip_map(other, Quantity::max)
    }

    /// Restricts every element to the interval `[lo, hi]`.
    ///
    /// # Panics
    ///
    /// Panics if any element of `lo` exceeds the corresponding element of
    /// `hi`, or if either is NaN.
    #[inline]
    pub fn clamp(self, lo: Self, hi: Self) -> Self {
        let [c0, c1, c2] = self.cols;
        let [l0, l1, l2] = lo.cols;
        let [h0, h1, h2] = hi.cols;
        Self::from_cols(c0.clamp(l0, h0), c1.clamp(l1, h1), c2.clamp(l2, h2))
    }

    /// Returns the element-wise midpoint of `self` and `other`.
    #[inline]
    pub fn midpoint(self, other: Self) -> Self {
        self.zip_map(other, Quantity::midpoint)
    }

    /// Returns the smallest of the nine elements.
    #[inline]
    pub fn min_element(self) -> Q {
        let [c0, c1, c2] = self.cols;
        c0.min_element().min(c1.min_element()).min(c2.min_element())
    }

    /// Returns the largest of the nine elements.
    #[inline]
    pub fn max_element(self) -> Q {
        let [c0, c1, c2] = self.cols;
        c0.max_element().max(c1.max_element()).max(c2.max_element())
    }

    /// Returns the element-wise sign, each `1.0`, `-1.0`, or NaN.
    #[inline]
    pub fn signum(self) -> Matrix3<Q::Value> {
        self.map(Quantity::signum)
    }

    /// Returns a matrix with the magnitudes of `self` and the element-wise
    /// signs of `sign`.
    #[inline]
    pub fn copysign(self, sign: Self) -> Self {
        self.zip_map(sign, Quantity::copysign)
    }

    /// Returns the element-wise floor.
    #[inline]
    pub fn floor(self) -> Self {
        self.map(Quantity::floor)
    }

    /// Returns the element-wise ceiling.
    #[inline]
    pub fn ceil(self) -> Self {
        self.map(Quantity::ceil)
    }

    /// Returns the element-wise nearest integer, rounding halves away from
    /// zero.
    #[inline]
    pub fn round(self) -> Self {
        self.map(Quantity::round)
    }

    /// Returns the element-wise nearest integer, rounding halves to even.
    #[inline]
    pub fn round_ties_even(self) -> Self {
        self.map(Quantity::round_ties_even)
    }

    /// Returns the element-wise truncation toward zero.
    #[inline]
    pub fn trunc(self) -> Self {
        self.map(Quantity::trunc)
    }

    /// Returns the element-wise fractional part.
    #[inline]
    pub fn fract(self) -> Self {
        self.map(Quantity::fract)
    }

    /// Returns the element-wise Euclidean quotient against `rhs`.
    #[inline]
    pub fn div_euclid(self, rhs: Self) -> Matrix3<Q::Value> {
        self.zip_map(rhs, Quantity::div_euclid)
    }

    /// Returns the element-wise least nonnegative remainder against `rhs`.
    #[inline]
    pub fn rem_euclid(self, rhs: Self) -> Self {
        self.zip_map(rhs, Quantity::rem_euclid)
    }

    /// Returns the element-wise fused multiply-add `self * a + b`, each computed with a
    /// single rounding error.
    #[inline]
    pub fn mul_add(self, a: Q::Value, b: Self) -> Self {
        self.zip_map(b, |factor, addend| factor.mul_add(a, addend))
    }

    /// Returns the element-wise hypotenuse of `self` and `other`, each computed without
    /// unnecessary overflow or underflow.
    #[inline]
    pub fn hypot(self, other: Self) -> Self {
        self.zip_map(other, Quantity::hypot)
    }

    /// Returns `true` if any element is NaN.
    #[inline]
    pub fn is_nan(self) -> bool {
        let [c0, c1, c2] = self.cols;
        c0.is_nan() || c1.is_nan() || c2.is_nan()
    }

    /// Returns `true` if any element is positive or negative infinity.
    #[inline]
    pub fn is_infinite(self) -> bool {
        let [c0, c1, c2] = self.cols;
        c0.is_infinite() || c1.is_infinite() || c2.is_infinite()
    }

    /// Returns `true` if every element is finite.
    #[inline]
    pub fn is_finite(self) -> bool {
        let [c0, c1, c2] = self.cols;
        c0.is_finite() && c1.is_finite() && c2.is_finite()
    }

    /// Returns `true` if every element is neither zero, subnormal, infinite,
    /// nor NaN.
    #[inline]
    pub fn is_normal(self) -> bool {
        let [c0, c1, c2] = self.cols;
        c0.is_normal() && c1.is_normal() && c2.is_normal()
    }

    /// Returns `true` if any element is subnormal.
    #[inline]
    pub fn is_subnormal(self) -> bool {
        let [c0, c1, c2] = self.cols;
        c0.is_subnormal() || c1.is_subnormal() || c2.is_subnormal()
    }
}

impl<V: Scalar> Matrix3<V> {
    /// The identity matrix.
    pub const IDENTITY: Self = Self::from_cols(Vector3::X, Vector3::Y, Vector3::Z);

    /// Constructs a right-handed rotation of `angle` radians about the `x`
    /// axis.
    #[inline]
    pub fn from_rotation_x(angle: V) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self::from_cols(
            Vector3::X,
            Vector3::new(V::ZERO, cos, sin),
            Vector3::new(V::ZERO, -sin, cos),
        )
    }

    /// Constructs a right-handed rotation of `angle` radians about the `y`
    /// axis.
    #[inline]
    pub fn from_rotation_y(angle: V) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self::from_cols(
            Vector3::new(cos, V::ZERO, -sin),
            Vector3::Y,
            Vector3::new(sin, V::ZERO, cos),
        )
    }

    /// Constructs a right-handed rotation of `angle` radians about the `z`
    /// axis.
    #[inline]
    pub fn from_rotation_z(angle: V) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self::from_cols(
            Vector3::new(cos, sin, V::ZERO),
            Vector3::new(-sin, cos, V::ZERO),
            Vector3::Z,
        )
    }

    /// Constructs a right-handed rotation of `angle` radians about the unit
    /// vector `axis` (Rodrigues' formula).
    ///
    /// The result is a rotation only when `axis` is normalized.
    #[inline]
    pub fn from_axis_angle(axis: Vector3<V>, angle: V) -> Self {
        let (sin, cos) = angle.sin_cos();
        let complement = V::ONE - cos;
        let (x, y, z) = (axis.x, axis.y, axis.z);
        Self::from_cols(
            Vector3::new(
                cos + x * x * complement,
                y * x * complement + z * sin,
                z * x * complement - y * sin,
            ),
            Vector3::new(
                x * y * complement - z * sin,
                cos + y * y * complement,
                z * y * complement + x * sin,
            ),
            Vector3::new(
                x * z * complement + y * sin,
                y * z * complement - x * sin,
                cos + z * z * complement,
            ),
        )
    }

    /// Constructs the outer (dyadic) product `a ⊗ b`, the 3×3 matrix whose
    /// `(i, j)` entry is `aᵢ · bⱼ`.
    #[inline]
    pub fn outer_product(a: Vector3<V>, b: Vector3<V>) -> Self {
        Self::from_cols(a * b.x, a * b.y, a * b.z)
    }

    /// Returns the determinant.
    #[inline]
    pub fn determinant(self) -> V {
        self.cols[0].dot(self.cols[1].cross(self.cols[2]))
    }

    /// Returns `true` if the matrix has a finite, non-zero determinant and is
    /// therefore invertible.
    #[inline]
    pub fn is_invertible(self) -> bool {
        let determinant = self.determinant();
        determinant != V::ZERO && determinant.is_finite()
    }

    /// Returns the inverse, or `None` if the matrix is not invertible.
    ///
    /// Returns `None` on the same boundary as [`is_invertible`](Self::is_invertible),
    /// including the matrices whose inverse is representable but whose determinant is not.
    #[inline]
    pub fn try_inverse(self) -> Option<Self> {
        let [c0, c1, c2] = self.cols;
        let (r0, r1, r2) = (c1.cross(c2), c2.cross(c0), c0.cross(c1));
        let determinant = c0.dot(r0);
        if determinant == V::ZERO || !determinant.is_finite() {
            return None;
        }
        Some(Self::from_rows(r0, r1, r2) / determinant)
    }

    /// Returns the inverse.
    ///
    /// # Panics
    ///
    /// Panics if the matrix is not invertible.
    #[inline]
    pub fn inverse(self) -> Self {
        self.try_inverse().expect("matrix is not invertible")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use core::f64::consts::FRAC_PI_2;

    use crate::units::length::{Angstrom, Length};

    const TOL: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL
    }

    fn vectors_close(a: Vector3<f64>, b: Vector3<f64>) -> bool {
        close(a.x, b.x) && close(a.y, b.y) && close(a.z, b.z)
    }

    fn matrices_close(a: Matrix3<f64>, b: Matrix3<f64>) -> bool {
        vectors_close(a.col(0), b.col(0))
            && vectors_close(a.col(1), b.col(1))
            && vectors_close(a.col(2), b.col(2))
    }

    fn length(value: f64) -> Length<f64, Angstrom> {
        Length::new(value)
    }

    fn uniform(value: f64) -> Matrix3<f64> {
        Matrix3::from_cols(
            Vector3::splat(value),
            Vector3::splat(value),
            Vector3::splat(value),
        )
    }

    fn matrix() -> Matrix3<f64> {
        Matrix3::from_cols(
            Vector3::new(1.0, 4.0, 7.0),
            Vector3::new(2.0, 5.0, 8.0),
            Vector3::new(3.0, 6.0, 10.0),
        )
    }

    fn singular() -> Matrix3<f64> {
        Matrix3::from_diagonal(Vector3::new(1.0, 1.0, 0.0))
    }

    fn length_matrix() -> Matrix3<Length<f64, Angstrom>> {
        Matrix3::from_diagonal(Vector3::new(length(2.0), length(3.0), length(4.0)))
    }

    #[test]
    fn default_is_the_zero_matrix() {
        assert_eq!(Matrix3::<f64>::default(), Matrix3::ZERO);
    }

    #[test]
    fn the_zero_matrix_has_zero_elements() {
        assert_eq!(Matrix3::<f64>::ZERO.to_cols_array(), [0.0; 9]);
    }

    #[test]
    fn trace_of_the_zero_matrix_is_zero() {
        assert_eq!(Matrix3::<f64>::ZERO.trace(), 0.0);
    }

    #[test]
    fn determinant_of_the_zero_matrix_is_zero() {
        assert_eq!(Matrix3::<f64>::ZERO.determinant(), 0.0);
    }

    #[test]
    fn try_inverse_of_the_zero_matrix_is_none() {
        assert!(Matrix3::<f64>::ZERO.try_inverse().is_none());
    }

    #[test]
    fn summing_no_matrices_yields_the_zero_matrix() {
        let none: [Matrix3<f64>; 0] = [];
        assert_eq!(none.into_iter().sum::<Matrix3<f64>>(), Matrix3::ZERO);
    }

    #[test]
    fn lerp_at_zero_yields_the_starting_matrix() {
        assert_eq!(matrix().lerp(Matrix3::ZERO, 0.0), matrix());
    }

    #[test]
    fn lerp_at_one_yields_the_ending_matrix() {
        assert_eq!(matrix().lerp(Matrix3::ZERO, 1.0), Matrix3::ZERO);
    }

    #[test]
    fn from_rotation_x_with_zero_angle_is_the_identity() {
        assert_eq!(Matrix3::from_rotation_x(0.0), Matrix3::<f64>::IDENTITY);
    }

    #[test]
    fn from_axis_angle_with_zero_angle_is_the_identity() {
        let rotation = Matrix3::from_axis_angle(Vector3::<f64>::Z, 0.0);
        assert!(matrices_close(rotation, Matrix3::IDENTITY));
    }

    #[test]
    fn debug_labels_the_three_columns() {
        assert_eq!(
            format!("{:?}", matrix()),
            "Matrix3 { \
             x_col: Vector3 { x: 1.0, y: 4.0, z: 7.0 }, \
             y_col: Vector3 { x: 2.0, y: 5.0, z: 8.0 }, \
             z_col: Vector3 { x: 3.0, y: 6.0, z: 10.0 } }"
        );
    }

    #[test]
    fn from_cols_takes_the_three_columns_in_order() {
        assert_eq!(matrix().col(0), Vector3::new(1.0, 4.0, 7.0));
    }

    #[test]
    fn map_applies_the_function_to_every_element() {
        assert_eq!(matrix().map(|e| e as i32).col(0), Vector3::new(1, 4, 7));
    }

    #[test]
    fn zip_map_combines_the_two_matrices_element_wise() {
        let combined = matrix().zip_map(matrix(), |a, b| (a + b) as i32);
        assert_eq!(combined.col(0), Vector3::new(2, 8, 14));
    }

    #[test]
    fn from_rows_takes_the_three_rows_in_order() {
        let rows = Matrix3::from_rows(
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(4.0, 5.0, 6.0),
            Vector3::new(7.0, 8.0, 10.0),
        );
        assert_eq!(rows, matrix());
    }

    #[test]
    fn from_cols_array_reads_the_elements_in_column_major_order() {
        let elements = [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0];
        assert_eq!(Matrix3::from_cols_array(&elements), matrix());
    }

    #[test]
    fn to_cols_array_writes_the_elements_in_column_major_order() {
        let elements = [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0];
        assert_eq!(matrix().to_cols_array(), elements);
    }

    #[test]
    fn col_yields_the_column_at_that_index() {
        assert_eq!(matrix().col(1), Vector3::new(2.0, 5.0, 8.0));
    }

    #[test]
    fn row_yields_the_row_at_that_index() {
        assert_eq!(matrix().row(1), Vector3::new(4.0, 5.0, 6.0));
    }

    #[test]
    fn diagonal_yields_the_main_diagonal() {
        assert_eq!(matrix().diagonal(), Vector3::new(1.0, 5.0, 10.0));
    }

    #[test]
    fn transpose_exchanges_rows_and_columns() {
        assert_eq!(matrix().transpose().col(0), matrix().row(0));
    }

    #[test]
    fn from_diagonal_places_the_vector_on_the_diagonal() {
        let diagonal = Matrix3::from_diagonal(Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(
            diagonal.to_cols_array(),
            [1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]
        );
    }

    #[test]
    fn negation_flips_the_sign_of_every_element() {
        assert_eq!((-matrix()).col(0), Vector3::new(-1.0, -4.0, -7.0));
    }

    #[test]
    fn addition_sums_the_elements() {
        assert_eq!((matrix() + matrix()).col(0), Vector3::new(2.0, 8.0, 14.0));
    }

    #[test]
    fn add_assign_sums_the_elements_in_place() {
        let mut m = matrix();
        m += matrix();
        assert_eq!(m.col(0), Vector3::new(2.0, 8.0, 14.0));
    }

    #[test]
    fn subtraction_differences_the_elements() {
        assert_eq!(matrix() - matrix(), Matrix3::ZERO);
    }

    #[test]
    fn sub_assign_differences_the_elements_in_place() {
        let mut m = matrix();
        m -= matrix();
        assert_eq!(m, Matrix3::ZERO);
    }

    #[test]
    fn multiplying_by_a_number_scales_every_element() {
        assert_eq!((matrix() * 2.0).col(0), Vector3::new(2.0, 8.0, 14.0));
    }

    #[test]
    fn mul_assign_scales_every_element_in_place() {
        let mut m = matrix();
        m *= 2.0;
        assert_eq!(m.col(0), Vector3::new(2.0, 8.0, 14.0));
    }

    #[test]
    fn a_double_on_the_left_scales_every_element() {
        assert_eq!((2.0 * matrix()).col(0), Vector3::new(2.0, 8.0, 14.0));
    }

    #[test]
    fn a_single_on_the_left_scales_every_element() {
        let m = Matrix3::from_diagonal(Vector3::new(1.0_f32, 2.0, 3.0));
        assert_eq!((2.0_f32 * m).diagonal(), Vector3::new(2.0_f32, 4.0, 6.0));
    }

    #[test]
    fn dividing_by_a_number_scales_every_element_down() {
        assert_eq!(((matrix() * 2.0) / 2.0), matrix());
    }

    #[test]
    fn div_assign_scales_every_element_down_in_place() {
        let mut m = matrix() * 2.0;
        m /= 2.0;
        assert_eq!(m, matrix());
    }

    #[test]
    fn multiplying_by_a_vector_combines_the_columns() {
        assert_eq!(matrix() * Vector3::new(1.0, 0.0, 0.0), matrix().col(0));
    }

    #[test]
    fn multiplying_by_a_matrix_composes_the_maps() {
        let a = Matrix3::from_diagonal(Vector3::new(2.0, 3.0, 4.0));
        let b = Matrix3::from_diagonal(Vector3::new(5.0, 6.0, 7.0));
        assert_eq!((a * b).diagonal(), Vector3::new(10.0, 18.0, 28.0));
    }

    #[test]
    fn mul_assign_by_a_matrix_composes_in_place() {
        let mut a = Matrix3::from_diagonal(Vector3::new(2.0, 3.0, 4.0));
        a *= Matrix3::from_diagonal(Vector3::new(5.0, 6.0, 7.0));
        assert_eq!(a.diagonal(), Vector3::new(10.0, 18.0, 28.0));
    }

    #[test]
    fn summing_owned_matrices_adds_them() {
        let sum = [matrix(), matrix()].into_iter().sum::<Matrix3<f64>>();
        assert_eq!(sum, matrix() * 2.0);
    }

    #[test]
    fn summing_borrowed_matrices_adds_them() {
        let sum = [matrix(), matrix()].iter().sum::<Matrix3<f64>>();
        assert_eq!(sum, matrix() * 2.0);
    }

    #[test]
    fn trace_adds_the_diagonal_elements() {
        assert_eq!(matrix().trace(), 16.0);
    }

    #[test]
    fn lerp_interpolates_between_the_matrices() {
        assert_eq!(Matrix3::ZERO.lerp(matrix() * 4.0, 0.25), matrix());
    }

    #[test]
    fn abs_takes_the_magnitude_of_every_element() {
        assert_eq!((-matrix()).abs(), matrix());
    }

    #[test]
    fn min_takes_the_smaller_of_each_pair_of_elements() {
        assert_eq!(matrix().min(Matrix3::ZERO), Matrix3::ZERO);
    }

    #[test]
    fn max_takes_the_larger_of_each_pair_of_elements() {
        assert_eq!(matrix().max(Matrix3::ZERO), matrix());
    }

    #[test]
    fn clamp_raises_an_element_below_the_interval_to_the_lower_bound() {
        let clamped = (-matrix()).clamp(Matrix3::ZERO, uniform(99.0));
        assert_eq!(clamped, Matrix3::ZERO);
    }

    #[test]
    fn clamp_leaves_an_element_inside_the_interval_unchanged() {
        assert_eq!(matrix().clamp(Matrix3::ZERO, uniform(99.0)), matrix());
    }

    #[test]
    fn clamp_lowers_an_element_above_the_interval_to_the_upper_bound() {
        let clamped = (matrix() * 100.0).clamp(Matrix3::ZERO, uniform(99.0));
        assert_eq!(clamped, uniform(99.0));
    }

    #[test]
    fn midpoint_is_halfway_between_the_matrices() {
        assert_eq!(Matrix3::ZERO.midpoint(matrix() * 2.0), matrix());
    }

    #[test]
    fn min_element_is_the_smallest_element() {
        assert_eq!(matrix().min_element(), 1.0);
    }

    #[test]
    fn max_element_is_the_largest_element() {
        assert_eq!(matrix().max_element(), 10.0);
    }

    #[test]
    fn signum_takes_the_sign_of_every_element() {
        assert_eq!((-matrix()).signum().col(0), Vector3::new(-1.0, -1.0, -1.0));
    }

    #[test]
    fn copysign_keeps_the_magnitudes_and_takes_the_signs_of_its_argument() {
        assert_eq!(matrix().copysign(-matrix()), -matrix());
    }

    #[test]
    fn floor_rounds_every_element_toward_negative_infinity() {
        assert_eq!((matrix() / 2.0).floor().col(0), Vector3::new(0.0, 2.0, 3.0));
    }

    #[test]
    fn ceil_rounds_every_element_toward_positive_infinity() {
        assert_eq!((matrix() / 2.0).ceil().col(0), Vector3::new(1.0, 2.0, 4.0));
    }

    #[test]
    fn round_sends_a_half_away_from_zero() {
        let halves = Matrix3::from_diagonal(Vector3::new(2.5, 3.5, -2.5));
        assert_eq!(halves.round().diagonal(), Vector3::new(3.0, 4.0, -3.0));
    }

    #[test]
    fn round_ties_even_sends_a_half_to_the_even_integer() {
        let halves = Matrix3::from_diagonal(Vector3::new(2.5, 3.5, -2.5));
        assert_eq!(
            halves.round_ties_even().diagonal(),
            Vector3::new(2.0, 4.0, -2.0)
        );
    }

    #[test]
    fn trunc_drops_the_fractional_part_of_every_element() {
        assert_eq!((matrix() / 2.0).trunc().col(0), Vector3::new(0.0, 2.0, 3.0));
    }

    #[test]
    fn fract_keeps_the_fractional_part_of_every_element() {
        assert_eq!((matrix() / 2.0).fract().col(0), Vector3::new(0.5, 0.0, 0.5));
    }

    #[test]
    fn div_euclid_is_the_euclidean_quotient_of_every_element() {
        let divisor = Matrix3::from_diagonal(Vector3::new(3.0, 3.0, 3.0));
        let dividend = Matrix3::from_diagonal(Vector3::new(7.0, -7.0, 0.0));
        assert_eq!(
            dividend.div_euclid(divisor).diagonal(),
            Vector3::new(2.0, -3.0, 0.0)
        );
    }

    #[test]
    fn rem_euclid_is_nonnegative_for_a_negative_element() {
        let divisor = Matrix3::from_diagonal(Vector3::new(3.0, 3.0, 3.0));
        let dividend = Matrix3::from_diagonal(Vector3::new(-1.0, 7.0, 5.0));
        assert_eq!(
            dividend.rem_euclid(divisor).diagonal(),
            Vector3::new(2.0, 1.0, 2.0)
        );
    }

    #[test]
    fn mul_add_scales_then_offsets_every_element() {
        let offset = Matrix3::from_diagonal(Vector3::splat(10.0));
        let expected = matrix() * 2.0 + offset;
        assert_eq!(matrix().mul_add(2.0, offset), expected);
    }

    #[test]
    fn hypot_combines_the_elements_pairwise() {
        let legs = Matrix3::from_diagonal(Vector3::new(4.0, 0.0, 5.0));
        let m = Matrix3::from_diagonal(Vector3::new(3.0, 4.0, 12.0));
        let expected = Matrix3::from_diagonal(Vector3::new(5.0, 4.0, 13.0));
        assert!(matrices_close(m.hypot(legs), expected));
    }

    #[test]
    fn is_nan_holds_when_an_element_is_not_a_number() {
        assert!(Matrix3::from_diagonal(Vector3::new(1.0, f64::NAN, 3.0)).is_nan());
    }

    #[test]
    fn is_infinite_holds_when_an_element_is_infinite() {
        assert!(Matrix3::from_diagonal(Vector3::new(1.0, f64::INFINITY, 3.0)).is_infinite());
    }

    #[test]
    fn is_finite_holds_when_every_element_is_finite() {
        assert!(matrix().is_finite());
    }

    #[test]
    fn is_normal_holds_when_every_element_is_normal() {
        assert!(matrix().is_normal());
    }

    #[test]
    fn is_subnormal_holds_when_an_element_is_subnormal() {
        let tiny = Matrix3::from_diagonal(Vector3::splat(f64::MIN_POSITIVE / 2.0));
        assert!(tiny.is_subnormal());
    }

    #[test]
    fn the_identity_is_one_on_the_diagonal_and_zero_elsewhere() {
        let elements = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert_eq!(Matrix3::<f64>::IDENTITY.to_cols_array(), elements);
    }

    #[test]
    fn from_rotation_x_turns_the_y_axis_toward_the_z_axis() {
        let rotation = Matrix3::from_rotation_x(FRAC_PI_2);
        assert!(vectors_close(rotation * Vector3::<f64>::Y, Vector3::Z));
    }

    #[test]
    fn from_rotation_y_turns_the_z_axis_toward_the_x_axis() {
        let rotation = Matrix3::from_rotation_y(FRAC_PI_2);
        assert!(vectors_close(rotation * Vector3::<f64>::Z, Vector3::X));
    }

    #[test]
    fn from_rotation_z_turns_the_x_axis_toward_the_y_axis() {
        let rotation = Matrix3::from_rotation_z(FRAC_PI_2);
        assert!(vectors_close(rotation * Vector3::<f64>::X, Vector3::Y));
    }

    #[test]
    fn from_axis_angle_about_an_axis_matches_the_axis_rotation() {
        let rotation = Matrix3::from_axis_angle(Vector3::<f64>::Z, FRAC_PI_2);
        assert!(matrices_close(
            rotation,
            Matrix3::from_rotation_z(FRAC_PI_2)
        ));
    }

    #[test]
    fn outer_product_multiplies_every_pair_of_components() {
        let product =
            Matrix3::outer_product(Vector3::new(1.0, 2.0, 3.0), Vector3::new(4.0, 5.0, 6.0));
        assert_eq!(product.col(0), Vector3::new(4.0, 8.0, 12.0));
    }

    #[test]
    fn determinant_is_the_signed_volume_of_the_columns() {
        assert_eq!(matrix().determinant(), -3.0);
    }

    #[test]
    fn is_invertible_holds_for_a_nonsingular_matrix() {
        assert!(matrix().is_invertible());
    }

    #[test]
    fn try_inverse_of_a_nonsingular_matrix_is_some() {
        assert!(matrix().try_inverse().is_some());
    }

    #[test]
    fn inverse_of_a_diagonal_matrix_inverts_each_element() {
        let diagonal = Matrix3::from_diagonal(Vector3::new(2.0, 4.0, 8.0));
        assert!(vectors_close(
            diagonal.inverse().diagonal(),
            Vector3::new(0.5, 0.25, 0.125)
        ));
    }

    #[test]
    fn signum_of_a_length_matrix_is_dimensionless() {
        assert_eq!(
            length_matrix().signum().diagonal(),
            Vector3::new(1.0, 1.0, 1.0)
        );
    }

    #[test]
    fn div_euclid_of_length_matrices_is_dimensionless() {
        let divisor = Matrix3::from_diagonal(Vector3::splat(length(2.0)));
        assert_eq!(
            length_matrix().div_euclid(divisor).diagonal(),
            Vector3::new(1.0, 1.0, 2.0)
        );
    }

    #[test]
    fn mul_add_of_length_matrices_takes_a_dimensionless_factor() {
        let offset = Matrix3::from_diagonal(Vector3::splat(length(10.0)));
        let expected = length_matrix() * 2.0 + offset;
        assert_eq!(length_matrix().mul_add(2.0, offset), expected);
    }

    #[test]
    fn hypot_of_length_matrices_stays_a_length() {
        let legs = Matrix3::from_diagonal(Vector3::new(length(4.0), length(0.0), length(5.0)));
        let m = Matrix3::from_diagonal(Vector3::new(length(3.0), length(4.0), length(12.0)));
        let expected = Matrix3::from_diagonal(Vector3::new(5.0, 4.0, 13.0));
        assert!(matrices_close(m.hypot(legs).map(Quantity::value), expected));
    }

    #[test]
    fn lerp_of_length_matrices_takes_a_dimensionless_factor() {
        let halved = Matrix3::from_diagonal(Vector3::new(length(1.0), length(1.5), length(2.0)));
        assert_eq!(Matrix3::ZERO.lerp(length_matrix(), 0.5), halved);
    }

    #[test]
    fn multiplying_a_dimensionless_matrix_by_a_length_carries_the_dimension() {
        let scaled = Matrix3::from_diagonal(Vector3::new(2.0, 3.0, 4.0)) * length(1.0);
        assert_eq!(scaled, length_matrix());
    }

    #[test]
    fn a_double_on_the_left_scales_a_length_matrix() {
        let halved = Matrix3::from_diagonal(Vector3::new(length(1.0), length(1.5), length(2.0)));
        assert_eq!(2.0 * halved, length_matrix());
    }

    #[test]
    fn dividing_a_length_matrix_by_a_length_leaves_it_dimensionless() {
        let plain = Matrix3::from_diagonal(Vector3::new(2.0, 3.0, 4.0));
        assert_eq!(length_matrix() / length(1.0), plain);
    }

    #[test]
    fn a_length_matrix_maps_a_dimensionless_vector_to_a_length_vector() {
        let expected = Vector3::new(length(2.0), length(3.0), length(4.0));
        assert_eq!(length_matrix() * Vector3::new(1.0, 1.0, 1.0), expected);
    }

    #[test]
    fn a_dimensionless_matrix_maps_a_length_vector_to_a_length_vector() {
        let plain = Matrix3::from_diagonal(Vector3::new(2.0, 3.0, 4.0));
        let expected = Vector3::new(length(2.0), length(3.0), length(4.0));
        assert_eq!(plain * Vector3::splat(length(1.0)), expected);
    }

    #[test]
    fn multiplying_a_length_matrix_by_a_dimensionless_matrix_keeps_the_length() {
        assert_eq!(length_matrix() * Matrix3::<f64>::IDENTITY, length_matrix());
    }

    #[test]
    fn multiplying_a_dimensionless_matrix_by_a_length_matrix_keeps_the_length() {
        assert_eq!(Matrix3::<f64>::IDENTITY * length_matrix(), length_matrix());
    }

    #[test]
    fn mul_assign_by_a_dimensionless_matrix_keeps_the_length() {
        let mut m = length_matrix();
        m *= Matrix3::<f64>::IDENTITY;
        assert_eq!(m, length_matrix());
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn col_out_of_bounds_panics() {
        let _ = matrix().col(3);
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn row_out_of_bounds_panics() {
        let _ = matrix().row(3);
    }

    #[test]
    #[should_panic(expected = "min > max")]
    fn clamp_with_an_inverted_interval_panics() {
        let _ = matrix().clamp(uniform(99.0), Matrix3::ZERO);
    }

    #[test]
    #[should_panic(expected = "matrix is not invertible")]
    fn inverse_of_a_singular_matrix_panics() {
        let _ = singular().inverse();
    }

    #[test]
    fn try_inverse_of_a_singular_matrix_is_none() {
        assert!(singular().try_inverse().is_none());
    }

    #[test]
    fn is_invertible_does_not_hold_for_a_singular_matrix() {
        assert!(!singular().is_invertible());
    }

    #[test]
    fn min_ignores_a_not_a_number_element() {
        let nan = Matrix3::from_diagonal(Vector3::splat(f64::NAN));
        assert_eq!(nan.min(Matrix3::ZERO), Matrix3::ZERO);
    }

    #[test]
    fn max_ignores_a_not_a_number_element() {
        let nan = Matrix3::from_diagonal(Vector3::splat(f64::NAN));
        assert_eq!(nan.max(Matrix3::ZERO), Matrix3::ZERO);
    }

    #[test]
    fn is_nan_does_not_hold_when_no_element_is_not_a_number() {
        assert!(!matrix().is_nan());
    }

    #[test]
    fn is_infinite_does_not_hold_when_every_element_is_finite() {
        assert!(!matrix().is_infinite());
    }

    #[test]
    fn is_finite_does_not_hold_when_an_element_is_infinite() {
        assert!(!Matrix3::from_diagonal(Vector3::splat(f64::INFINITY)).is_finite());
    }

    #[test]
    fn is_normal_does_not_hold_when_an_element_is_zero() {
        assert!(!Matrix3::<f64>::IDENTITY.is_normal());
    }

    #[test]
    fn is_subnormal_does_not_hold_when_every_element_is_normal() {
        assert!(!matrix().is_subnormal());
    }

    #[test]
    fn floor_leaves_an_integer_element_unchanged() {
        assert_eq!(matrix().floor(), matrix());
    }

    #[test]
    fn ceil_leaves_an_integer_element_unchanged() {
        assert_eq!(matrix().ceil(), matrix());
    }

    #[test]
    fn fract_of_an_integer_element_is_zero() {
        assert_eq!(matrix().fract(), Matrix3::ZERO);
    }

    #[test]
    fn signum_of_a_negative_zero_element_is_minus_one() {
        let negative_zero = Matrix3::from_diagonal(Vector3::splat(-0.0));
        assert_eq!(negative_zero.signum().diagonal().x, -1.0);
    }

    #[test]
    fn lerp_beyond_one_extrapolates_past_the_ending_matrix() {
        assert_eq!(Matrix3::ZERO.lerp(matrix(), 2.0), matrix() * 2.0);
    }

    #[test]
    fn lerp_holds_when_the_difference_between_the_ends_overflows() {
        let start = Matrix3::from_diagonal(Vector3::splat(-f64::MAX));
        let end = Matrix3::from_diagonal(Vector3::splat(f64::MAX));
        assert_eq!(start.lerp(end, 0.5), Matrix3::ZERO);
    }

    #[test]
    fn array_roundtrip_preserves_the_matrix() {
        assert_eq!(
            Matrix3::from_cols_array(&matrix().to_cols_array()),
            matrix()
        );
    }

    #[test]
    fn diagonal_roundtrip_preserves_the_diagonal() {
        let diagonal = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(Matrix3::from_diagonal(diagonal).diagonal(), diagonal);
    }

    #[test]
    fn transpose_is_its_own_inverse() {
        assert_eq!(matrix().transpose().transpose(), matrix());
    }

    #[test]
    fn from_rows_is_the_transpose_of_from_cols() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        let c = Vector3::new(7.0, 8.0, 10.0);
        assert_eq!(
            Matrix3::from_rows(a, b, c),
            Matrix3::from_cols(a, b, c).transpose()
        );
    }

    #[test]
    fn the_identity_is_the_multiplicative_identity() {
        assert_eq!(matrix() * Matrix3::<f64>::IDENTITY, matrix());
    }

    #[test]
    fn the_identity_maps_a_vector_to_itself() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(Matrix3::<f64>::IDENTITY * v, v);
    }

    #[test]
    fn a_matrix_times_its_inverse_is_the_identity() {
        assert!(matrices_close(
            matrix() * matrix().inverse(),
            Matrix3::IDENTITY
        ));
    }

    #[test]
    fn the_determinant_is_multiplicative() {
        let other = Matrix3::from_diagonal(Vector3::new(2.0, 3.0, 4.0));
        let product = (matrix() * other).determinant();
        assert!(close(product, matrix().determinant() * other.determinant()));
    }

    #[test]
    fn the_determinant_is_unchanged_by_transposition() {
        assert!(close(
            matrix().transpose().determinant(),
            matrix().determinant()
        ));
    }

    #[test]
    fn the_trace_of_an_outer_product_is_the_dot_product() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(Matrix3::outer_product(a, b).trace(), a.dot(b));
    }

    #[test]
    fn rotating_twice_about_an_axis_adds_the_angles() {
        let twice = Matrix3::from_rotation_z(FRAC_PI_2) * Matrix3::from_rotation_z(FRAC_PI_2);
        assert!(matrices_close(
            twice,
            Matrix3::from_rotation_z(FRAC_PI_2 * 2.0)
        ));
    }

    #[test]
    fn scaling_on_either_side_gives_the_same_matrix() {
        assert_eq!(2.0 * matrix(), matrix() * 2.0);
    }

    #[test]
    fn truncating_and_the_fractional_part_reconstruct_the_matrix() {
        let m = matrix() / 2.0;
        assert_eq!(m.trunc() + m.fract(), m);
    }

    #[test]
    fn min_of_a_matrix_with_itself_is_that_matrix() {
        assert_eq!(matrix().min(matrix()), matrix());
    }

    #[test]
    fn max_of_a_matrix_with_itself_is_that_matrix() {
        assert_eq!(matrix().max(matrix()), matrix());
    }

    #[test]
    fn clamp_agrees_with_taking_the_maximum_then_the_minimum() {
        let (lo, hi) = (uniform(1.5), uniform(2.5));
        assert_eq!(matrix().clamp(lo, hi), matrix().max(lo).min(hi));
    }

    #[test]
    fn midpoint_agrees_with_interpolating_halfway() {
        assert_eq!(
            matrix().midpoint(Matrix3::ZERO),
            matrix().lerp(Matrix3::ZERO, 0.5)
        );
    }

    #[test]
    fn matrices_are_equal_exactly_when_all_columns_match() {
        assert_eq!(matrix(), matrix());
        assert_ne!(matrix(), matrix().transpose());
    }
}

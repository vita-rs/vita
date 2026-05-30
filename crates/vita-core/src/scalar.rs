use core::fmt;
use core::iter::{Product, Sum};
use core::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

/// Numeric scalar type for values in physical-quantity and tensor types.
///
/// Implemented for [`f32`] and [`f64`].
pub trait Scalar:
    Copy
    + PartialEq
    + PartialOrd
    + Default
    + fmt::Debug
    + fmt::Display
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Div<Output = Self>
    + DivAssign
    + Neg<Output = Self>
    + Rem<Output = Self>
    + RemAssign
    + Sum
    + for<'a> Sum<&'a Self>
    + Product
    + for<'a> Product<&'a Self>
{
    /// The additive identity, `0.0`.
    const ZERO: Self;

    /// The multiplicative identity, `1.0`.
    const ONE: Self;

    /// Machine epsilon — the difference between `1` and the next larger representable value.
    const EPSILON: Self;

    /// Returns the absolute value.
    fn abs(self) -> Self;

    /// Returns the minimum of `self` and `other`, ignoring NaN.
    ///
    /// If one argument is NaN, then the other argument is returned.
    fn min(self, other: Self) -> Self;

    /// Returns the maximum of `self` and `other`, ignoring NaN.
    ///
    /// If one argument is NaN, then the other argument is returned.
    fn max(self, other: Self) -> Self;

    /// Restricts `self` to the interval `[lo, hi]`.
    ///
    /// # Panics
    ///
    /// Panics if `lo > hi`, `lo` is NaN, or `hi` is NaN.
    fn clamp(self, lo: Self, hi: Self) -> Self;

    /// Returns `1.0`, `-1.0`, or NaN based on the sign of `self`.
    fn signum(self) -> Self;

    /// Returns a number composed of the magnitude of `self` and the sign of `sign`.
    fn copysign(self, sign: Self) -> Self;

    /// Returns the largest integer less than or equal to `self`.
    fn floor(self) -> Self;

    /// Returns the smallest integer greater than or equal to `self`.
    fn ceil(self) -> Self;

    /// Returns the nearest integer to `self`, with halves rounded away from zero.
    fn round(self) -> Self;

    /// Returns the nearest integer to `self`, with halves rounded to even.
    fn round_ties_even(self) -> Self;

    /// Returns the integer part of `self`, truncated toward zero.
    fn trunc(self) -> Self;

    /// Returns the fractional part of `self`.
    fn fract(self) -> Self;

    /// Euclidean division, the matching counterpart of [`Self::rem_euclid`].
    fn div_euclid(self, rhs: Self) -> Self;

    /// Least nonnegative remainder of `self` divided by `rhs`.
    fn rem_euclid(self, rhs: Self) -> Self;

    /// Returns the reciprocal `1 / self`.
    fn recip(self) -> Self;

    /// Computes `(self * a) + b` with a single rounding error (fused multiply-add).
    fn mul_add(self, a: Self, b: Self) -> Self;

    /// Returns the square root of `self`.
    ///
    /// Returns NaN if `self` is negative (other than `-0.0`).
    fn sqrt(self) -> Self;

    /// Returns the cube root of `self`.
    fn cbrt(self) -> Self;

    /// Raises `self` to an integer power.
    fn powi(self, n: i32) -> Self;

    /// Raises `self` to a floating-point power.
    fn powf(self, n: Self) -> Self;

    /// Returns `sqrt(self * self + other * other)` without unnecessary overflow or underflow.
    fn hypot(self, other: Self) -> Self;

    /// Returns `e^self`.
    fn exp(self) -> Self;

    /// Returns `2^self`.
    fn exp2(self) -> Self;

    /// Returns `e^self - 1`, accurate when `self` is near zero.
    fn exp_m1(self) -> Self;

    /// Returns the natural logarithm of `self`.
    ///
    /// Returns NaN when `self` is negative and `-∞` when `self` is zero.
    fn ln(self) -> Self;

    /// Returns `ln(1 + self)`, accurate when `self` is near zero.
    fn ln_1p(self) -> Self;

    /// Returns the logarithm of `self` with respect to `base`.
    fn log(self, base: Self) -> Self;

    /// Returns the base-2 logarithm of `self`.
    fn log2(self) -> Self;

    /// Returns the base-10 logarithm of `self`.
    fn log10(self) -> Self;

    /// Computes the sine of `self` (in radians).
    fn sin(self) -> Self;

    /// Computes the cosine of `self` (in radians).
    fn cos(self) -> Self;

    /// Computes the tangent of `self` (in radians).
    fn tan(self) -> Self;

    /// Simultaneously computes the sine and cosine of `self`, returning `(sin, cos)`.
    fn sin_cos(self) -> (Self, Self);

    /// Computes the arcsine of `self`, in radians within `[-π/2, π/2]`.
    ///
    /// Returns NaN if `self` is outside `[-1, 1]`.
    fn asin(self) -> Self;

    /// Computes the arccosine of `self`, in radians within `[0, π]`.
    ///
    /// Returns NaN if `self` is outside `[-1, 1]`.
    fn acos(self) -> Self;

    /// Computes the arctangent of `self`, in radians within `[-π/2, π/2]`.
    fn atan(self) -> Self;

    /// Computes the four-quadrant arctangent of `self` (`y`) and `other` (`x`).
    ///
    /// Returns the angle in radians within `[-π, π]`.
    fn atan2(self, other: Self) -> Self;

    /// Hyperbolic sine of `self`.
    fn sinh(self) -> Self;

    /// Hyperbolic cosine of `self`.
    fn cosh(self) -> Self;

    /// Hyperbolic tangent of `self`.
    fn tanh(self) -> Self;

    /// Inverse hyperbolic sine of `self`.
    fn asinh(self) -> Self;

    /// Inverse hyperbolic cosine of `self`.
    fn acosh(self) -> Self;

    /// Inverse hyperbolic tangent of `self`.
    fn atanh(self) -> Self;

    /// Converts radians to degrees.
    fn to_degrees(self) -> Self;

    /// Converts degrees to radians.
    fn to_radians(self) -> Self;

    /// Returns `true` if `self` is NaN.
    fn is_nan(self) -> bool;

    /// Returns `true` if `self` is positive or negative infinity.
    fn is_infinite(self) -> bool;

    /// Returns `true` if `self` is neither NaN nor infinite.
    fn is_finite(self) -> bool;

    /// Returns `true` if `self` has a positive sign, including `+0.0`, `+∞`, and positive NaN.
    fn is_sign_positive(self) -> bool;

    /// Returns `true` if `self` has a negative sign, including `-0.0`, `-∞`, and negative NaN.
    fn is_sign_negative(self) -> bool;

    /// Converts `v` from `f64` to `Self`.
    ///
    /// May lose precision when `Self` is narrower than `f64` (e.g. [`f32`]).
    fn from_f64(v: f64) -> Self;

    /// Converts `self` to `f64`.
    ///
    /// May lose precision when `Self` has a wider precision range than `f64`.
    fn to_f64(self) -> f64;
}

macro_rules! impl_scalar {
    ($($T:ty),+ $(,)?) => {
        $(
            impl Scalar for $T {
                const ZERO: Self = 0.0;
                const ONE: Self = 1.0;
                const EPSILON: Self = <$T>::EPSILON;

                #[inline]
                fn abs(self) -> Self { <$T>::abs(self) }

                #[inline]
                fn min(self, other: Self) -> Self { <$T>::min(self, other) }

                #[inline]
                fn max(self, other: Self) -> Self { <$T>::max(self, other) }

                #[inline]
                fn clamp(self, lo: Self, hi: Self) -> Self { <$T>::clamp(self, lo, hi) }

                #[inline]
                fn signum(self) -> Self { <$T>::signum(self) }

                #[inline]
                fn copysign(self, sign: Self) -> Self { <$T>::copysign(self, sign) }

                #[inline]
                fn floor(self) -> Self { <$T>::floor(self) }

                #[inline]
                fn ceil(self) -> Self { <$T>::ceil(self) }

                #[inline]
                fn round(self) -> Self { <$T>::round(self) }

                #[inline]
                fn round_ties_even(self) -> Self { <$T>::round_ties_even(self) }

                #[inline]
                fn trunc(self) -> Self { <$T>::trunc(self) }

                #[inline]
                fn fract(self) -> Self { <$T>::fract(self) }

                #[inline]
                fn div_euclid(self, rhs: Self) -> Self { <$T>::div_euclid(self, rhs) }

                #[inline]
                fn rem_euclid(self, rhs: Self) -> Self { <$T>::rem_euclid(self, rhs) }

                #[inline]
                fn recip(self) -> Self { <$T>::recip(self) }

                #[inline]
                fn mul_add(self, a: Self, b: Self) -> Self { <$T>::mul_add(self, a, b) }

                #[inline]
                fn sqrt(self) -> Self { <$T>::sqrt(self) }

                #[inline]
                fn cbrt(self) -> Self { <$T>::cbrt(self) }

                #[inline]
                fn powi(self, n: i32) -> Self { <$T>::powi(self, n) }

                #[inline]
                fn powf(self, n: Self) -> Self { <$T>::powf(self, n) }

                #[inline]
                fn hypot(self, other: Self) -> Self { <$T>::hypot(self, other) }

                #[inline]
                fn exp(self) -> Self { <$T>::exp(self) }

                #[inline]
                fn exp2(self) -> Self { <$T>::exp2(self) }

                #[inline]
                fn exp_m1(self) -> Self { <$T>::exp_m1(self) }

                #[inline]
                fn ln(self) -> Self { <$T>::ln(self) }

                #[inline]
                fn ln_1p(self) -> Self { <$T>::ln_1p(self) }

                #[inline]
                fn log(self, base: Self) -> Self { <$T>::log(self, base) }

                #[inline]
                fn log2(self) -> Self { <$T>::log2(self) }

                #[inline]
                fn log10(self) -> Self { <$T>::log10(self) }

                #[inline]
                fn sin(self) -> Self { <$T>::sin(self) }

                #[inline]
                fn cos(self) -> Self { <$T>::cos(self) }

                #[inline]
                fn tan(self) -> Self { <$T>::tan(self) }

                #[inline]
                fn sin_cos(self) -> (Self, Self) { <$T>::sin_cos(self) }

                #[inline]
                fn asin(self) -> Self { <$T>::asin(self) }

                #[inline]
                fn acos(self) -> Self { <$T>::acos(self) }

                #[inline]
                fn atan(self) -> Self { <$T>::atan(self) }

                #[inline]
                fn atan2(self, other: Self) -> Self { <$T>::atan2(self, other) }

                #[inline]
                fn sinh(self) -> Self { <$T>::sinh(self) }

                #[inline]
                fn cosh(self) -> Self { <$T>::cosh(self) }

                #[inline]
                fn tanh(self) -> Self { <$T>::tanh(self) }

                #[inline]
                fn asinh(self) -> Self { <$T>::asinh(self) }

                #[inline]
                fn acosh(self) -> Self { <$T>::acosh(self) }

                #[inline]
                fn atanh(self) -> Self { <$T>::atanh(self) }

                #[inline]
                fn to_degrees(self) -> Self { <$T>::to_degrees(self) }

                #[inline]
                fn to_radians(self) -> Self { <$T>::to_radians(self) }

                #[inline]
                fn is_nan(self) -> bool { <$T>::is_nan(self) }

                #[inline]
                fn is_infinite(self) -> bool { <$T>::is_infinite(self) }

                #[inline]
                fn is_finite(self) -> bool { <$T>::is_finite(self) }

                #[inline]
                fn is_sign_positive(self) -> bool { <$T>::is_sign_positive(self) }

                #[inline]
                fn is_sign_negative(self) -> bool { <$T>::is_sign_negative(self) }

                #[inline]
                fn from_f64(v: f64) -> Self { v as $T }

                #[inline]
                fn to_f64(self) -> f64 { self as f64 }
            }
        )+
    };
}

impl_scalar!(f32, f64);

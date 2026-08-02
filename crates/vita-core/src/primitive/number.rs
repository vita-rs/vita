use core::fmt;
use core::iter::{Product, Sum};
use core::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

/// Numeric type carrying a physical dimension, with the operations that leave it unchanged.
///
/// Implemented for [`f32`], [`f64`], and every quantity type in [`units`](crate::units).
pub trait Quantity:
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
    + Neg<Output = Self>
    + Rem<Output = Self>
    + RemAssign
    + Mul<<Self as Quantity>::Value, Output = Self>
    + MulAssign<<Self as Quantity>::Value>
    + Div<<Self as Quantity>::Value, Output = Self>
    + DivAssign<<Self as Quantity>::Value>
    + Rem<<Self as Quantity>::Value, Output = Self>
    + RemAssign<<Self as Quantity>::Value>
    + Div<Self, Output = <Self as Quantity>::Value>
    + Sum
    + for<'a> Sum<&'a Self>
{
    /// The number underneath, read in `Self`'s own unit.
    type Value: Scalar;

    /// The additive identity, `0`.
    const ZERO: Self;

    /// Positive infinity, `∞`.
    const INFINITY: Self;

    /// Negative infinity, `-∞`.
    const NEG_INFINITY: Self;

    /// Not a Number, NaN, which compares unequal to every value including itself.
    const NAN: Self;

    /// Returns the underlying number, expressed in `Self`'s own unit.
    fn value(self) -> Self::Value;

    /// Reads `value` as a quantity in `Self`'s own unit.
    ///
    /// The inverse of [`value`](Self::value).
    fn from_value(value: Self::Value) -> Self;

    /// Returns the absolute value.
    #[inline]
    fn abs(self) -> Self {
        Self::from_value(self.value().abs())
    }

    /// Returns the minimum of `self` and `other`, ignoring NaN.
    ///
    /// If one argument is NaN, then the other argument is returned.
    #[inline]
    fn min(self, other: Self) -> Self {
        Self::from_value(self.value().min(other.value()))
    }

    /// Returns the maximum of `self` and `other`, ignoring NaN.
    ///
    /// If one argument is NaN, then the other argument is returned.
    #[inline]
    fn max(self, other: Self) -> Self {
        Self::from_value(self.value().max(other.value()))
    }

    /// Restricts `self` to the interval `[lo, hi]`.
    ///
    /// # Panics
    ///
    /// Panics if `lo > hi`, `lo` is NaN, or `hi` is NaN.
    #[inline]
    fn clamp(self, lo: Self, hi: Self) -> Self {
        Self::from_value(self.value().clamp(lo.value(), hi.value()))
    }

    /// Returns the midpoint of `self` and `other`, without intermediate overflow.
    #[inline]
    fn midpoint(self, other: Self) -> Self {
        Self::from_value(self.value().midpoint(other.value()))
    }

    /// Returns `1.0`, `-1.0`, or NaN based on the sign of `self`.
    ///
    /// A sign carries no dimension, so the result is a bare [`Value`](Self::Value).
    #[inline]
    fn signum(self) -> Self::Value {
        self.value().signum()
    }

    /// Returns a value with the magnitude of `self` and the sign of `sign`.
    #[inline]
    fn copysign(self, sign: Self) -> Self {
        Self::from_value(self.value().copysign(sign.value()))
    }

    /// Returns the largest integer less than or equal to `self`.
    #[inline]
    fn floor(self) -> Self {
        Self::from_value(self.value().floor())
    }

    /// Returns the smallest integer greater than or equal to `self`.
    #[inline]
    fn ceil(self) -> Self {
        Self::from_value(self.value().ceil())
    }

    /// Returns the nearest integer to `self`, with halves rounded away from zero.
    #[inline]
    fn round(self) -> Self {
        Self::from_value(self.value().round())
    }

    /// Returns the nearest integer to `self`, with halves rounded to even.
    #[inline]
    fn round_ties_even(self) -> Self {
        Self::from_value(self.value().round_ties_even())
    }

    /// Returns the integer part of `self`, truncated toward zero.
    #[inline]
    fn trunc(self) -> Self {
        Self::from_value(self.value().trunc())
    }

    /// Returns the fractional part of `self`.
    #[inline]
    fn fract(self) -> Self {
        Self::from_value(self.value().fract())
    }

    /// Euclidean division, the matching counterpart of [`rem_euclid`](Self::rem_euclid).
    ///
    /// A quotient of like quantities carries no dimension, so the result is a bare
    /// [`Value`](Self::Value).
    #[inline]
    fn div_euclid(self, rhs: Self) -> Self::Value {
        self.value().div_euclid(rhs.value())
    }

    /// Least nonnegative remainder of `self` divided by `rhs`.
    #[inline]
    fn rem_euclid(self, rhs: Self) -> Self {
        Self::from_value(self.value().rem_euclid(rhs.value()))
    }

    /// Computes `(self * a) + b` with a single rounding error (fused multiply-add).
    #[inline]
    fn mul_add(self, a: Self::Value, b: Self) -> Self {
        Self::from_value(self.value().mul_add(a, b.value()))
    }

    /// Returns the hypotenuse of a right triangle with legs `self` and `other`, without
    /// unnecessary overflow or underflow.
    #[inline]
    fn hypot(self, other: Self) -> Self {
        Self::from_value(self.value().hypot(other.value()))
    }

    /// Returns `true` if `self` is NaN.
    #[inline]
    fn is_nan(self) -> bool {
        self.value().is_nan()
    }

    /// Returns `true` if `self` is positive or negative infinity.
    #[inline]
    fn is_infinite(self) -> bool {
        self.value().is_infinite()
    }

    /// Returns `true` if `self` is neither NaN nor infinite.
    #[inline]
    fn is_finite(self) -> bool {
        self.value().is_finite()
    }

    /// Returns `true` if `self` is neither zero, subnormal, infinite, nor NaN.
    #[inline]
    fn is_normal(self) -> bool {
        self.value().is_normal()
    }

    /// Returns `true` if `self` is subnormal — nonzero, but held with reduced precision.
    #[inline]
    fn is_subnormal(self) -> bool {
        self.value().is_subnormal()
    }

    /// Returns `true` if `self` has a positive sign, including `+0`, `+∞`, and positive NaN.
    #[inline]
    fn is_sign_positive(self) -> bool {
        self.value().is_sign_positive()
    }

    /// Returns `true` if `self` has a negative sign, including `-0`, `-∞`, and negative NaN.
    #[inline]
    fn is_sign_negative(self) -> bool {
        self.value().is_sign_negative()
    }
}

/// Numeric scalar type for values in physical-quantity and tensor types.
///
/// Implemented for [`f32`] and [`f64`].
pub trait Scalar: Quantity<Value = Self> + Product + for<'a> Product<&'a Self> {
    /// The multiplicative identity, `1.0`.
    const ONE: Self;

    /// Machine epsilon — the difference between `1` and the next larger representable value.
    const EPSILON: Self;

    /// The smallest (most negative) finite value.
    const MIN: Self;

    /// The largest finite value.
    const MAX: Self;

    /// The smallest positive value held at full precision.
    ///
    /// Below it the representation turns subnormal and loses significant digits, so a
    /// squared magnitude smaller than this is no longer a reliable comparison.
    const MIN_POSITIVE: Self;

    /// Returns the reciprocal `1 / self`.
    fn recip(self) -> Self;

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
            impl Quantity for $T {
                type Value = Self;

                const ZERO: Self = 0.0;
                const INFINITY: Self = <$T>::INFINITY;
                const NEG_INFINITY: Self = <$T>::NEG_INFINITY;
                const NAN: Self = <$T>::NAN;

                #[inline]
                fn value(self) -> Self { self }

                #[inline]
                fn from_value(value: Self) -> Self { value }

                #[inline]
                fn abs(self) -> Self { <$T>::abs(self) }

                #[inline]
                fn min(self, other: Self) -> Self { <$T>::min(self, other) }

                #[inline]
                fn max(self, other: Self) -> Self { <$T>::max(self, other) }

                #[inline]
                fn clamp(self, lo: Self, hi: Self) -> Self { <$T>::clamp(self, lo, hi) }

                #[inline]
                fn midpoint(self, other: Self) -> Self { <$T>::midpoint(self, other) }

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
                fn mul_add(self, a: Self, b: Self) -> Self { <$T>::mul_add(self, a, b) }

                #[inline]
                fn hypot(self, other: Self) -> Self { <$T>::hypot(self, other) }

                #[inline]
                fn is_nan(self) -> bool { <$T>::is_nan(self) }

                #[inline]
                fn is_infinite(self) -> bool { <$T>::is_infinite(self) }

                #[inline]
                fn is_finite(self) -> bool { <$T>::is_finite(self) }

                #[inline]
                fn is_normal(self) -> bool { <$T>::is_normal(self) }

                #[inline]
                fn is_subnormal(self) -> bool { <$T>::is_subnormal(self) }

                #[inline]
                fn is_sign_positive(self) -> bool { <$T>::is_sign_positive(self) }

                #[inline]
                fn is_sign_negative(self) -> bool { <$T>::is_sign_negative(self) }
            }

            impl Scalar for $T {
                const ONE: Self = 1.0;
                const EPSILON: Self = <$T>::EPSILON;
                const MIN: Self = <$T>::MIN;
                const MAX: Self = <$T>::MAX;
                const MIN_POSITIVE: Self = <$T>::MIN_POSITIVE;

                #[inline]
                fn recip(self) -> Self { <$T>::recip(self) }

                #[inline]
                fn sqrt(self) -> Self { <$T>::sqrt(self) }

                #[inline]
                fn cbrt(self) -> Self { <$T>::cbrt(self) }

                #[inline]
                fn powi(self, n: i32) -> Self { <$T>::powi(self, n) }

                #[inline]
                fn powf(self, n: Self) -> Self { <$T>::powf(self, n) }

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
                fn from_f64(v: f64) -> Self { v as $T }

                #[inline]
                fn to_f64(self) -> f64 { self as f64 }
            }
        )+
    };
}

impl_scalar!(f32, f64);

#[cfg(test)]
mod tests {
    use super::*;

    use core::f64::consts::{E, FRAC_PI_2, FRAC_PI_4, PI};

    use crate::units::length::{Angstrom, Length, Nanometer};

    const TOL: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL
    }

    fn length(value: f64) -> Length<f64, Angstrom> {
        Length::new(value)
    }

    #[test]
    fn a_scalar_is_its_own_value() {
        assert_eq!(Quantity::value(2.5_f64), 2.5);
    }

    #[test]
    fn the_value_of_a_quantity_is_its_number_in_its_own_unit() {
        assert_eq!(length(2.5).value(), 2.5);
    }

    #[test]
    fn zero_leaves_a_quantity_unchanged_under_addition() {
        assert_eq!(length(2.5) + Length::<f64, Angstrom>::ZERO, length(2.5));
    }

    #[test]
    fn one_leaves_a_scalar_unchanged_under_multiplication() {
        assert_eq!(2.5 * <f64 as Scalar>::ONE, 2.5);
    }

    #[test]
    fn a_double_on_the_left_scales_a_quantity() {
        assert_eq!(2.0 * length(2.5), length(5.0));
    }

    #[test]
    fn a_single_on_the_left_scales_a_quantity() {
        let single = Length::<f32, Angstrom>::new(2.5);
        assert_eq!(2.0_f32 * single, Length::new(5.0_f32));
    }

    #[test]
    fn abs_of_a_negative_quantity_is_its_magnitude() {
        assert_eq!(length(-2.5).abs(), length(2.5));
    }

    #[test]
    fn min_returns_the_smaller_of_two_quantities() {
        assert_eq!(length(3.0).min(length(1.0)), length(1.0));
    }

    #[test]
    fn max_returns_the_larger_of_two_quantities() {
        assert_eq!(length(3.0).max(length(1.0)), length(3.0));
    }

    #[test]
    fn clamp_raises_a_quantity_below_the_interval_to_its_lower_bound() {
        assert_eq!(length(-1.0).clamp(length(0.0), length(3.0)), length(0.0));
    }

    #[test]
    fn clamp_leaves_a_quantity_inside_the_interval_unchanged() {
        assert_eq!(length(1.5).clamp(length(0.0), length(3.0)), length(1.5));
    }

    #[test]
    fn clamp_lowers_a_quantity_above_the_interval_to_its_upper_bound() {
        assert_eq!(length(5.0).clamp(length(0.0), length(3.0)), length(3.0));
    }

    #[test]
    fn midpoint_of_two_quantities_lies_halfway_between_them() {
        assert_eq!(length(1.0).midpoint(length(4.0)), length(2.5));
    }

    #[test]
    fn signum_of_a_positive_quantity_is_one() {
        assert_eq!(length(2.0).signum(), 1.0);
    }

    #[test]
    fn signum_of_a_negative_quantity_is_minus_one() {
        assert_eq!(length(-2.0).signum(), -1.0);
    }

    #[test]
    fn copysign_keeps_the_magnitude_and_takes_the_sign_of_its_argument() {
        assert_eq!(length(3.0).copysign(length(-1.0)), length(-3.0));
    }

    #[test]
    fn floor_rounds_a_quantity_toward_negative_infinity() {
        assert_eq!(length(-2.5).floor(), length(-3.0));
    }

    #[test]
    fn ceil_rounds_a_quantity_toward_positive_infinity() {
        assert_eq!(length(2.1).ceil(), length(3.0));
    }

    #[test]
    fn round_sends_a_half_away_from_zero() {
        assert_eq!(length(2.5).round(), length(3.0));
    }

    #[test]
    fn round_ties_even_sends_a_half_to_the_even_integer() {
        assert_eq!(length(2.5).round_ties_even(), length(2.0));
    }

    #[test]
    fn trunc_drops_the_fractional_part_toward_zero() {
        assert_eq!(length(-2.7).trunc(), length(-2.0));
    }

    #[test]
    fn fract_keeps_only_the_fractional_part() {
        assert_eq!(length(2.75).fract(), length(0.75));
    }

    #[test]
    fn div_euclid_of_two_quantities_is_their_dimensionless_quotient() {
        assert_eq!(length(7.0).div_euclid(length(3.0)), 2.0);
    }

    #[test]
    fn rem_euclid_of_a_negative_dividend_is_nonnegative() {
        assert_eq!(length(-1.0).rem_euclid(length(3.0)), length(2.0));
    }

    #[test]
    fn mul_add_scales_a_quantity_by_a_number_then_adds_a_quantity() {
        assert_eq!(length(2.0).mul_add(3.0, length(1.0)), length(7.0));
    }

    #[test]
    fn hypot_of_two_quantities_is_the_hypotenuse_they_span() {
        assert!(close(length(3.0).hypot(length(4.0)).value(), 5.0));
    }

    #[test]
    fn is_nan_holds_for_a_not_a_number_quantity() {
        assert!(length(f64::NAN).is_nan());
    }

    #[test]
    fn is_infinite_holds_for_an_infinite_quantity() {
        assert!(length(f64::INFINITY).is_infinite());
    }

    #[test]
    fn is_finite_holds_for_an_ordinary_quantity() {
        assert!(length(2.0).is_finite());
    }

    #[test]
    fn is_normal_holds_for_an_ordinary_quantity() {
        assert!(length(2.0).is_normal());
    }

    #[test]
    fn is_subnormal_holds_for_a_quantity_below_full_precision() {
        assert!(length(f64::MIN_POSITIVE / 2.0).is_subnormal());
    }

    #[test]
    fn is_sign_positive_holds_for_a_positive_quantity() {
        assert!(length(2.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative_holds_for_a_negative_quantity() {
        assert!(length(-2.0).is_sign_negative());
    }

    #[test]
    fn infinity_exceeds_every_finite_quantity() {
        assert!(Length::<f64, Angstrom>::INFINITY > length(f64::MAX));
    }

    #[test]
    fn negative_infinity_falls_below_every_finite_quantity() {
        assert!(Length::<f64, Angstrom>::NEG_INFINITY < length(f64::MIN));
    }

    #[test]
    fn the_not_a_number_quantity_is_not_a_number() {
        assert!(Length::<f64, Angstrom>::NAN.is_nan());
    }

    #[test]
    fn clamp_confines_a_scalar_to_the_interval() {
        assert_eq!(Quantity::clamp(5.0_f64, 0.0, 3.0), 3.0);
    }

    #[test]
    fn copysign_gives_a_scalar_the_sign_of_its_argument() {
        assert_eq!(Quantity::copysign(3.0_f64, -1.0), -3.0);
    }

    #[test]
    fn div_euclid_of_two_scalars_is_their_euclidean_quotient() {
        assert_eq!(Quantity::div_euclid(7.0_f64, 3.0), 2.0);
    }

    #[test]
    fn rem_euclid_of_a_negative_scalar_is_nonnegative() {
        assert_eq!(Quantity::rem_euclid(-1.0_f64, 3.0), 2.0);
    }

    #[test]
    fn mul_add_scales_a_scalar_then_adds() {
        assert_eq!(Quantity::mul_add(2.0_f64, 3.0, 1.0), 7.0);
    }

    #[test]
    fn recip_of_a_scalar_is_one_divided_by_it() {
        assert_eq!(Scalar::recip(4.0_f64), 0.25);
    }

    #[test]
    fn sqrt_of_a_square_is_its_root() {
        assert_eq!(Scalar::sqrt(9.0_f64), 3.0);
    }

    #[test]
    fn cbrt_of_a_cube_is_its_root() {
        assert!(close(Scalar::cbrt(27.0_f64), 3.0));
    }

    #[test]
    fn powi_raises_a_scalar_to_an_integer_power() {
        assert!(close(Scalar::powi(2.0_f64, 3), 8.0));
    }

    #[test]
    fn powf_raises_a_scalar_to_a_fractional_power() {
        assert!(close(Scalar::powf(9.0_f64, 0.5), 3.0));
    }

    #[test]
    fn exp_raises_e_to_a_scalar() {
        assert!(close(Scalar::exp(1.0_f64), E));
    }

    #[test]
    fn exp2_raises_two_to_a_scalar() {
        assert!(close(Scalar::exp2(3.0_f64), 8.0));
    }

    #[test]
    fn exp_m1_subtracts_one_from_the_exponential() {
        assert!(close(Scalar::exp_m1(1.0_f64), E - 1.0));
    }

    #[test]
    fn ln_of_e_is_one() {
        assert!(close(Scalar::ln(E), 1.0));
    }

    #[test]
    fn ln_1p_adds_one_before_taking_the_logarithm() {
        assert!(close(Scalar::ln_1p(E - 1.0), 1.0));
    }

    #[test]
    fn log_takes_the_logarithm_in_the_given_base() {
        assert!(close(Scalar::log(8.0_f64, 2.0), 3.0));
    }

    #[test]
    fn log2_takes_the_base_two_logarithm() {
        assert!(close(Scalar::log2(8.0_f64), 3.0));
    }

    #[test]
    fn log10_takes_the_base_ten_logarithm() {
        assert!(close(Scalar::log10(1000.0_f64), 3.0));
    }

    #[test]
    fn sin_of_a_quarter_turn_is_one() {
        assert!(close(Scalar::sin(FRAC_PI_2), 1.0));
    }

    #[test]
    fn cos_of_a_half_turn_is_minus_one() {
        assert!(close(Scalar::cos(PI), -1.0));
    }

    #[test]
    fn tan_of_an_eighth_turn_is_one() {
        assert!(close(Scalar::tan(FRAC_PI_4), 1.0));
    }

    #[test]
    fn sin_cos_returns_the_sine_and_the_cosine_together() {
        let (sine, cosine) = Scalar::sin_cos(FRAC_PI_4);
        assert!(close(sine, Scalar::sin(FRAC_PI_4)) && close(cosine, Scalar::cos(FRAC_PI_4)));
    }

    #[test]
    fn asin_of_one_is_a_quarter_turn() {
        assert!(close(Scalar::asin(1.0_f64), FRAC_PI_2));
    }

    #[test]
    fn acos_of_zero_is_a_quarter_turn() {
        assert!(close(Scalar::acos(0.0_f64), FRAC_PI_2));
    }

    #[test]
    fn atan_of_one_is_an_eighth_turn() {
        assert!(close(Scalar::atan(1.0_f64), FRAC_PI_4));
    }

    #[test]
    fn atan2_resolves_the_quadrant_from_both_coordinates() {
        assert!(close(Scalar::atan2(1.0_f64, -1.0), 3.0 * FRAC_PI_4));
    }

    #[test]
    fn sinh_is_half_the_difference_of_the_exponentials() {
        assert!(close(Scalar::sinh(1.0_f64), (E - Scalar::recip(E)) / 2.0));
    }

    #[test]
    fn cosh_is_half_the_sum_of_the_exponentials() {
        assert!(close(Scalar::cosh(1.0_f64), (E + Scalar::recip(E)) / 2.0));
    }

    #[test]
    fn tanh_is_the_hyperbolic_sine_over_the_hyperbolic_cosine() {
        let ratio = Scalar::sinh(1.0_f64) / Scalar::cosh(1.0_f64);
        assert!(close(Scalar::tanh(1.0_f64), ratio));
    }

    #[test]
    fn to_degrees_turns_a_half_turn_into_a_hundred_and_eighty() {
        assert!(close(Scalar::to_degrees(PI), 180.0));
    }

    #[test]
    fn from_f64_narrows_a_double_to_the_scalar_type() {
        assert_eq!(<f32 as Scalar>::from_f64(0.5), 0.5_f32);
    }

    #[test]
    fn to_f64_widens_a_scalar_to_a_double() {
        assert_eq!(Scalar::to_f64(0.5_f32), 0.5);
    }

    #[test]
    fn the_largest_scalar_is_finite() {
        assert!(Quantity::is_finite(<f64 as Scalar>::MAX));
    }

    #[test]
    fn the_smallest_positive_scalar_is_normal() {
        assert!(Quantity::is_normal(<f64 as Scalar>::MIN_POSITIVE));
    }

    #[test]
    fn is_nan_does_not_hold_for_an_ordinary_quantity() {
        assert!(!length(2.0).is_nan());
    }

    #[test]
    fn is_infinite_does_not_hold_for_an_ordinary_quantity() {
        assert!(!length(2.0).is_infinite());
    }

    #[test]
    fn is_finite_does_not_hold_for_an_infinite_quantity() {
        assert!(!length(f64::INFINITY).is_finite());
    }

    #[test]
    fn is_normal_does_not_hold_for_a_quantity_below_full_precision() {
        assert!(!length(f64::MIN_POSITIVE / 2.0).is_normal());
    }

    #[test]
    fn is_subnormal_does_not_hold_for_an_ordinary_quantity() {
        assert!(!length(2.0).is_subnormal());
    }

    #[test]
    fn is_sign_positive_does_not_hold_for_a_negative_quantity() {
        assert!(!length(-2.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative_does_not_hold_for_a_positive_quantity() {
        assert!(!length(2.0).is_sign_negative());
    }

    #[test]
    fn min_ignores_a_not_a_number_argument() {
        assert_eq!(length(f64::NAN).min(length(1.0)), length(1.0));
    }

    #[test]
    fn max_ignores_a_not_a_number_argument() {
        assert_eq!(length(f64::NAN).max(length(1.0)), length(1.0));
    }

    #[test]
    #[should_panic]
    fn clamp_with_an_inverted_interval_panics() {
        let _ = length(1.0).clamp(length(3.0), length(0.0));
    }

    #[test]
    fn sqrt_of_a_negative_scalar_is_not_a_number() {
        assert!(Quantity::is_nan(Scalar::sqrt(-1.0_f64)));
    }

    #[test]
    fn ln_of_zero_is_negative_infinity() {
        assert_eq!(Scalar::ln(0.0_f64), f64::NEG_INFINITY);
    }

    #[test]
    fn ln_of_a_negative_scalar_is_not_a_number() {
        assert!(Quantity::is_nan(Scalar::ln(-1.0_f64)));
    }

    #[test]
    fn asin_outside_the_unit_interval_is_not_a_number() {
        assert!(Quantity::is_nan(Scalar::asin(2.0_f64)));
    }

    #[test]
    fn acos_outside_the_unit_interval_is_not_a_number() {
        assert!(Quantity::is_nan(Scalar::acos(2.0_f64)));
    }

    #[test]
    fn from_f64_beyond_the_scalar_range_is_infinite() {
        assert!(Quantity::is_infinite(<f32 as Scalar>::from_f64(f64::MAX)));
    }

    #[test]
    fn one_plus_epsilon_differs_from_one() {
        assert_ne!(1.0 + <f64 as Scalar>::EPSILON, 1.0);
    }

    #[test]
    fn one_plus_half_of_epsilon_is_one() {
        assert_eq!(1.0 + <f64 as Scalar>::EPSILON / 2.0, 1.0);
    }

    #[test]
    fn twice_the_largest_scalar_is_infinite() {
        assert!(Quantity::is_infinite(<f64 as Scalar>::MAX * 2.0));
    }

    #[test]
    fn half_of_the_smallest_positive_scalar_is_subnormal() {
        assert!(Quantity::is_subnormal(<f64 as Scalar>::MIN_POSITIVE / 2.0));
    }

    #[test]
    fn signum_of_a_negative_zero_quantity_is_minus_one() {
        assert_eq!(length(-0.0).signum(), -1.0);
    }

    #[test]
    fn is_sign_negative_holds_for_a_negative_zero_quantity() {
        assert!(length(-0.0).is_sign_negative());
    }

    #[test]
    fn zero_is_the_same_quantity_in_every_unit() {
        let converted = Length::<f64, Angstrom>::ZERO.to::<Nanometer>();
        assert_eq!(converted, Length::<f64, Nanometer>::ZERO);
    }

    #[test]
    fn infinity_is_the_same_quantity_in_every_unit() {
        let converted = Length::<f64, Angstrom>::INFINITY.to::<Nanometer>();
        assert_eq!(converted, Length::<f64, Nanometer>::INFINITY);
    }

    #[test]
    fn negative_infinity_is_the_same_quantity_in_every_unit() {
        let converted = Length::<f64, Angstrom>::NEG_INFINITY.to::<Nanometer>();
        assert_eq!(converted, Length::<f64, Nanometer>::NEG_INFINITY);
    }

    #[test]
    fn a_not_a_number_quantity_stays_not_a_number_in_every_unit() {
        assert!(Length::<f64, Angstrom>::NAN.to::<Nanometer>().is_nan());
    }

    #[test]
    fn flooring_a_quantity_depends_on_the_unit_it_is_read_in() {
        let angstroms = length(27.0);
        assert_eq!(angstroms.floor(), length(27.0));
        assert_eq!(
            angstroms.to::<Nanometer>().floor().to::<Angstrom>(),
            length(20.0)
        );
    }

    #[test]
    fn a_quantity_is_recovered_from_its_own_value() {
        assert_eq!(
            Length::<f64, Angstrom>::from_value(length(2.5).value()),
            length(2.5)
        );
    }

    #[test]
    fn the_smallest_scalar_is_the_negation_of_the_largest() {
        assert_eq!(<f64 as Scalar>::MIN, -<f64 as Scalar>::MAX);
    }

    #[test]
    fn recip_is_its_own_inverse() {
        assert!(close(Scalar::recip(Scalar::recip(4.0_f64)), 4.0));
    }

    #[test]
    fn asinh_is_the_inverse_of_the_hyperbolic_sine() {
        assert!(close(Scalar::asinh(Scalar::sinh(0.75_f64)), 0.75));
    }

    #[test]
    fn acosh_is_the_inverse_of_the_hyperbolic_cosine() {
        assert!(close(Scalar::acosh(Scalar::cosh(0.75_f64)), 0.75));
    }

    #[test]
    fn atanh_is_the_inverse_of_the_hyperbolic_tangent() {
        assert!(close(Scalar::atanh(Scalar::tanh(0.75_f64)), 0.75));
    }

    #[test]
    fn to_radians_is_the_inverse_of_to_degrees() {
        assert!(close(
            Scalar::to_radians(Scalar::to_degrees(0.75_f64)),
            0.75
        ));
    }

    #[test]
    fn from_f64_and_to_f64_round_trip_a_double() {
        assert_eq!(Scalar::to_f64(<f64 as Scalar>::from_f64(2.5)), 2.5);
    }

    #[test]
    fn hypot_is_symmetric_in_its_arguments() {
        assert!(close(
            length(3.0).hypot(length(4.0)).value(),
            length(4.0).hypot(length(3.0)).value()
        ));
    }
}

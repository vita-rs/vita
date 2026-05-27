use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Numeric scalar type for values in physical-quantity types.
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
{
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
                #[inline]
                fn abs(self) -> Self { <$T>::abs(self) }

                #[inline]
                fn min(self, other: Self) -> Self { <$T>::min(self, other) }

                #[inline]
                fn max(self, other: Self) -> Self { <$T>::max(self, other) }

                #[inline]
                fn clamp(self, lo: Self, hi: Self) -> Self { <$T>::clamp(self, lo, hi) }

                #[inline]
                fn from_f64(v: f64) -> Self { v as $T }

                #[inline]
                fn to_f64(self) -> f64 { self as f64 }
            }
        )+
    };
}

impl_scalar!(f32, f64);

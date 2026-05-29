//! Length quantities and unit markers.
//!
//! The canonical unit is the **ångström** (Å).
//!
//! | Type | Symbol | Å per unit |
//! |---|---|---|
//! | [`Angstrom`] | Å | 1 |
//! | [`Bohr`] | a₀ | 0.529177210544 |
//! | [`Nanometer`] | nm | 10 |
//! | [`Picometer`] | pm | 0.01 |
//! | [`Meter`] | m | 1e10 |

use crate::units::quantity::define_quantity;

/// Marker trait for length units.
///
/// Implement this on a zero-sized type to define a new length unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of ångströms
/// per one unit of `Self`.
pub trait LengthUnit {
    /// Ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å"`, `"nm"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A length parameterised by scalar type `V` and unit marker `U`.
    Length,
    LengthUnit
);

/// The ångström (Å) — canonical length unit.
///
/// 1 Å = 1e-10 m.
pub struct Angstrom;

impl LengthUnit for Angstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å";
}

/// The bohr (a₀) — atomic unit of length (CODATA 2022).
///
/// 1 a₀ ≈ 0.529177210544 Å.
pub struct Bohr;

impl LengthUnit for Bohr {
    const TO_CANONICAL: f64 = 0.529_177_210_544;
    const SYMBOL: &'static str = "a₀";
}

/// The nanometre (nm).
///
/// 1 nm = 10 Å.
pub struct Nanometer;

impl LengthUnit for Nanometer {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "nm";
}

/// The picometre (pm).
///
/// 1 pm = 0.01 Å.
pub struct Picometer;

impl LengthUnit for Picometer {
    const TO_CANONICAL: f64 = 0.01;
    const SYMBOL: &'static str = "pm";
}

/// The metre (m) — SI base unit of length.
///
/// 1 m = 1e10 Å.
pub struct Meter;

impl LengthUnit for Meter {
    const TO_CANONICAL: f64 = 1e10;
    const SYMBOL: &'static str = "m";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Length::<f64, Angstrom>::new(1.52).value(), 1.52);
    }

    #[test]
    fn from_scalar() {
        let l: Length<f64, Nanometer> = Length::from(3.0);
        assert_eq!(l.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Length::<f64, Angstrom>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Length::<f64, Angstrom>::new(2.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn angstrom_to_nanometer() {
        let nm: Length<f64, Nanometer> = Length::<f64, Angstrom>::new(10.0).to();
        assert!((nm.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn nanometer_to_angstrom() {
        let a: Length<f64, Angstrom> = Length::<f64, Nanometer>::new(1.0).to();
        assert!((a.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn angstrom_to_picometer() {
        let pm: Length<f64, Picometer> = Length::<f64, Angstrom>::new(1.0).to();
        assert!((pm.value() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn angstrom_to_meter() {
        let m: Length<f64, Meter> = Length::<f64, Angstrom>::new(1e10).to();
        assert!((m.value() - 1.0).abs() < 1e-4);
    }

    #[test]
    fn bohr_to_angstrom() {
        let a: Length<f64, Angstrom> = Length::<f64, Bohr>::new(1.0).to();
        assert!((a.value() - 0.529_177_210_544).abs() < 1e-12);
    }

    #[test]
    fn roundtrip_nanometer_bohr_nanometer() {
        let orig = Length::<f64, Nanometer>::new(0.5);
        let back: Length<f64, Nanometer> = orig.to::<Bohr>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Length::<f64, Angstrom>::new(1.0) + Length::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut l = Length::<f64, Angstrom>::new(1.0);
        l += Length::new(0.5);
        assert_eq!(l.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Length::<f64, Angstrom>::new(3.0) - Length::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut l = Length::<f64, Angstrom>::new(3.0);
        l -= Length::new(1.0);
        assert_eq!(l.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = Length::<f64, Angstrom>::new(7.0) % Length::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut l = Length::<f64, Angstrom>::new(7.0);
        l %= Length::new(3.0);
        assert_eq!(l.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Length::<f64, Angstrom>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Length::<f64, Angstrom>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut l = Length::<f64, Angstrom>::new(2.0);
        l *= 3.0;
        assert_eq!(l.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Length::<f64, Angstrom>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut l = Length::<f64, Angstrom>::new(6.0);
        l /= 2.0;
        assert_eq!(l.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Length::<f64, Angstrom>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut l = Length::<f64, Angstrom>::new(7.0);
        l %= 3.0;
        assert_eq!(l.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Length::<f64, Angstrom>::new(6.0) / Length::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Length::<f64, Angstrom>::new(1.0);
        assert_eq!(a, Length::new(1.0));
        assert_ne!(a, Length::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Length::<f64, Angstrom>::new(1.0);
        let b = Length::<f64, Angstrom>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Length::<f64, Angstrom>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Length::<f64, Angstrom>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let l = Length::<f64, Angstrom>::new(1.0);
        let nan = Length::<f64, Angstrom>::new(f64::NAN);
        assert_eq!(l.min(nan).value(), 1.0);
        assert_eq!(nan.min(l).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let l = Length::<f64, Angstrom>::new(1.0);
        let nan = Length::<f64, Angstrom>::new(f64::NAN);
        assert_eq!(l.max(nan).value(), 1.0);
        assert_eq!(nan.max(l).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Length::<f64, Angstrom>::new(1.0);
        let hi = Length::<f64, Angstrom>::new(2.0);
        assert_eq!(Length::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Length::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Length::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Length::<f64, Angstrom>::new(2.0);
        let hi = Length::<f64, Angstrom>::new(1.0);
        Length::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(Length::<f64, Angstrom>::new(3.0).signum(), 1.0);
        assert_eq!(Length::<f64, Angstrom>::new(-3.0).signum(), -1.0);
    }

    #[test]
    fn copysign() {
        let l = Length::<f64, Angstrom>::new(3.0);
        let sign = Length::<f64, Angstrom>::new(-1.0);
        assert_eq!(l.copysign(sign).value(), -3.0);
        assert_eq!((-l).copysign(l).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(Length::<f64, Angstrom>::new(2.7).floor().value(), 2.0);
        assert_eq!(Length::<f64, Angstrom>::new(-2.3).floor().value(), -3.0);
    }

    #[test]
    fn ceil() {
        assert_eq!(Length::<f64, Angstrom>::new(2.3).ceil().value(), 3.0);
        assert_eq!(Length::<f64, Angstrom>::new(-2.7).ceil().value(), -2.0);
    }

    #[test]
    fn round() {
        assert_eq!(Length::<f64, Angstrom>::new(2.5).round().value(), 3.0);
        assert_eq!(Length::<f64, Angstrom>::new(-2.5).round().value(), -3.0);
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            Length::<f64, Angstrom>::new(2.5).round_ties_even().value(),
            2.0
        );
        assert_eq!(
            Length::<f64, Angstrom>::new(3.5).round_ties_even().value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(Length::<f64, Angstrom>::new(2.7).trunc().value(), 2.0);
        assert_eq!(Length::<f64, Angstrom>::new(-2.7).trunc().value(), -2.0);
    }

    #[test]
    fn fract() {
        assert!((Length::<f64, Angstrom>::new(2.75).fract().value() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn div_euclid() {
        let q = Length::<f64, Angstrom>::new(7.0).div_euclid(Length::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Length::<f64, Angstrom>::new(-7.0).rem_euclid(Length::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Length::<f64, Angstrom>::new(2.0).mul_add(3.0, Length::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Length::<f64, Angstrom>::new(3.0).hypot(Length::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Length::<f64, Angstrom>::new(f64::NAN).is_nan());
        assert!(!Length::<f64, Angstrom>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Length::<f64, Angstrom>::new(f64::INFINITY).is_infinite());
        assert!(!Length::<f64, Angstrom>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Length::<f64, Angstrom>::new(1.0).is_finite());
        assert!(!Length::<f64, Angstrom>::new(f64::INFINITY).is_finite());
        assert!(!Length::<f64, Angstrom>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Length::<f64, Angstrom>::new(1.0).is_sign_positive());
        assert!(!Length::<f64, Angstrom>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Length::<f64, Angstrom>::new(-1.0).is_sign_negative());
        assert!(!Length::<f64, Angstrom>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let v = [
            Length::<f64, Angstrom>::new(1.0),
            Length::new(2.0),
            Length::new(3.0),
        ];
        let total: Length<f64, Angstrom> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Length::<f64, Angstrom>::new(1.0),
            Length::new(2.0),
            Length::new(3.0),
        ];
        let total: Length<f64, Angstrom> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Length<f64, Angstrom> = iter::empty::<Length<f64, Angstrom>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Length::<f64, Nanometer>::new(1.5).to_string(), "1.5 nm");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Length::<f64, Angstrom>::new(1.0)),
            "Length(1.0)"
        );
    }

    #[test]
    fn f32_angstrom_to_nanometer() {
        let nm: Length<f32, Nanometer> = Length::<f32, Angstrom>::new(10.0_f32).to();
        assert!((nm.value() - 1.0_f32).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Length::<f32, Angstrom>::new(1.0_f32) + Length::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

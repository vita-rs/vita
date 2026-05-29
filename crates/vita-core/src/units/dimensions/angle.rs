//! Angle quantities and unit markers.
//!
//! The canonical unit is the **radian** (rad).
//!
//! | Type | Symbol | rad per unit |
//! |---|---|---|
//! | [`Radian`] | rad | 1 |
//! | [`Degree`] | ° | π / 180 |
//! | [`Milliradian`] | mrad | 1e-3 |
//! | [`Revolution`] | rev | 2π |

use crate::units::quantity::define_quantity;

/// Marker trait for angle units.
///
/// Implement this on a zero-sized type to define a new angle unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of radians
/// per one unit of `Self`.
pub trait AngleUnit {
    /// Radians per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"rad"`, `"°"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An angle parameterized by scalar type `V` and unit marker `U`.
    Angle,
    AngleUnit
);

/// The radian (rad) — canonical angle unit.
///
/// 1 rad is the SI unit of plane angle.
pub struct Radian;

impl AngleUnit for Radian {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "rad";
}

/// The degree (°).
///
/// 1° = π / 180 rad.
pub struct Degree;

impl AngleUnit for Degree {
    const TO_CANONICAL: f64 = core::f64::consts::PI / 180.0;
    const SYMBOL: &'static str = "°";
}

/// The milliradian (mrad).
///
/// 1 mrad = 1e-3 rad.
pub struct Milliradian;

impl AngleUnit for Milliradian {
    const TO_CANONICAL: f64 = 1e-3;
    const SYMBOL: &'static str = "mrad";
}

/// The revolution (rev) — full turn.
///
/// 1 rev = 2π rad.
pub struct Revolution;

impl AngleUnit for Revolution {
    const TO_CANONICAL: f64 = 2.0 * core::f64::consts::PI;
    const SYMBOL: &'static str = "rev";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Angle::<f64, Radian>::new(1.52).value(), 1.52);
    }

    #[test]
    fn from_scalar() {
        let a: Angle<f64, Degree> = Angle::from(3.0);
        assert_eq!(a.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Angle::<f64, Radian>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Angle::<f64, Radian>::new(2.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn radian_to_degree() {
        let d: Angle<f64, Degree> = Angle::<f64, Radian>::new(1.0).to();
        assert!((d.value() - 57.295_779_513_082_32).abs() < 1e-12);
    }

    #[test]
    fn degree_to_radian() {
        let r: Angle<f64, Radian> = Angle::<f64, Degree>::new(1.0).to();
        assert!((r.value() - 0.017_453_292_519_943_295).abs() < 1e-17);
    }

    #[test]
    fn radian_to_milliradian() {
        let mr: Angle<f64, Milliradian> = Angle::<f64, Radian>::new(1.0).to();
        assert!((mr.value() - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn degree_to_revolution() {
        let rev: Angle<f64, Revolution> = Angle::<f64, Degree>::new(360.0).to();
        assert!((rev.value() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn revolution_to_radian() {
        let r: Angle<f64, Radian> = Angle::<f64, Revolution>::new(1.0).to();
        assert!((r.value() - ::core::f64::consts::TAU).abs() < 1e-14);
    }

    #[test]
    fn roundtrip_degree_revolution_degree() {
        let orig = Angle::<f64, Degree>::new(90.0);
        let back: Angle<f64, Degree> = orig.to::<Revolution>().to();
        assert!((back.value() - 90.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Angle::<f64, Radian>::new(1.0) + Angle::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut a = Angle::<f64, Radian>::new(1.0);
        a += Angle::new(0.5);
        assert_eq!(a.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Angle::<f64, Radian>::new(3.0) - Angle::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut a = Angle::<f64, Radian>::new(3.0);
        a -= Angle::new(1.0);
        assert_eq!(a.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = Angle::<f64, Radian>::new(7.0) % Angle::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut a = Angle::<f64, Radian>::new(7.0);
        a %= Angle::new(3.0);
        assert_eq!(a.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Angle::<f64, Radian>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Angle::<f64, Radian>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut a = Angle::<f64, Radian>::new(2.0);
        a *= 3.0;
        assert_eq!(a.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Angle::<f64, Radian>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut a = Angle::<f64, Radian>::new(6.0);
        a /= 2.0;
        assert_eq!(a.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Angle::<f64, Radian>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut a = Angle::<f64, Radian>::new(7.0);
        a %= 3.0;
        assert_eq!(a.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Angle::<f64, Radian>::new(6.0) / Angle::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Angle::<f64, Radian>::new(1.0);
        assert_eq!(a, Angle::new(1.0));
        assert_ne!(a, Angle::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Angle::<f64, Radian>::new(1.0);
        let b = Angle::<f64, Radian>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Angle::<f64, Radian>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Angle::<f64, Radian>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let a = Angle::<f64, Radian>::new(1.0);
        let nan = Angle::<f64, Radian>::new(f64::NAN);
        assert_eq!(a.min(nan).value(), 1.0);
        assert_eq!(nan.min(a).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let a = Angle::<f64, Radian>::new(1.0);
        let nan = Angle::<f64, Radian>::new(f64::NAN);
        assert_eq!(a.max(nan).value(), 1.0);
        assert_eq!(nan.max(a).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Angle::<f64, Radian>::new(1.0);
        let hi = Angle::<f64, Radian>::new(2.0);
        assert_eq!(Angle::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Angle::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Angle::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Angle::<f64, Radian>::new(2.0);
        let hi = Angle::<f64, Radian>::new(1.0);
        Angle::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(Angle::<f64, Radian>::new(3.0).signum(), 1.0);
        assert_eq!(Angle::<f64, Radian>::new(-3.0).signum(), -1.0);
    }

    #[test]
    fn copysign() {
        let a = Angle::<f64, Radian>::new(3.0);
        let sign = Angle::<f64, Radian>::new(-1.0);
        assert_eq!(a.copysign(sign).value(), -3.0);
        assert_eq!((-a).copysign(a).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(Angle::<f64, Radian>::new(2.7).floor().value(), 2.0);
        assert_eq!(Angle::<f64, Radian>::new(-2.3).floor().value(), -3.0);
    }

    #[test]
    fn ceil() {
        assert_eq!(Angle::<f64, Radian>::new(2.3).ceil().value(), 3.0);
        assert_eq!(Angle::<f64, Radian>::new(-2.7).ceil().value(), -2.0);
    }

    #[test]
    fn round() {
        assert_eq!(Angle::<f64, Radian>::new(2.5).round().value(), 3.0);
        assert_eq!(Angle::<f64, Radian>::new(-2.5).round().value(), -3.0);
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            Angle::<f64, Radian>::new(2.5).round_ties_even().value(),
            2.0
        );
        assert_eq!(
            Angle::<f64, Radian>::new(3.5).round_ties_even().value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(Angle::<f64, Radian>::new(2.7).trunc().value(), 2.0);
        assert_eq!(Angle::<f64, Radian>::new(-2.7).trunc().value(), -2.0);
    }

    #[test]
    fn fract() {
        assert!((Angle::<f64, Radian>::new(2.75).fract().value() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn div_euclid() {
        let q = Angle::<f64, Radian>::new(7.0).div_euclid(Angle::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Angle::<f64, Radian>::new(-7.0).rem_euclid(Angle::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Angle::<f64, Radian>::new(2.0).mul_add(3.0, Angle::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Angle::<f64, Radian>::new(3.0).hypot(Angle::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Angle::<f64, Radian>::new(f64::NAN).is_nan());
        assert!(!Angle::<f64, Radian>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Angle::<f64, Radian>::new(f64::INFINITY).is_infinite());
        assert!(!Angle::<f64, Radian>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Angle::<f64, Radian>::new(1.0).is_finite());
        assert!(!Angle::<f64, Radian>::new(f64::INFINITY).is_finite());
        assert!(!Angle::<f64, Radian>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Angle::<f64, Radian>::new(1.0).is_sign_positive());
        assert!(!Angle::<f64, Radian>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Angle::<f64, Radian>::new(-1.0).is_sign_negative());
        assert!(!Angle::<f64, Radian>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let v = [
            Angle::<f64, Radian>::new(1.0),
            Angle::new(2.0),
            Angle::new(3.0),
        ];
        let total: Angle<f64, Radian> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Angle::<f64, Radian>::new(1.0),
            Angle::new(2.0),
            Angle::new(3.0),
        ];
        let total: Angle<f64, Radian> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Angle<f64, Radian> = iter::empty::<Angle<f64, Radian>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Angle::<f64, Degree>::new(1.5).to_string(), "1.5 °");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Angle::<f64, Radian>::new(1.0)),
            "Angle(1.0)"
        );
    }

    #[test]
    fn f32_degree_to_radian() {
        let r: Angle<f32, Radian> = Angle::<f32, Degree>::new(180.0_f32).to();
        assert!((r.value() - core::f32::consts::PI).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Angle::<f32, Radian>::new(1.0_f32) + Angle::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

//! Temperature quantities and unit markers.
//!
//! The canonical unit is the **kelvin** (K).
//!
//! | Type | Symbol | K per unit |
//! |---|---|---|
//! | [`Kelvin`] | K | 1 |
//! | [`Rankine`] | °R | 5/9 |

use crate::units::quantity::define_quantity;

/// Marker trait for temperature units.
///
/// Implement this on a zero-sized type to define a new temperature unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of kelvins
/// per one unit of `Self`.
pub trait TemperatureUnit {
    /// Kelvins per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"K"`, `"°R"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A temperature parameterised by scalar type `V` and unit marker `U`.
    Temperature,
    TemperatureUnit
);

/// The kelvin (K) — canonical temperature unit.
///
/// SI base unit of thermodynamic temperature.
pub struct Kelvin;

impl TemperatureUnit for Kelvin {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "K";
}

/// The rankine (°R) — absolute temperature scale with Fahrenheit-sized degrees.
///
/// 1 °R = 5/9 K exactly.
pub struct Rankine;

impl TemperatureUnit for Rankine {
    const TO_CANONICAL: f64 = 5.0 / 9.0;
    const SYMBOL: &'static str = "°R";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Temperature::<f64, Kelvin>::new(300.0).value(), 300.0);
    }

    #[test]
    fn from_scalar() {
        let t: Temperature<f64, Kelvin> = Temperature::from(300.0);
        assert_eq!(t.value(), 300.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Temperature::<f64, Kelvin>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Temperature::<f64, Kelvin>::new(300.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn kelvin_to_rankine() {
        let r: Temperature<f64, Rankine> = Temperature::<f64, Kelvin>::new(1.0).to();
        assert!((r.value() - 1.8).abs() < 1e-12);
    }

    #[test]
    fn rankine_to_kelvin() {
        let k: Temperature<f64, Kelvin> = Temperature::<f64, Rankine>::new(1.8).to();
        assert!((k.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn roundtrip_rankine_kelvin_rankine() {
        let orig = Temperature::<f64, Rankine>::new(540.0);
        let back: Temperature<f64, Rankine> = orig.to::<Kelvin>().to();
        assert!((back.value() - 540.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Temperature::<f64, Kelvin>::new(100.0) + Temperature::new(200.0);
        assert_eq!(sum.value(), 300.0);
    }

    #[test]
    fn add_assign() {
        let mut t = Temperature::<f64, Kelvin>::new(200.0);
        t += Temperature::new(100.0);
        assert_eq!(t.value(), 300.0);
    }

    #[test]
    fn sub() {
        let diff = Temperature::<f64, Kelvin>::new(300.0) - Temperature::new(100.0);
        assert_eq!(diff.value(), 200.0);
    }

    #[test]
    fn sub_assign() {
        let mut t = Temperature::<f64, Kelvin>::new(300.0);
        t -= Temperature::new(100.0);
        assert_eq!(t.value(), 200.0);
    }

    #[test]
    fn rem() {
        let r = Temperature::<f64, Kelvin>::new(7.0) % Temperature::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut t = Temperature::<f64, Kelvin>::new(7.0);
        t %= Temperature::new(3.0);
        assert_eq!(t.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Temperature::<f64, Kelvin>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(
            (Temperature::<f64, Kelvin>::new(200.0) * 1.5).value(),
            300.0
        );
    }

    #[test]
    fn mul_assign_scalar() {
        let mut t = Temperature::<f64, Kelvin>::new(200.0);
        t *= 1.5;
        assert_eq!(t.value(), 300.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!(
            (Temperature::<f64, Kelvin>::new(300.0) / 3.0).value(),
            100.0
        );
    }

    #[test]
    fn div_assign_scalar() {
        let mut t = Temperature::<f64, Kelvin>::new(300.0);
        t /= 3.0;
        assert_eq!(t.value(), 100.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Temperature::<f64, Kelvin>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut t = Temperature::<f64, Kelvin>::new(7.0);
        t %= 3.0;
        assert_eq!(t.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Temperature::<f64, Kelvin>::new(600.0) / Temperature::new(200.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Temperature::<f64, Kelvin>::new(300.0);
        assert_eq!(a, Temperature::new(300.0));
        assert_ne!(a, Temperature::new(200.0));
    }

    #[test]
    fn ord() {
        let a = Temperature::<f64, Kelvin>::new(200.0);
        let b = Temperature::<f64, Kelvin>::new(300.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Temperature::<f64, Kelvin>::new(-300.0).abs().value(), 300.0);
        assert_eq!(Temperature::<f64, Kelvin>::new(300.0).abs().value(), 300.0);
    }

    #[test]
    fn min_ignores_nan() {
        let t = Temperature::<f64, Kelvin>::new(300.0);
        let nan = Temperature::<f64, Kelvin>::new(f64::NAN);
        assert_eq!(t.min(nan).value(), 300.0);
        assert_eq!(nan.min(t).value(), 300.0);
    }

    #[test]
    fn max_ignores_nan() {
        let t = Temperature::<f64, Kelvin>::new(300.0);
        let nan = Temperature::<f64, Kelvin>::new(f64::NAN);
        assert_eq!(t.max(nan).value(), 300.0);
        assert_eq!(nan.max(t).value(), 300.0);
    }

    #[test]
    fn clamp() {
        let lo = Temperature::<f64, Kelvin>::new(200.0);
        let hi = Temperature::<f64, Kelvin>::new(400.0);
        assert_eq!(Temperature::new(300.0_f64).clamp(lo, hi).value(), 300.0);
        assert_eq!(Temperature::new(100.0_f64).clamp(lo, hi).value(), 200.0);
        assert_eq!(Temperature::new(500.0_f64).clamp(lo, hi).value(), 400.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Temperature::<f64, Kelvin>::new(400.0);
        let hi = Temperature::<f64, Kelvin>::new(200.0);
        Temperature::new(300.0_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(Temperature::<f64, Kelvin>::new(3.0).signum(), 1.0);
        assert_eq!(Temperature::<f64, Kelvin>::new(-3.0).signum(), -1.0);
    }

    #[test]
    fn copysign() {
        let t = Temperature::<f64, Kelvin>::new(3.0);
        let sign = Temperature::<f64, Kelvin>::new(-1.0);
        assert_eq!(t.copysign(sign).value(), -3.0);
        assert_eq!((-t).copysign(t).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(Temperature::<f64, Kelvin>::new(2.7).floor().value(), 2.0);
        assert_eq!(Temperature::<f64, Kelvin>::new(-2.3).floor().value(), -3.0);
    }

    #[test]
    fn ceil() {
        assert_eq!(Temperature::<f64, Kelvin>::new(2.3).ceil().value(), 3.0);
        assert_eq!(Temperature::<f64, Kelvin>::new(-2.7).ceil().value(), -2.0);
    }

    #[test]
    fn round() {
        assert_eq!(Temperature::<f64, Kelvin>::new(2.5).round().value(), 3.0);
        assert_eq!(Temperature::<f64, Kelvin>::new(-2.5).round().value(), -3.0);
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            Temperature::<f64, Kelvin>::new(2.5)
                .round_ties_even()
                .value(),
            2.0
        );
        assert_eq!(
            Temperature::<f64, Kelvin>::new(3.5)
                .round_ties_even()
                .value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(Temperature::<f64, Kelvin>::new(2.7).trunc().value(), 2.0);
        assert_eq!(Temperature::<f64, Kelvin>::new(-2.7).trunc().value(), -2.0);
    }

    #[test]
    fn fract() {
        assert!((Temperature::<f64, Kelvin>::new(2.75).fract().value() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn div_euclid() {
        let q = Temperature::<f64, Kelvin>::new(7.0).div_euclid(Temperature::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Temperature::<f64, Kelvin>::new(-7.0).rem_euclid(Temperature::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Temperature::<f64, Kelvin>::new(2.0).mul_add(3.0, Temperature::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Temperature::<f64, Kelvin>::new(3.0).hypot(Temperature::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Temperature::<f64, Kelvin>::new(f64::NAN).is_nan());
        assert!(!Temperature::<f64, Kelvin>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Temperature::<f64, Kelvin>::new(f64::INFINITY).is_infinite());
        assert!(!Temperature::<f64, Kelvin>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Temperature::<f64, Kelvin>::new(1.0).is_finite());
        assert!(!Temperature::<f64, Kelvin>::new(f64::INFINITY).is_finite());
        assert!(!Temperature::<f64, Kelvin>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Temperature::<f64, Kelvin>::new(1.0).is_sign_positive());
        assert!(!Temperature::<f64, Kelvin>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Temperature::<f64, Kelvin>::new(-1.0).is_sign_negative());
        assert!(!Temperature::<f64, Kelvin>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let v = [
            Temperature::<f64, Kelvin>::new(100.0),
            Temperature::new(200.0),
            Temperature::new(300.0),
        ];
        let total: Temperature<f64, Kelvin> = v.iter().copied().sum();
        assert_eq!(total.value(), 600.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Temperature::<f64, Kelvin>::new(100.0),
            Temperature::new(200.0),
            Temperature::new(300.0),
        ];
        let total: Temperature<f64, Kelvin> = v.iter().sum();
        assert_eq!(total.value(), 600.0);
    }

    #[test]
    fn sum_empty() {
        let total: Temperature<f64, Kelvin> = iter::empty::<Temperature<f64, Kelvin>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Temperature::<f64, Kelvin>::new(300.0).to_string(), "300 K");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Temperature::<f64, Kelvin>::new(1.0)),
            "Temperature(1.0)"
        );
    }

    #[test]
    fn f32_kelvin_to_rankine() {
        let r: Temperature<f32, Rankine> = Temperature::<f32, Kelvin>::new(1.0_f32).to();
        assert!((r.value() - 1.8_f32).abs() < 1e-5_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Temperature::<f32, Kelvin>::new(100.0_f32) + Temperature::new(200.0_f32);
        assert_eq!(sum.value(), 300.0_f32);
    }
}

//! Mass quantities and unit markers.
//!
//! The canonical unit is the **dalton** (Da).
//!
//! | Type | Symbol | Da per unit |
//! |---|---|---|
//! | [`Dalton`] | Da | 1 |
//! | [`ElectronMass`] | mₑ | 5.485799090441e-4 |
//! | [`ProtonMass`] | mₚ | 1.0072764665789 |
//! | [`Gram`] | g | 6.0221407537e23 |
//! | [`Kilogram`] | kg | 6.0221407537e26 |

use crate::units::quantity::define_quantity;

/// Marker trait for mass units.
///
/// Implement this on a zero-sized type to define a new mass unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of daltons
/// per one unit of `Self`.
pub trait MassUnit {
    /// Daltons per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Da"`, `"kg"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A mass parameterized by scalar type `V` and unit marker `U`.
    Mass,
    MassUnit
);

/// The dalton (Da) — canonical mass unit.
///
/// 1 Da ≈ 1.66053906892e-27 kg (CODATA 2022).
pub struct Dalton;

impl MassUnit for Dalton {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Da";
}

/// The electron mass (mₑ) — atomic unit of mass (CODATA 2022).
///
/// 1 mₑ ≈ 5.485799090441e-4 Da.
pub struct ElectronMass;

impl MassUnit for ElectronMass {
    const TO_CANONICAL: f64 = 5.485_799_090_441e-4;
    const SYMBOL: &'static str = "mₑ";
}

/// The proton mass (mₚ) (CODATA 2022).
///
/// 1 mₚ ≈ 1.0072764665789 Da.
pub struct ProtonMass;

impl MassUnit for ProtonMass {
    const TO_CANONICAL: f64 = 1.007_276_466_578_9;
    const SYMBOL: &'static str = "mₚ";
}

/// The gram (g) (CODATA 2022).
///
/// 1 g ≈ 6.0221407537e23 Da.
pub struct Gram;

impl MassUnit for Gram {
    const TO_CANONICAL: f64 = 6.022_140_753_7e23;
    const SYMBOL: &'static str = "g";
}

/// The kilogram (kg) — SI base unit of mass (CODATA 2022).
///
/// 1 kg ≈ 6.0221407537e26 Da.
pub struct Kilogram;

impl MassUnit for Kilogram {
    const TO_CANONICAL: f64 = 6.022_140_753_7e26;
    const SYMBOL: &'static str = "kg";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Mass::<f64, Dalton>::new(1.52).value(), 1.52);
    }

    #[test]
    fn from_scalar() {
        let m: Mass<f64, ElectronMass> = Mass::from(3.0);
        assert_eq!(m.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Mass::<f64, Dalton>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Mass::<f64, Dalton>::new(2.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn electron_mass_to_dalton() {
        let da: Mass<f64, Dalton> = Mass::<f64, ElectronMass>::new(1.0).to();
        assert!((da.value() - 5.485_799_090_441e-4).abs() < 1e-18);
    }

    #[test]
    fn dalton_to_electron_mass() {
        let me: Mass<f64, ElectronMass> = Mass::<f64, Dalton>::new(1.0).to();
        assert!((me.value() - 1_822.888_486_278).abs() < 1e-9);
    }

    #[test]
    fn proton_mass_to_dalton() {
        let da: Mass<f64, Dalton> = Mass::<f64, ProtonMass>::new(1.0).to();
        assert!((da.value() - 1.007_276_466_578_9).abs() < 1e-15);
    }

    #[test]
    fn dalton_to_gram() {
        let g: Mass<f64, Gram> = Mass::<f64, Dalton>::new(1.0).to();
        assert!((g.value() - 1.660_539_068_92e-24).abs() < 1e-34);
    }

    #[test]
    fn gram_to_dalton() {
        let da: Mass<f64, Dalton> = Mass::<f64, Gram>::new(1.0).to();
        assert!((da.value() - 6.022_140_753_7e23).abs() < 1e11);
    }

    #[test]
    fn dalton_to_kilogram() {
        let kg: Mass<f64, Kilogram> = Mass::<f64, Dalton>::new(1.0).to();
        assert!((kg.value() - 1.660_539_068_92e-27).abs() < 1e-37);
    }

    #[test]
    fn kilogram_to_dalton() {
        let da: Mass<f64, Dalton> = Mass::<f64, Kilogram>::new(1.0).to();
        assert!((da.value() - 6.022_140_753_7e26).abs() < 1e14);
    }

    #[test]
    fn roundtrip_electron_mass_kilogram_electron_mass() {
        let orig = Mass::<f64, ElectronMass>::new(10.0);
        let back: Mass<f64, ElectronMass> = orig.to::<Kilogram>().to();
        assert!((back.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Mass::<f64, Dalton>::new(1.0) + Mass::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut m = Mass::<f64, Dalton>::new(1.0);
        m += Mass::new(0.5);
        assert_eq!(m.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Mass::<f64, Dalton>::new(3.0) - Mass::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut m = Mass::<f64, Dalton>::new(3.0);
        m -= Mass::new(1.0);
        assert_eq!(m.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = Mass::<f64, Dalton>::new(7.0) % Mass::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut m = Mass::<f64, Dalton>::new(7.0);
        m %= Mass::new(3.0);
        assert_eq!(m.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Mass::<f64, Dalton>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Mass::<f64, Dalton>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut m = Mass::<f64, Dalton>::new(2.0);
        m *= 3.0;
        assert_eq!(m.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Mass::<f64, Dalton>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut m = Mass::<f64, Dalton>::new(6.0);
        m /= 2.0;
        assert_eq!(m.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Mass::<f64, Dalton>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut m = Mass::<f64, Dalton>::new(7.0);
        m %= 3.0;
        assert_eq!(m.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Mass::<f64, Dalton>::new(6.0) / Mass::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Mass::<f64, Dalton>::new(1.0);
        assert_eq!(a, Mass::new(1.0));
        assert_ne!(a, Mass::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Mass::<f64, Dalton>::new(1.0);
        let b = Mass::<f64, Dalton>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Mass::<f64, Dalton>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Mass::<f64, Dalton>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let m = Mass::<f64, Dalton>::new(1.0);
        let nan = Mass::<f64, Dalton>::new(f64::NAN);
        assert_eq!(m.min(nan).value(), 1.0);
        assert_eq!(nan.min(m).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let m = Mass::<f64, Dalton>::new(1.0);
        let nan = Mass::<f64, Dalton>::new(f64::NAN);
        assert_eq!(m.max(nan).value(), 1.0);
        assert_eq!(nan.max(m).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Mass::<f64, Dalton>::new(1.0);
        let hi = Mass::<f64, Dalton>::new(2.0);
        assert_eq!(Mass::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Mass::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Mass::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Mass::<f64, Dalton>::new(2.0);
        let hi = Mass::<f64, Dalton>::new(1.0);
        Mass::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(Mass::<f64, Dalton>::new(3.0).signum(), 1.0);
        assert_eq!(Mass::<f64, Dalton>::new(-3.0).signum(), -1.0);
    }

    #[test]
    fn copysign() {
        let m = Mass::<f64, Dalton>::new(3.0);
        let sign = Mass::<f64, Dalton>::new(-1.0);
        assert_eq!(m.copysign(sign).value(), -3.0);
        assert_eq!((-m).copysign(m).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(Mass::<f64, Dalton>::new(2.7).floor().value(), 2.0);
        assert_eq!(Mass::<f64, Dalton>::new(-2.3).floor().value(), -3.0);
    }

    #[test]
    fn ceil() {
        assert_eq!(Mass::<f64, Dalton>::new(2.3).ceil().value(), 3.0);
        assert_eq!(Mass::<f64, Dalton>::new(-2.7).ceil().value(), -2.0);
    }

    #[test]
    fn round() {
        assert_eq!(Mass::<f64, Dalton>::new(2.5).round().value(), 3.0);
        assert_eq!(Mass::<f64, Dalton>::new(-2.5).round().value(), -3.0);
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(Mass::<f64, Dalton>::new(2.5).round_ties_even().value(), 2.0);
        assert_eq!(Mass::<f64, Dalton>::new(3.5).round_ties_even().value(), 4.0);
    }

    #[test]
    fn trunc() {
        assert_eq!(Mass::<f64, Dalton>::new(2.7).trunc().value(), 2.0);
        assert_eq!(Mass::<f64, Dalton>::new(-2.7).trunc().value(), -2.0);
    }

    #[test]
    fn fract() {
        assert!((Mass::<f64, Dalton>::new(2.75).fract().value() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn div_euclid() {
        let q = Mass::<f64, Dalton>::new(7.0).div_euclid(Mass::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Mass::<f64, Dalton>::new(-7.0).rem_euclid(Mass::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Mass::<f64, Dalton>::new(2.0).mul_add(3.0, Mass::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Mass::<f64, Dalton>::new(3.0).hypot(Mass::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Mass::<f64, Dalton>::new(f64::NAN).is_nan());
        assert!(!Mass::<f64, Dalton>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Mass::<f64, Dalton>::new(f64::INFINITY).is_infinite());
        assert!(!Mass::<f64, Dalton>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Mass::<f64, Dalton>::new(1.0).is_finite());
        assert!(!Mass::<f64, Dalton>::new(f64::INFINITY).is_finite());
        assert!(!Mass::<f64, Dalton>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Mass::<f64, Dalton>::new(1.0).is_sign_positive());
        assert!(!Mass::<f64, Dalton>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Mass::<f64, Dalton>::new(-1.0).is_sign_negative());
        assert!(!Mass::<f64, Dalton>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let v = [
            Mass::<f64, Dalton>::new(1.0),
            Mass::new(2.0),
            Mass::new(3.0),
        ];
        let total: Mass<f64, Dalton> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Mass::<f64, Dalton>::new(1.0),
            Mass::new(2.0),
            Mass::new(3.0),
        ];
        let total: Mass<f64, Dalton> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Mass<f64, Dalton> = iter::empty::<Mass<f64, Dalton>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Mass::<f64, ElectronMass>::new(1.5).to_string(), "1.5 mₑ");
    }

    #[test]
    fn debug() {
        assert_eq!(format!("{:?}", Mass::<f64, Dalton>::new(1.0)), "Mass(1.0)");
    }

    #[test]
    fn f32_electron_mass_to_dalton() {
        let da: Mass<f32, Dalton> = Mass::<f32, ElectronMass>::new(1.0_f32).to();
        assert!((da.value() - 5.485_799e-4_f32).abs() < 1e-10_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Mass::<f32, Dalton>::new(1.0_f32) + Mass::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

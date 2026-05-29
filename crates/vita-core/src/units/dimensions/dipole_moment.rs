//! Dipole-moment quantities and unit markers.
//!
//! The canonical unit is the **debye** (D).
//!
//! | Type | Symbol | D per unit |
//! |---|---|---|
//! | [`Debye`] | D | 1 |
//! | [`ElectronBohr`] | ea₀ | 8.4783536198e-30 / 3.335640951981521e-30 |
//! | [`CoulombMeter`] | C·m | 1 / 3.335640951981521e-30 |

use crate::units::quantity::define_quantity;

/// Marker trait for dipole-moment units.
///
/// Implement this on a zero-sized type to define a new dipole-moment unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of debyes
/// per one unit of `Self`.
pub trait DipoleMomentUnit {
    /// Debyes per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"D"`, `"ea₀"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A dipole moment parameterized by scalar type `V` and unit marker `U`.
    DipoleMoment,
    DipoleMomentUnit
);

/// The debye (D) — canonical dipole-moment unit.
///
/// 1 D = 3.335640951981521e-30 C·m (exact, derived from c = 299 792 458 m s⁻¹).
pub struct Debye;

impl DipoleMomentUnit for Debye {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "D";
}

/// The atomic unit of electric dipole moment (ea₀) (CODATA 2022).
///
/// 1 ea₀ = 8.4783536198e-30 C·m ≈ 2.5417464715 D.
pub struct ElectronBohr;

impl DipoleMomentUnit for ElectronBohr {
    const TO_CANONICAL: f64 = 8.478_353_619_8e-30 / 3.335_640_951_981_521e-30;
    const SYMBOL: &'static str = "ea₀";
}

/// The coulomb meter (C·m) — SI unit of electric dipole moment.
///
/// 1 C·m = 2.997924580e29 D (exact, derived from c = 299 792 458 m s⁻¹).
pub struct CoulombMeter;

impl DipoleMomentUnit for CoulombMeter {
    const TO_CANONICAL: f64 = 1.0 / 3.335_640_951_981_521e-30;
    const SYMBOL: &'static str = "C·m";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(DipoleMoment::<f64, Debye>::new(1.85).value(), 1.85);
    }

    #[test]
    fn from_scalar() {
        let d: DipoleMoment<f64, ElectronBohr> = DipoleMoment::from(3.0);
        assert_eq!(d.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(DipoleMoment::<f64, Debye>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = DipoleMoment::<f64, Debye>::new(1.85);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn electron_bohr_to_debye() {
        let d: DipoleMoment<f64, Debye> = DipoleMoment::<f64, ElectronBohr>::new(1.0).to();
        assert!((d.value() - 2.541_746_471_47).abs() < 1e-10);
    }

    #[test]
    fn debye_to_electron_bohr() {
        let au: DipoleMoment<f64, ElectronBohr> = DipoleMoment::<f64, Debye>::new(1.0).to();
        assert!((au.value() - 0.393_430_269_79).abs() < 1e-10);
    }

    #[test]
    fn debye_to_coulomb_meter() {
        let cm: DipoleMoment<f64, CoulombMeter> = DipoleMoment::<f64, Debye>::new(1.0).to();
        assert!((cm.value() - 3.335_640_951_981_521e-30).abs() < 1e-43);
    }

    #[test]
    fn coulomb_meter_to_debye() {
        let d: DipoleMoment<f64, Debye> =
            DipoleMoment::<f64, CoulombMeter>::new(3.335_640_951_981_521e-30).to();
        assert!((d.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn roundtrip_electron_bohr_coulomb_meter_electron_bohr() {
        let orig = DipoleMoment::<f64, ElectronBohr>::new(1.0);
        let back: DipoleMoment<f64, ElectronBohr> = orig.to::<CoulombMeter>().to();
        assert!((back.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = DipoleMoment::<f64, Debye>::new(1.0) + DipoleMoment::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut d = DipoleMoment::<f64, Debye>::new(1.0);
        d += DipoleMoment::new(0.5);
        assert_eq!(d.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = DipoleMoment::<f64, Debye>::new(3.0) - DipoleMoment::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut d = DipoleMoment::<f64, Debye>::new(3.0);
        d -= DipoleMoment::new(1.0);
        assert_eq!(d.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = DipoleMoment::<f64, Debye>::new(7.0) % DipoleMoment::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut d = DipoleMoment::<f64, Debye>::new(7.0);
        d %= DipoleMoment::new(3.0);
        assert_eq!(d.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-DipoleMoment::<f64, Debye>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((DipoleMoment::<f64, Debye>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut d = DipoleMoment::<f64, Debye>::new(2.0);
        d *= 3.0;
        assert_eq!(d.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((DipoleMoment::<f64, Debye>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut d = DipoleMoment::<f64, Debye>::new(6.0);
        d /= 2.0;
        assert_eq!(d.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = DipoleMoment::<f64, Debye>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut d = DipoleMoment::<f64, Debye>::new(7.0);
        d %= 3.0;
        assert_eq!(d.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = DipoleMoment::<f64, Debye>::new(6.0) / DipoleMoment::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = DipoleMoment::<f64, Debye>::new(1.0);
        assert_eq!(a, DipoleMoment::new(1.0));
        assert_ne!(a, DipoleMoment::new(2.0));
    }

    #[test]
    fn ord() {
        let a = DipoleMoment::<f64, Debye>::new(1.0);
        let b = DipoleMoment::<f64, Debye>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(DipoleMoment::<f64, Debye>::new(-3.0).abs().value(), 3.0);
        assert_eq!(DipoleMoment::<f64, Debye>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let d = DipoleMoment::<f64, Debye>::new(1.0);
        let nan = DipoleMoment::<f64, Debye>::new(f64::NAN);
        assert_eq!(d.min(nan).value(), 1.0);
        assert_eq!(nan.min(d).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let d = DipoleMoment::<f64, Debye>::new(1.0);
        let nan = DipoleMoment::<f64, Debye>::new(f64::NAN);
        assert_eq!(d.max(nan).value(), 1.0);
        assert_eq!(nan.max(d).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = DipoleMoment::<f64, Debye>::new(1.0);
        let hi = DipoleMoment::<f64, Debye>::new(2.0);
        assert_eq!(DipoleMoment::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(DipoleMoment::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(DipoleMoment::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = DipoleMoment::<f64, Debye>::new(2.0);
        let hi = DipoleMoment::<f64, Debye>::new(1.0);
        DipoleMoment::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(DipoleMoment::<f64, Debye>::new(3.0).signum(), 1.0);
        assert_eq!(DipoleMoment::<f64, Debye>::new(-3.0).signum(), -1.0);
    }

    #[test]
    fn copysign() {
        let d = DipoleMoment::<f64, Debye>::new(3.0);
        let sign = DipoleMoment::<f64, Debye>::new(-1.0);
        assert_eq!(d.copysign(sign).value(), -3.0);
        assert_eq!((-d).copysign(d).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(DipoleMoment::<f64, Debye>::new(2.7).floor().value(), 2.0);
        assert_eq!(DipoleMoment::<f64, Debye>::new(-2.3).floor().value(), -3.0);
    }

    #[test]
    fn ceil() {
        assert_eq!(DipoleMoment::<f64, Debye>::new(2.3).ceil().value(), 3.0);
        assert_eq!(DipoleMoment::<f64, Debye>::new(-2.7).ceil().value(), -2.0);
    }

    #[test]
    fn round() {
        assert_eq!(DipoleMoment::<f64, Debye>::new(2.5).round().value(), 3.0);
        assert_eq!(DipoleMoment::<f64, Debye>::new(-2.5).round().value(), -3.0);
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            DipoleMoment::<f64, Debye>::new(2.5)
                .round_ties_even()
                .value(),
            2.0
        );
        assert_eq!(
            DipoleMoment::<f64, Debye>::new(3.5)
                .round_ties_even()
                .value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(DipoleMoment::<f64, Debye>::new(2.7).trunc().value(), 2.0);
        assert_eq!(DipoleMoment::<f64, Debye>::new(-2.7).trunc().value(), -2.0);
    }

    #[test]
    fn fract() {
        assert!((DipoleMoment::<f64, Debye>::new(2.75).fract().value() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn div_euclid() {
        let q = DipoleMoment::<f64, Debye>::new(7.0).div_euclid(DipoleMoment::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = DipoleMoment::<f64, Debye>::new(-7.0).rem_euclid(DipoleMoment::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = DipoleMoment::<f64, Debye>::new(2.0).mul_add(3.0, DipoleMoment::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = DipoleMoment::<f64, Debye>::new(3.0).hypot(DipoleMoment::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(DipoleMoment::<f64, Debye>::new(f64::NAN).is_nan());
        assert!(!DipoleMoment::<f64, Debye>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(DipoleMoment::<f64, Debye>::new(f64::INFINITY).is_infinite());
        assert!(!DipoleMoment::<f64, Debye>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(DipoleMoment::<f64, Debye>::new(1.0).is_finite());
        assert!(!DipoleMoment::<f64, Debye>::new(f64::INFINITY).is_finite());
        assert!(!DipoleMoment::<f64, Debye>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(DipoleMoment::<f64, Debye>::new(1.0).is_sign_positive());
        assert!(!DipoleMoment::<f64, Debye>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(DipoleMoment::<f64, Debye>::new(-1.0).is_sign_negative());
        assert!(!DipoleMoment::<f64, Debye>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let v = [
            DipoleMoment::<f64, Debye>::new(1.0),
            DipoleMoment::new(2.0),
            DipoleMoment::new(3.0),
        ];
        let total: DipoleMoment<f64, Debye> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            DipoleMoment::<f64, Debye>::new(1.0),
            DipoleMoment::new(2.0),
            DipoleMoment::new(3.0),
        ];
        let total: DipoleMoment<f64, Debye> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: DipoleMoment<f64, Debye> = iter::empty::<DipoleMoment<f64, Debye>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(DipoleMoment::<f64, Debye>::new(1.85).to_string(), "1.85 D");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", DipoleMoment::<f64, Debye>::new(1.0)),
            "DipoleMoment(1.0)"
        );
    }

    #[test]
    fn f32_debye_to_electron_bohr() {
        let au: DipoleMoment<f32, ElectronBohr> = DipoleMoment::<f32, Debye>::new(1.0_f32).to();
        assert!((au.value() - 0.393_430_27_f32).abs() < 1e-4_f32);
    }

    #[test]
    fn f32_add() {
        let sum = DipoleMoment::<f32, Debye>::new(1.0_f32) + DipoleMoment::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

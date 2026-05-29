//! Frequency quantities and unit markers.
//!
//! The canonical unit is the **wavenumber** (cm⁻¹).
//!
//! | Type | Symbol | cm⁻¹ per unit |
//! |---|---|---|
//! | [`Wavenumber`] | cm⁻¹ | 1 |
//! | [`AtomicFrequency`] | a.u. | 219474.63136314 |
//! | [`Terahertz`] | THz | 33.3564095198152 |
//! | [`Hertz`] | Hz | 1 / 29979245800 |

use crate::units::quantity::define_quantity;

/// Marker trait for frequency units.
///
/// Implement this on a zero-sized type to define a new frequency unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of wavenumbers
/// (cm⁻¹) per one unit of `Self`.
pub trait FrequencyUnit {
    /// Wavenumbers (cm⁻¹) per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"cm⁻¹"`, `"THz"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A frequency parameterised by scalar type `V` and unit marker `U`.
    Frequency,
    FrequencyUnit
);

/// The wavenumber (cm⁻¹) — canonical frequency unit.
///
/// 1 cm⁻¹ = 29979245800 Hz (exact).
pub struct Wavenumber;

impl FrequencyUnit for Wavenumber {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "cm⁻¹";
}

/// The atomic unit of frequency (a.u.) (CODATA 2022).
///
/// 1 a.u. ≈ 219474.63136314 cm⁻¹.
pub struct AtomicFrequency;

impl FrequencyUnit for AtomicFrequency {
    const TO_CANONICAL: f64 = 219_474.631_363_14;
    const SYMBOL: &'static str = "a.u.";
}

/// The terahertz (THz).
///
/// 1 THz ≈ 33.3564095198152 cm⁻¹ (exact, c = 299792458 m s⁻¹).
pub struct Terahertz;

impl FrequencyUnit for Terahertz {
    const TO_CANONICAL: f64 = 33.356_409_519_815_2;
    const SYMBOL: &'static str = "THz";
}

/// The hertz (Hz) — SI base unit of frequency.
///
/// 1 Hz = 1 / 29979245800 cm⁻¹ (exact).
pub struct Hertz;

impl FrequencyUnit for Hertz {
    const TO_CANONICAL: f64 = 1.0 / 29_979_245_800.0;
    const SYMBOL: &'static str = "Hz";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Frequency::<f64, Wavenumber>::new(1000.0).value(), 1000.0);
    }

    #[test]
    fn from_scalar() {
        let f: Frequency<f64, Terahertz> = Frequency::from(3.0);
        assert_eq!(f.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Frequency::<f64, Wavenumber>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Frequency::<f64, Wavenumber>::new(1.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn terahertz_to_wavenumber() {
        let wn: Frequency<f64, Wavenumber> = Frequency::<f64, Terahertz>::new(1.0).to();
        assert!((wn.value() - 33.356_409_519_815_2).abs() < 1e-12);
    }

    #[test]
    fn wavenumber_to_terahertz() {
        let thz: Frequency<f64, Terahertz> =
            Frequency::<f64, Wavenumber>::new(33.356_409_519_815_2).to();
        assert!((thz.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn atomic_frequency_to_wavenumber() {
        let wn: Frequency<f64, Wavenumber> = Frequency::<f64, AtomicFrequency>::new(1.0).to();
        assert!((wn.value() - 219_474.631_363_14).abs() < 1e-6);
    }

    #[test]
    fn wavenumber_to_atomic_frequency() {
        let au: Frequency<f64, AtomicFrequency> =
            Frequency::<f64, Wavenumber>::new(219_474.631_363_14).to();
        assert!((au.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn hertz_to_wavenumber() {
        let wn: Frequency<f64, Wavenumber> = Frequency::<f64, Hertz>::new(29_979_245_800.0).to();
        assert!((wn.value() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn wavenumber_to_hertz() {
        let hz: Frequency<f64, Hertz> = Frequency::<f64, Wavenumber>::new(1.0).to();
        assert!((hz.value() - 29_979_245_800.0).abs() < 0.01);
    }

    #[test]
    fn roundtrip_terahertz_wavenumber_terahertz() {
        let orig = Frequency::<f64, Terahertz>::new(10.0);
        let back: Frequency<f64, Terahertz> = orig.to::<Wavenumber>().to();
        assert!((back.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Frequency::<f64, Wavenumber>::new(1000.0) + Frequency::new(500.0);
        assert_eq!(sum.value(), 1500.0);
    }

    #[test]
    fn add_assign() {
        let mut f = Frequency::<f64, Wavenumber>::new(1000.0);
        f += Frequency::new(500.0);
        assert_eq!(f.value(), 1500.0);
    }

    #[test]
    fn sub() {
        let diff = Frequency::<f64, Wavenumber>::new(3000.0) - Frequency::new(1000.0);
        assert_eq!(diff.value(), 2000.0);
    }

    #[test]
    fn sub_assign() {
        let mut f = Frequency::<f64, Wavenumber>::new(3000.0);
        f -= Frequency::new(1000.0);
        assert_eq!(f.value(), 2000.0);
    }

    #[test]
    fn neg() {
        assert_eq!(
            (-Frequency::<f64, Wavenumber>::new(1500.0)).value(),
            -1500.0
        );
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Frequency::<f64, Wavenumber>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut f = Frequency::<f64, Wavenumber>::new(2.0);
        f *= 3.0;
        assert_eq!(f.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Frequency::<f64, Wavenumber>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut f = Frequency::<f64, Wavenumber>::new(6.0);
        f /= 2.0;
        assert_eq!(f.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Frequency::<f64, Wavenumber>::new(6.0) / Frequency::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Frequency::<f64, Wavenumber>::new(1000.0);
        assert_eq!(a, Frequency::new(1000.0));
        assert_ne!(a, Frequency::new(2000.0));
    }

    #[test]
    fn ord() {
        let a = Frequency::<f64, Wavenumber>::new(1000.0);
        let b = Frequency::<f64, Wavenumber>::new(2000.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(
            Frequency::<f64, Wavenumber>::new(-1000.0).abs().value(),
            1000.0
        );
        assert_eq!(
            Frequency::<f64, Wavenumber>::new(1000.0).abs().value(),
            1000.0
        );
    }

    #[test]
    fn min_ignores_nan() {
        let f = Frequency::<f64, Wavenumber>::new(1000.0);
        let nan = Frequency::<f64, Wavenumber>::new(f64::NAN);
        assert_eq!(f.min(nan).value(), 1000.0);
        assert_eq!(nan.min(f).value(), 1000.0);
    }

    #[test]
    fn max_ignores_nan() {
        let f = Frequency::<f64, Wavenumber>::new(1000.0);
        let nan = Frequency::<f64, Wavenumber>::new(f64::NAN);
        assert_eq!(f.max(nan).value(), 1000.0);
        assert_eq!(nan.max(f).value(), 1000.0);
    }

    #[test]
    fn clamp() {
        let lo = Frequency::<f64, Wavenumber>::new(500.0);
        let hi = Frequency::<f64, Wavenumber>::new(3000.0);
        assert_eq!(Frequency::new(1000.0_f64).clamp(lo, hi).value(), 1000.0);
        assert_eq!(Frequency::new(100.0_f64).clamp(lo, hi).value(), 500.0);
        assert_eq!(Frequency::new(5000.0_f64).clamp(lo, hi).value(), 3000.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Frequency::<f64, Wavenumber>::new(3000.0);
        let hi = Frequency::<f64, Wavenumber>::new(500.0);
        Frequency::new(1000.0_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            Frequency::<f64, Wavenumber>::new(1000.0),
            Frequency::new(2000.0),
            Frequency::new(3000.0),
        ];
        let total: Frequency<f64, Wavenumber> = v.iter().copied().sum();
        assert_eq!(total.value(), 6000.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Frequency::<f64, Wavenumber>::new(1000.0),
            Frequency::new(2000.0),
            Frequency::new(3000.0),
        ];
        let total: Frequency<f64, Wavenumber> = v.iter().sum();
        assert_eq!(total.value(), 6000.0);
    }

    #[test]
    fn sum_empty() {
        let total: Frequency<f64, Wavenumber> = iter::empty::<Frequency<f64, Wavenumber>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Frequency::<f64, Terahertz>::new(1.5).to_string(), "1.5 THz");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Frequency::<f64, Wavenumber>::new(1.0)),
            "Frequency(1.0)"
        );
    }

    #[test]
    fn f32_terahertz_to_wavenumber() {
        let wn: Frequency<f32, Wavenumber> = Frequency::<f32, Terahertz>::new(1.0_f32).to();
        assert!((wn.value() - 33.356_41_f32).abs() < 1e-3_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Frequency::<f32, Wavenumber>::new(1000.0_f32) + Frequency::new(500.0_f32);
        assert_eq!(sum.value(), 1500.0_f32);
    }
}

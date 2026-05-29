//! Volume quantities and unit markers.
//!
//! The canonical unit is the **cubic ångström** (Å³).
//!
//! | Type | Symbol | Å³ per unit |
//! |---|---|---|
//! | [`CubicAngstrom`] | Å³ | 1 |
//! | [`CubicBohr`] | a₀³ | 0.148184711171 |
//! | [`CubicNanometer`] | nm³ | 1000 |
//! | [`CubicPicometer`] | pm³ | 1e-6 |
//! | [`CubicMeter`] | m³ | 1e30 |
//! | [`Liter`] | L | 1e27 |
//! | [`Milliliter`] | mL | 1e24 |

use crate::units::quantity::define_quantity;

/// Marker trait for volume units.
///
/// Implement this on a zero-sized type to define a new volume unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of cubic ångströms
/// per one unit of `Self`.
pub trait VolumeUnit {
    /// Cubic ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å³"`, `"nm³"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A volume parameterised by scalar type `V` and unit marker `U`.
    Volume,
    VolumeUnit
);

/// The cubic ångström (Å³) — canonical volume unit.
///
/// 1 Å³ = 1e-30 m³.
pub struct CubicAngstrom;

impl VolumeUnit for CubicAngstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å³";
}

/// The cubic bohr (a₀³) — atomic unit of volume (CODATA 2022, derived).
///
/// 1 a₀³ ≈ 0.148184711171 Å³.
pub struct CubicBohr;

impl VolumeUnit for CubicBohr {
    const TO_CANONICAL: f64 = 0.148_184_711_171;
    const SYMBOL: &'static str = "a₀³";
}

/// The cubic nanometre (nm³).
///
/// 1 nm³ = 1000 Å³.
pub struct CubicNanometer;

impl VolumeUnit for CubicNanometer {
    const TO_CANONICAL: f64 = 1_000.0;
    const SYMBOL: &'static str = "nm³";
}

/// The cubic picometre (pm³).
///
/// 1 pm³ = 1e-6 Å³.
pub struct CubicPicometer;

impl VolumeUnit for CubicPicometer {
    const TO_CANONICAL: f64 = 1e-6;
    const SYMBOL: &'static str = "pm³";
}

/// The cubic metre (m³) — SI base unit of volume.
///
/// 1 m³ = 1e30 Å³.
pub struct CubicMeter;

impl VolumeUnit for CubicMeter {
    const TO_CANONICAL: f64 = 1e30;
    const SYMBOL: &'static str = "m³";
}

/// The litre (L).
///
/// 1 L = 1e27 Å³.
pub struct Liter;

impl VolumeUnit for Liter {
    const TO_CANONICAL: f64 = 1e27;
    const SYMBOL: &'static str = "L";
}

/// The millilitre (mL).
///
/// 1 mL = 1e24 Å³.
pub struct Milliliter;

impl VolumeUnit for Milliliter {
    const TO_CANONICAL: f64 = 1e24;
    const SYMBOL: &'static str = "mL";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Volume::<f64, CubicAngstrom>::new(8.0).value(), 8.0);
    }

    #[test]
    fn from_scalar() {
        let v: Volume<f64, CubicNanometer> = Volume::from(3.0);
        assert_eq!(v.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Volume::<f64, CubicAngstrom>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Volume::<f64, CubicAngstrom>::new(2.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn cubic_angstrom_to_cubic_nanometer() {
        let nm3: Volume<f64, CubicNanometer> = Volume::<f64, CubicAngstrom>::new(1_000.0).to();
        assert!((nm3.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cubic_nanometer_to_cubic_angstrom() {
        let a3: Volume<f64, CubicAngstrom> = Volume::<f64, CubicNanometer>::new(1.0).to();
        assert!((a3.value() - 1_000.0).abs() < 1e-9);
    }

    #[test]
    fn cubic_bohr_to_cubic_angstrom() {
        let a3: Volume<f64, CubicAngstrom> = Volume::<f64, CubicBohr>::new(1.0).to();
        assert!((a3.value() - 0.148_184_711_171).abs() < 1e-12);
    }

    #[test]
    fn cubic_angstrom_to_cubic_picometer() {
        let pm3: Volume<f64, CubicPicometer> = Volume::<f64, CubicAngstrom>::new(1.0).to();
        assert!((pm3.value() - 1e6).abs() < 1e-6);
    }

    #[test]
    fn liter_to_cubic_angstrom() {
        let a3: Volume<f64, CubicAngstrom> = Volume::<f64, Liter>::new(1e-27).to();
        assert!((a3.value() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn roundtrip_cubic_nanometer_cubic_bohr_cubic_nanometer() {
        let orig = Volume::<f64, CubicNanometer>::new(0.5);
        let back: Volume<f64, CubicNanometer> = orig.to::<CubicBohr>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Volume::<f64, CubicAngstrom>::new(1.0) + Volume::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut v = Volume::<f64, CubicAngstrom>::new(1.0);
        v += Volume::new(0.5);
        assert_eq!(v.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Volume::<f64, CubicAngstrom>::new(3.0) - Volume::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut v = Volume::<f64, CubicAngstrom>::new(3.0);
        v -= Volume::new(1.0);
        assert_eq!(v.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Volume::<f64, CubicAngstrom>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Volume::<f64, CubicAngstrom>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut v = Volume::<f64, CubicAngstrom>::new(2.0);
        v *= 3.0;
        assert_eq!(v.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Volume::<f64, CubicAngstrom>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut v = Volume::<f64, CubicAngstrom>::new(6.0);
        v /= 2.0;
        assert_eq!(v.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Volume::<f64, CubicAngstrom>::new(6.0) / Volume::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Volume::<f64, CubicAngstrom>::new(1.0);
        assert_eq!(a, Volume::new(1.0));
        assert_ne!(a, Volume::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Volume::<f64, CubicAngstrom>::new(1.0);
        let b = Volume::<f64, CubicAngstrom>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Volume::<f64, CubicAngstrom>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Volume::<f64, CubicAngstrom>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let v = Volume::<f64, CubicAngstrom>::new(1.0);
        let nan = Volume::<f64, CubicAngstrom>::new(f64::NAN);
        assert_eq!(v.min(nan).value(), 1.0);
        assert_eq!(nan.min(v).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let v = Volume::<f64, CubicAngstrom>::new(1.0);
        let nan = Volume::<f64, CubicAngstrom>::new(f64::NAN);
        assert_eq!(v.max(nan).value(), 1.0);
        assert_eq!(nan.max(v).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Volume::<f64, CubicAngstrom>::new(1.0);
        let hi = Volume::<f64, CubicAngstrom>::new(2.0);
        assert_eq!(Volume::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Volume::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Volume::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Volume::<f64, CubicAngstrom>::new(2.0);
        let hi = Volume::<f64, CubicAngstrom>::new(1.0);
        Volume::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            Volume::<f64, CubicAngstrom>::new(1.0),
            Volume::new(2.0),
            Volume::new(3.0),
        ];
        let total: Volume<f64, CubicAngstrom> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Volume::<f64, CubicAngstrom>::new(1.0),
            Volume::new(2.0),
            Volume::new(3.0),
        ];
        let total: Volume<f64, CubicAngstrom> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Volume<f64, CubicAngstrom> = iter::empty::<Volume<f64, CubicAngstrom>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            Volume::<f64, CubicNanometer>::new(1.5).to_string(),
            "1.5 nm³"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Volume::<f64, CubicAngstrom>::new(1.0)),
            "Volume(1.0)"
        );
    }

    #[test]
    fn f32_cubic_angstrom_to_cubic_nanometer() {
        let nm3: Volume<f32, CubicNanometer> = Volume::<f32, CubicAngstrom>::new(1_000.0_f32).to();
        assert!((nm3.value() - 1.0_f32).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Volume::<f32, CubicAngstrom>::new(1.0_f32) + Volume::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

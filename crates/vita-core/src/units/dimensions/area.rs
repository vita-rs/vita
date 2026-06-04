//! Area quantities and unit markers.
//!
//! The canonical unit is the **square ångström** (Å²).
//!
//! | Type | Symbol | Å² per unit |
//! |---|---|---|
//! | [`SquareAngstrom`] | Å² | 1 |
//! | [`SquareBohr`] | a₀² | 0.280028520159 |
//! | [`SquareNanometer`] | nm² | 100 |
//! | [`SquarePicometer`] | pm² | 1e-4 |
//! | [`SquareMeter`] | m² | 1e20 |

use crate::units::quantity::define_quantity;

/// Marker trait for area units.
///
/// Implement this on a zero-sized type to define a new area unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of square
/// ångströms per one unit of `Self`.
pub trait AreaUnit {
    /// Square ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å²"`, `"nm²"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An area parameterized by scalar type `V` and unit marker `U`.
    Area,
    AreaUnit
);

/// The square ångström (Å²) — canonical area unit.
///
/// 1 Å² = 1e-20 m².
pub struct SquareAngstrom;

impl AreaUnit for SquareAngstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å²";
}

/// The square bohr (a₀²) — atomic unit of area (CODATA 2022, derived).
///
/// 1 a₀² ≈ 0.280028520159 Å².
pub struct SquareBohr;

impl AreaUnit for SquareBohr {
    const TO_CANONICAL: f64 = 0.280_028_520_159;
    const SYMBOL: &'static str = "a₀²";
}

/// The square nanometer (nm²).
///
/// 1 nm² = 100 Å².
pub struct SquareNanometer;

impl AreaUnit for SquareNanometer {
    const TO_CANONICAL: f64 = 100.0;
    const SYMBOL: &'static str = "nm²";
}

/// The square picometer (pm²).
///
/// 1 pm² = 1e-4 Å².
pub struct SquarePicometer;

impl AreaUnit for SquarePicometer {
    const TO_CANONICAL: f64 = 1e-4;
    const SYMBOL: &'static str = "pm²";
}

/// The square meter (m²) — SI base unit of area.
///
/// 1 m² = 1e20 Å².
pub struct SquareMeter;

impl AreaUnit for SquareMeter {
    const TO_CANONICAL: f64 = 1e20;
    const SYMBOL: &'static str = "m²";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Area::<f64, SquareAngstrom>::new(4.5).value(), 4.5);
    }

    #[test]
    fn from_scalar() {
        let a: Area<f64, SquareNanometer> = Area::from(3.0);
        assert_eq!(a.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Area::<f64, SquareAngstrom>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Area::<f64, SquareAngstrom>::new(1.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn square_angstrom_to_square_nanometer() {
        let nm2: Area<f64, SquareNanometer> = Area::<f64, SquareAngstrom>::new(100.0).to();
        assert!((nm2.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn square_nanometer_to_square_angstrom() {
        let a2: Area<f64, SquareAngstrom> = Area::<f64, SquareNanometer>::new(1.0).to();
        assert!((a2.value() - 100.0).abs() < 1e-12);
    }

    #[test]
    fn square_bohr_to_square_angstrom() {
        let a2: Area<f64, SquareAngstrom> = Area::<f64, SquareBohr>::new(1.0).to();
        assert!((a2.value() - 0.280_028_520_159).abs() < 1e-12);
    }

    #[test]
    fn square_angstrom_to_square_picometer() {
        let pm2: Area<f64, SquarePicometer> = Area::<f64, SquareAngstrom>::new(1.0).to();
        assert!((pm2.value() - 1e4).abs() < 1e-8);
    }

    #[test]
    fn square_angstrom_to_square_meter() {
        let m2: Area<f64, SquareMeter> = Area::<f64, SquareAngstrom>::new(1e20).to();
        assert!((m2.value() - 1.0).abs() < 1e-4);
    }

    #[test]
    fn roundtrip_square_nanometer_square_bohr_square_nanometer() {
        let orig = Area::<f64, SquareNanometer>::new(0.5);
        let back: Area<f64, SquareNanometer> = orig.to::<SquareBohr>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Area::<f64, SquareAngstrom>::new(1.0) + Area::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut a = Area::<f64, SquareAngstrom>::new(1.0);
        a += Area::new(0.5);
        assert_eq!(a.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Area::<f64, SquareAngstrom>::new(3.0) - Area::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut a = Area::<f64, SquareAngstrom>::new(3.0);
        a -= Area::new(1.0);
        assert_eq!(a.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = Area::<f64, SquareAngstrom>::new(7.0) % Area::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut a = Area::<f64, SquareAngstrom>::new(7.0);
        a %= Area::new(3.0);
        assert_eq!(a.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Area::<f64, SquareAngstrom>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Area::<f64, SquareAngstrom>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut a = Area::<f64, SquareAngstrom>::new(2.0);
        a *= 3.0;
        assert_eq!(a.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Area::<f64, SquareAngstrom>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut a = Area::<f64, SquareAngstrom>::new(6.0);
        a /= 2.0;
        assert_eq!(a.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Area::<f64, SquareAngstrom>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut a = Area::<f64, SquareAngstrom>::new(7.0);
        a %= 3.0;
        assert_eq!(a.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Area::<f64, SquareAngstrom>::new(6.0) / Area::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Area::<f64, SquareAngstrom>::new(1.0);
        assert_eq!(a, Area::new(1.0));
        assert_ne!(a, Area::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Area::<f64, SquareAngstrom>::new(1.0);
        let b = Area::<f64, SquareAngstrom>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Area::<f64, SquareAngstrom>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Area::<f64, SquareAngstrom>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let a = Area::<f64, SquareAngstrom>::new(1.0);
        let nan = Area::<f64, SquareAngstrom>::new(f64::NAN);
        assert_eq!(a.min(nan).value(), 1.0);
        assert_eq!(nan.min(a).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let a = Area::<f64, SquareAngstrom>::new(1.0);
        let nan = Area::<f64, SquareAngstrom>::new(f64::NAN);
        assert_eq!(a.max(nan).value(), 1.0);
        assert_eq!(nan.max(a).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Area::<f64, SquareAngstrom>::new(1.0);
        let hi = Area::<f64, SquareAngstrom>::new(2.0);
        assert_eq!(Area::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Area::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Area::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Area::<f64, SquareAngstrom>::new(2.0);
        let hi = Area::<f64, SquareAngstrom>::new(1.0);
        Area::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(Area::<f64, SquareAngstrom>::new(3.0).signum(), 1.0);
        assert_eq!(Area::<f64, SquareAngstrom>::new(-3.0).signum(), -1.0);
    }

    #[test]
    fn copysign() {
        let a = Area::<f64, SquareAngstrom>::new(3.0);
        let sign = Area::<f64, SquareAngstrom>::new(-1.0);
        assert_eq!(a.copysign(sign).value(), -3.0);
        assert_eq!((-a).copysign(a).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(Area::<f64, SquareAngstrom>::new(2.7).floor().value(), 2.0);
        assert_eq!(Area::<f64, SquareAngstrom>::new(-2.3).floor().value(), -3.0);
    }

    #[test]
    fn ceil() {
        assert_eq!(Area::<f64, SquareAngstrom>::new(2.3).ceil().value(), 3.0);
        assert_eq!(Area::<f64, SquareAngstrom>::new(-2.7).ceil().value(), -2.0);
    }

    #[test]
    fn round() {
        assert_eq!(Area::<f64, SquareAngstrom>::new(2.5).round().value(), 3.0);
        assert_eq!(Area::<f64, SquareAngstrom>::new(-2.5).round().value(), -3.0);
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            Area::<f64, SquareAngstrom>::new(2.5)
                .round_ties_even()
                .value(),
            2.0
        );
        assert_eq!(
            Area::<f64, SquareAngstrom>::new(3.5)
                .round_ties_even()
                .value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(Area::<f64, SquareAngstrom>::new(2.7).trunc().value(), 2.0);
        assert_eq!(Area::<f64, SquareAngstrom>::new(-2.7).trunc().value(), -2.0);
    }

    #[test]
    fn fract() {
        assert!((Area::<f64, SquareAngstrom>::new(2.75).fract().value() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn div_euclid() {
        let q = Area::<f64, SquareAngstrom>::new(7.0).div_euclid(Area::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Area::<f64, SquareAngstrom>::new(-7.0).rem_euclid(Area::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Area::<f64, SquareAngstrom>::new(2.0).mul_add(3.0, Area::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Area::<f64, SquareAngstrom>::new(3.0).hypot(Area::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Area::<f64, SquareAngstrom>::new(f64::NAN).is_nan());
        assert!(!Area::<f64, SquareAngstrom>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Area::<f64, SquareAngstrom>::new(f64::INFINITY).is_infinite());
        assert!(!Area::<f64, SquareAngstrom>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Area::<f64, SquareAngstrom>::new(1.0).is_finite());
        assert!(!Area::<f64, SquareAngstrom>::new(f64::INFINITY).is_finite());
        assert!(!Area::<f64, SquareAngstrom>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Area::<f64, SquareAngstrom>::new(1.0).is_sign_positive());
        assert!(!Area::<f64, SquareAngstrom>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Area::<f64, SquareAngstrom>::new(-1.0).is_sign_negative());
        assert!(!Area::<f64, SquareAngstrom>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let v = [
            Area::<f64, SquareAngstrom>::new(1.0),
            Area::new(2.0),
            Area::new(3.0),
        ];
        let total: Area<f64, SquareAngstrom> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Area::<f64, SquareAngstrom>::new(1.0),
            Area::new(2.0),
            Area::new(3.0),
        ];
        let total: Area<f64, SquareAngstrom> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Area<f64, SquareAngstrom> = iter::empty::<Area<f64, SquareAngstrom>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            Area::<f64, SquareNanometer>::new(1.5).to_string(),
            "1.5 nm²"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Area::<f64, SquareAngstrom>::new(1.0)),
            "Area(1.0)"
        );
    }

    #[test]
    fn f32_square_angstrom_to_square_nanometer() {
        let nm2: Area<f32, SquareNanometer> = Area::<f32, SquareAngstrom>::new(100.0_f32).to();
        assert!((nm2.value() - 1.0_f32).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Area::<f32, SquareAngstrom>::new(1.0_f32) + Area::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

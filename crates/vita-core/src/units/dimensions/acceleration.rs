//! Acceleration quantities and unit markers.
//!
//! The canonical unit is the **ångström per square picosecond** (Å ps⁻²).
//!
//! | Type | Symbol | Å ps⁻² per unit |
//! |---|---|---|
//! | [`AngstromPerSquarePicosecond`] | Å ps⁻² | 1 |
//! | [`NanometerPerSquarePicosecond`] | nm ps⁻² | 10 |
//! | [`AngstromPerSquareFemtosecond`] | Å fs⁻² | 1e6 |
//! | [`MeterPerSquareSecond`] | m s⁻² | 1e-14 |
//! | [`AtomicAcceleration`] | a₀ atu⁻² | 9.04421612109e8 |

use crate::units::quantity::define_quantity;

/// Marker trait for acceleration units.
///
/// Implement this on a zero-sized type to define a new acceleration unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of Å ps⁻²
/// per one unit of `Self`.
pub trait AccelerationUnit {
    /// Å ps⁻² per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å ps⁻²"`, `"nm ps⁻²"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An acceleration parameterized by scalar type `V` and unit marker `U`.
    Acceleration,
    AccelerationUnit
);

/// The ångström per square picosecond (Å ps⁻²) — canonical acceleration unit.
///
/// 1 Å ps⁻² = 1e14 m s⁻².
pub struct AngstromPerSquarePicosecond;

impl AccelerationUnit for AngstromPerSquarePicosecond {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å ps⁻²";
}

/// The nanometer per square picosecond (nm ps⁻²).
///
/// 1 nm ps⁻² = 10 Å ps⁻².
pub struct NanometerPerSquarePicosecond;

impl AccelerationUnit for NanometerPerSquarePicosecond {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "nm ps⁻²";
}

/// The ångström per square femtosecond (Å fs⁻²).
///
/// 1 Å fs⁻² = 1e6 Å ps⁻².
pub struct AngstromPerSquareFemtosecond;

impl AccelerationUnit for AngstromPerSquareFemtosecond {
    const TO_CANONICAL: f64 = 1e6;
    const SYMBOL: &'static str = "Å fs⁻²";
}

/// The meter per square second (m s⁻²) — SI unit of acceleration.
///
/// 1 m s⁻² = 1e-14 Å ps⁻².
pub struct MeterPerSquareSecond;

impl AccelerationUnit for MeterPerSquareSecond {
    const TO_CANONICAL: f64 = 1e-14;
    const SYMBOL: &'static str = "m s⁻²";
}

/// The atomic acceleration unit (a₀ atu⁻²) — atomic unit of acceleration (CODATA 2022, derived).
///
/// 1 a₀ atu⁻² ≈ 9.04421612109e8 Å ps⁻².
pub struct AtomicAcceleration;

impl AccelerationUnit for AtomicAcceleration {
    const TO_CANONICAL: f64 = 9.044_216_121_09e8;
    const SYMBOL: &'static str = "a₀ atu⁻²";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.52).value(),
            1.52
        );
    }

    #[test]
    fn from_scalar() {
        let a: Acceleration<f64, NanometerPerSquarePicosecond> = Acceleration::from(3.0);
        assert_eq!(a.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::default().value(),
            0.0_f64
        );
    }

    #[test]
    fn copy_and_clone() {
        let a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn angstrom_per_square_picosecond_to_nanometer_per_square_picosecond() {
        let nm: Acceleration<f64, NanometerPerSquarePicosecond> =
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(10.0).to();
        assert!((nm.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn nanometer_per_square_picosecond_to_angstrom_per_square_picosecond() {
        let a: Acceleration<f64, AngstromPerSquarePicosecond> =
            Acceleration::<f64, NanometerPerSquarePicosecond>::new(1.0).to();
        assert!((a.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn angstrom_per_square_picosecond_to_angstrom_per_square_femtosecond() {
        let afs: Acceleration<f64, AngstromPerSquareFemtosecond> =
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(1e6).to();
        assert!((afs.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn angstrom_per_square_picosecond_to_meter_per_square_second() {
        let ms2: Acceleration<f64, MeterPerSquareSecond> =
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(1e-14).to();
        assert!((ms2.value() - 1.0).abs() < 1e-4);
    }

    #[test]
    fn atomic_acceleration_to_angstrom_per_square_picosecond() {
        let a: Acceleration<f64, AngstromPerSquarePicosecond> =
            Acceleration::<f64, AtomicAcceleration>::new(1.0).to();
        assert!((a.value() - 9.044_216_121_09e8).abs() < 1.0);
    }

    #[test]
    fn roundtrip_nanometer_per_square_picosecond_atomic_acceleration_nanometer_per_square_picosecond()
     {
        let orig = Acceleration::<f64, NanometerPerSquarePicosecond>::new(0.5);
        let back: Acceleration<f64, NanometerPerSquarePicosecond> =
            orig.to::<AtomicAcceleration>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum =
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0) + Acceleration::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0);
        a += Acceleration::new(0.5);
        assert_eq!(a.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff =
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(3.0) - Acceleration::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(3.0);
        a -= Acceleration::new(1.0);
        assert_eq!(a.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = Acceleration::<f64, AngstromPerSquarePicosecond>::new(7.0) % Acceleration::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(7.0);
        a %= Acceleration::new(3.0);
        assert_eq!(a.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!(
            (-Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.5)).value(),
            -1.5
        );
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(
            (Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.0) * 3.0).value(),
            6.0
        );
    }

    #[test]
    fn mul_assign_scalar() {
        let mut a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.0);
        a *= 3.0;
        assert_eq!(a.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!(
            (Acceleration::<f64, AngstromPerSquarePicosecond>::new(6.0) / 2.0).value(),
            3.0
        );
    }

    #[test]
    fn div_assign_scalar() {
        let mut a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(6.0);
        a /= 2.0;
        assert_eq!(a.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Acceleration::<f64, AngstromPerSquarePicosecond>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(7.0);
        a %= 3.0;
        assert_eq!(a.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio =
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(6.0) / Acceleration::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0);
        assert_eq!(a, Acceleration::new(1.0));
        assert_ne!(a, Acceleration::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0);
        let b = Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(-3.0)
                .abs()
                .value(),
            3.0
        );
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(3.0)
                .abs()
                .value(),
            3.0
        );
    }

    #[test]
    fn min_ignores_nan() {
        let a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0);
        let nan = Acceleration::<f64, AngstromPerSquarePicosecond>::new(f64::NAN);
        assert_eq!(a.min(nan).value(), 1.0);
        assert_eq!(nan.min(a).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0);
        let nan = Acceleration::<f64, AngstromPerSquarePicosecond>::new(f64::NAN);
        assert_eq!(a.max(nan).value(), 1.0);
        assert_eq!(nan.max(a).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0);
        let hi = Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.0);
        assert_eq!(Acceleration::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Acceleration::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Acceleration::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.0);
        let hi = Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0);
        Acceleration::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(3.0).signum(),
            1.0
        );
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(-3.0).signum(),
            -1.0
        );
    }

    #[test]
    fn copysign() {
        let a = Acceleration::<f64, AngstromPerSquarePicosecond>::new(3.0);
        let sign = Acceleration::<f64, AngstromPerSquarePicosecond>::new(-1.0);
        assert_eq!(a.copysign(sign).value(), -3.0);
        assert_eq!((-a).copysign(a).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.7)
                .floor()
                .value(),
            2.0
        );
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(-2.3)
                .floor()
                .value(),
            -3.0
        );
    }

    #[test]
    fn ceil() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.3)
                .ceil()
                .value(),
            3.0
        );
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(-2.7)
                .ceil()
                .value(),
            -2.0
        );
    }

    #[test]
    fn round() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.5)
                .round()
                .value(),
            3.0
        );
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(-2.5)
                .round()
                .value(),
            -3.0
        );
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.5)
                .round_ties_even()
                .value(),
            2.0
        );
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(3.5)
                .round_ties_even()
                .value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.7)
                .trunc()
                .value(),
            2.0
        );
        assert_eq!(
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(-2.7)
                .trunc()
                .value(),
            -2.0
        );
    }

    #[test]
    fn fract() {
        assert!(
            (Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.75)
                .fract()
                .value()
                - 0.75)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn div_euclid() {
        let q = Acceleration::<f64, AngstromPerSquarePicosecond>::new(7.0)
            .div_euclid(Acceleration::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Acceleration::<f64, AngstromPerSquarePicosecond>::new(-7.0)
            .rem_euclid(Acceleration::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Acceleration::<f64, AngstromPerSquarePicosecond>::new(2.0)
            .mul_add(3.0, Acceleration::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Acceleration::<f64, AngstromPerSquarePicosecond>::new(3.0)
            .hypot(Acceleration::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Acceleration::<f64, AngstromPerSquarePicosecond>::new(f64::NAN).is_nan());
        assert!(!Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Acceleration::<f64, AngstromPerSquarePicosecond>::new(f64::INFINITY).is_infinite());
        assert!(!Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0).is_finite());
        assert!(!Acceleration::<f64, AngstromPerSquarePicosecond>::new(f64::INFINITY).is_finite());
        assert!(!Acceleration::<f64, AngstromPerSquarePicosecond>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0).is_sign_positive());
        assert!(!Acceleration::<f64, AngstromPerSquarePicosecond>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Acceleration::<f64, AngstromPerSquarePicosecond>::new(-1.0).is_sign_negative());
        assert!(!Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let accs = [
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0),
            Acceleration::new(2.0),
            Acceleration::new(3.0),
        ];
        let total: Acceleration<f64, AngstromPerSquarePicosecond> = accs.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let accs = [
            Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0),
            Acceleration::new(2.0),
            Acceleration::new(3.0),
        ];
        let total: Acceleration<f64, AngstromPerSquarePicosecond> = accs.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Acceleration<f64, AngstromPerSquarePicosecond> =
            iter::empty::<Acceleration<f64, AngstromPerSquarePicosecond>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            Acceleration::<f64, NanometerPerSquarePicosecond>::new(1.5).to_string(),
            "1.5 nm ps⁻²"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!(
                "{:?}",
                Acceleration::<f64, AngstromPerSquarePicosecond>::new(1.0)
            ),
            "Acceleration(1.0)"
        );
    }

    #[test]
    fn f32_angstrom_per_square_picosecond_to_nanometer_per_square_picosecond() {
        let nm: Acceleration<f32, NanometerPerSquarePicosecond> =
            Acceleration::<f32, AngstromPerSquarePicosecond>::new(10.0_f32).to();
        assert!((nm.value() - 1.0_f32).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Acceleration::<f32, AngstromPerSquarePicosecond>::new(1.0_f32)
            + Acceleration::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

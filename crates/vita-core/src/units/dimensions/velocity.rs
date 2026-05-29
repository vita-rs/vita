//! Velocity quantities and unit markers.
//!
//! The canonical unit is the **ångström per picosecond** (Å ps⁻¹).
//!
//! | Type | Symbol | Å ps⁻¹ per unit |
//! |---|---|---|
//! | [`AngstromPerPicosecond`] | Å ps⁻¹ | 1 |
//! | [`NanometerPerPicosecond`] | nm ps⁻¹ | 10 |
//! | [`AngstromPerFemtosecond`] | Å fs⁻¹ | 1000 |
//! | [`MeterPerSecond`] | m s⁻¹ | 0.01 |
//! | [`AtomicVelocity`] | a₀ atu⁻¹ | 2.18769126216e4 |

use crate::units::quantity::define_quantity;

/// Marker trait for velocity units.
///
/// Implement this on a zero-sized type to define a new velocity unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of Å ps⁻¹
/// per one unit of `Self`.
pub trait VelocityUnit {
    /// Å ps⁻¹ per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å ps⁻¹"`, `"nm ps⁻¹"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A velocity parameterised by scalar type `V` and unit marker `U`.
    Velocity,
    VelocityUnit
);

/// The ångström per picosecond (Å ps⁻¹) — canonical velocity unit.
///
/// 1 Å ps⁻¹ = 100 m s⁻¹.
pub struct AngstromPerPicosecond;

impl VelocityUnit for AngstromPerPicosecond {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å ps⁻¹";
}

/// The nanometre per picosecond (nm ps⁻¹).
///
/// 1 nm ps⁻¹ = 10 Å ps⁻¹.
pub struct NanometerPerPicosecond;

impl VelocityUnit for NanometerPerPicosecond {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "nm ps⁻¹";
}

/// The ångström per femtosecond (Å fs⁻¹).
///
/// 1 Å fs⁻¹ = 1000 Å ps⁻¹.
pub struct AngstromPerFemtosecond;

impl VelocityUnit for AngstromPerFemtosecond {
    const TO_CANONICAL: f64 = 1000.0;
    const SYMBOL: &'static str = "Å fs⁻¹";
}

/// The metre per second (m s⁻¹) — SI unit of velocity.
///
/// 1 m s⁻¹ = 0.01 Å ps⁻¹.
pub struct MeterPerSecond;

impl VelocityUnit for MeterPerSecond {
    const TO_CANONICAL: f64 = 0.01;
    const SYMBOL: &'static str = "m s⁻¹";
}

/// The atomic velocity unit (a₀ atu⁻¹) — atomic unit of velocity (CODATA 2022).
///
/// 1 a₀ atu⁻¹ ≈ 2.18769126216e4 Å ps⁻¹.
pub struct AtomicVelocity;

impl VelocityUnit for AtomicVelocity {
    const TO_CANONICAL: f64 = 2.187_691_262_16e4;
    const SYMBOL: &'static str = "a₀ atu⁻¹";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(1.52).value(),
            1.52
        );
    }

    #[test]
    fn from_scalar() {
        let v: Velocity<f64, NanometerPerPicosecond> = Velocity::from(3.0);
        assert_eq!(v.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::default().value(),
            0.0_f64
        );
    }

    #[test]
    fn copy_and_clone() {
        let a = Velocity::<f64, AngstromPerPicosecond>::new(2.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn angstrom_per_ps_to_nanometer_per_ps() {
        let nm: Velocity<f64, NanometerPerPicosecond> =
            Velocity::<f64, AngstromPerPicosecond>::new(10.0).to();
        assert!((nm.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn nanometer_per_ps_to_angstrom_per_ps() {
        let a: Velocity<f64, AngstromPerPicosecond> =
            Velocity::<f64, NanometerPerPicosecond>::new(1.0).to();
        assert!((a.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn angstrom_per_ps_to_angstrom_per_fs() {
        let afs: Velocity<f64, AngstromPerFemtosecond> =
            Velocity::<f64, AngstromPerPicosecond>::new(1000.0).to();
        assert!((afs.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn angstrom_per_ps_to_meter_per_second() {
        let ms: Velocity<f64, MeterPerSecond> =
            Velocity::<f64, AngstromPerPicosecond>::new(1.0).to();
        assert!((ms.value() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn atomic_velocity_to_angstrom_per_ps() {
        let a: Velocity<f64, AngstromPerPicosecond> =
            Velocity::<f64, AtomicVelocity>::new(1.0).to();
        assert!((a.value() - 2.187_691_262_16e4).abs() < 1e-7);
    }

    #[test]
    fn roundtrip_nanometer_per_ps_atomic_velocity_nanometer_per_ps() {
        let orig = Velocity::<f64, NanometerPerPicosecond>::new(0.5);
        let back: Velocity<f64, NanometerPerPicosecond> = orig.to::<AtomicVelocity>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Velocity::<f64, AngstromPerPicosecond>::new(1.0) + Velocity::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut v = Velocity::<f64, AngstromPerPicosecond>::new(1.0);
        v += Velocity::new(0.5);
        assert_eq!(v.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Velocity::<f64, AngstromPerPicosecond>::new(3.0) - Velocity::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut v = Velocity::<f64, AngstromPerPicosecond>::new(3.0);
        v -= Velocity::new(1.0);
        assert_eq!(v.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = Velocity::<f64, AngstromPerPicosecond>::new(7.0) % Velocity::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut v = Velocity::<f64, AngstromPerPicosecond>::new(7.0);
        v %= Velocity::new(3.0);
        assert_eq!(v.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!(
            (-Velocity::<f64, AngstromPerPicosecond>::new(1.5)).value(),
            -1.5
        );
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(
            (Velocity::<f64, AngstromPerPicosecond>::new(2.0) * 3.0).value(),
            6.0
        );
    }

    #[test]
    fn mul_assign_scalar() {
        let mut v = Velocity::<f64, AngstromPerPicosecond>::new(2.0);
        v *= 3.0;
        assert_eq!(v.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!(
            (Velocity::<f64, AngstromPerPicosecond>::new(6.0) / 2.0).value(),
            3.0
        );
    }

    #[test]
    fn div_assign_scalar() {
        let mut v = Velocity::<f64, AngstromPerPicosecond>::new(6.0);
        v /= 2.0;
        assert_eq!(v.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Velocity::<f64, AngstromPerPicosecond>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut v = Velocity::<f64, AngstromPerPicosecond>::new(7.0);
        v %= 3.0;
        assert_eq!(v.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Velocity::<f64, AngstromPerPicosecond>::new(6.0) / Velocity::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Velocity::<f64, AngstromPerPicosecond>::new(1.0);
        assert_eq!(a, Velocity::new(1.0));
        assert_ne!(a, Velocity::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Velocity::<f64, AngstromPerPicosecond>::new(1.0);
        let b = Velocity::<f64, AngstromPerPicosecond>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(-3.0)
                .abs()
                .value(),
            3.0
        );
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(3.0)
                .abs()
                .value(),
            3.0
        );
    }

    #[test]
    fn min_ignores_nan() {
        let v = Velocity::<f64, AngstromPerPicosecond>::new(1.0);
        let nan = Velocity::<f64, AngstromPerPicosecond>::new(f64::NAN);
        assert_eq!(v.min(nan).value(), 1.0);
        assert_eq!(nan.min(v).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let v = Velocity::<f64, AngstromPerPicosecond>::new(1.0);
        let nan = Velocity::<f64, AngstromPerPicosecond>::new(f64::NAN);
        assert_eq!(v.max(nan).value(), 1.0);
        assert_eq!(nan.max(v).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Velocity::<f64, AngstromPerPicosecond>::new(1.0);
        let hi = Velocity::<f64, AngstromPerPicosecond>::new(2.0);
        assert_eq!(Velocity::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Velocity::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Velocity::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Velocity::<f64, AngstromPerPicosecond>::new(2.0);
        let hi = Velocity::<f64, AngstromPerPicosecond>::new(1.0);
        Velocity::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(3.0).signum(),
            1.0
        );
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(-3.0).signum(),
            -1.0
        );
    }

    #[test]
    fn copysign() {
        let v = Velocity::<f64, AngstromPerPicosecond>::new(3.0);
        let sign = Velocity::<f64, AngstromPerPicosecond>::new(-1.0);
        assert_eq!(v.copysign(sign).value(), -3.0);
        assert_eq!((-v).copysign(v).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(2.7)
                .floor()
                .value(),
            2.0
        );
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(-2.3)
                .floor()
                .value(),
            -3.0
        );
    }

    #[test]
    fn ceil() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(2.3)
                .ceil()
                .value(),
            3.0
        );
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(-2.7)
                .ceil()
                .value(),
            -2.0
        );
    }

    #[test]
    fn round() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(2.5)
                .round()
                .value(),
            3.0
        );
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(-2.5)
                .round()
                .value(),
            -3.0
        );
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(2.5)
                .round_ties_even()
                .value(),
            2.0
        );
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(3.5)
                .round_ties_even()
                .value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(2.7)
                .trunc()
                .value(),
            2.0
        );
        assert_eq!(
            Velocity::<f64, AngstromPerPicosecond>::new(-2.7)
                .trunc()
                .value(),
            -2.0
        );
    }

    #[test]
    fn fract() {
        assert!(
            (Velocity::<f64, AngstromPerPicosecond>::new(2.75)
                .fract()
                .value()
                - 0.75)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn div_euclid() {
        let q = Velocity::<f64, AngstromPerPicosecond>::new(7.0).div_euclid(Velocity::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Velocity::<f64, AngstromPerPicosecond>::new(-7.0).rem_euclid(Velocity::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Velocity::<f64, AngstromPerPicosecond>::new(2.0).mul_add(3.0, Velocity::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Velocity::<f64, AngstromPerPicosecond>::new(3.0).hypot(Velocity::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Velocity::<f64, AngstromPerPicosecond>::new(f64::NAN).is_nan());
        assert!(!Velocity::<f64, AngstromPerPicosecond>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Velocity::<f64, AngstromPerPicosecond>::new(f64::INFINITY).is_infinite());
        assert!(!Velocity::<f64, AngstromPerPicosecond>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Velocity::<f64, AngstromPerPicosecond>::new(1.0).is_finite());
        assert!(!Velocity::<f64, AngstromPerPicosecond>::new(f64::INFINITY).is_finite());
        assert!(!Velocity::<f64, AngstromPerPicosecond>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Velocity::<f64, AngstromPerPicosecond>::new(1.0).is_sign_positive());
        assert!(!Velocity::<f64, AngstromPerPicosecond>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Velocity::<f64, AngstromPerPicosecond>::new(-1.0).is_sign_negative());
        assert!(!Velocity::<f64, AngstromPerPicosecond>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let vs = [
            Velocity::<f64, AngstromPerPicosecond>::new(1.0),
            Velocity::new(2.0),
            Velocity::new(3.0),
        ];
        let total: Velocity<f64, AngstromPerPicosecond> = vs.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let vs = [
            Velocity::<f64, AngstromPerPicosecond>::new(1.0),
            Velocity::new(2.0),
            Velocity::new(3.0),
        ];
        let total: Velocity<f64, AngstromPerPicosecond> = vs.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Velocity<f64, AngstromPerPicosecond> =
            iter::empty::<Velocity<f64, AngstromPerPicosecond>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            Velocity::<f64, NanometerPerPicosecond>::new(1.5).to_string(),
            "1.5 nm ps⁻¹"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Velocity::<f64, AngstromPerPicosecond>::new(1.0)),
            "Velocity(1.0)"
        );
    }

    #[test]
    fn f32_angstrom_per_ps_to_nanometer_per_ps() {
        let nm: Velocity<f32, NanometerPerPicosecond> =
            Velocity::<f32, AngstromPerPicosecond>::new(10.0_f32).to();
        assert!((nm.value() - 1.0_f32).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Velocity::<f32, AngstromPerPicosecond>::new(1.0_f32) + Velocity::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

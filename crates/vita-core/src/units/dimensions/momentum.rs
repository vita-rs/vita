//! Momentum quantities and unit markers.
//!
//! The canonical unit is the **dalton-ångström per picosecond** (Da Å ps⁻¹).
//!
//! | Type | Symbol | Da Å ps⁻¹ per unit |
//! |---|---|---|
//! | [`DaltonAngstromPerPicosecond`] | Da Å ps⁻¹ | 1 |
//! | [`DaltonNanometerPerPicosecond`] | Da nm ps⁻¹ | 10 |
//! | [`DaltonAngstromPerFemtosecond`] | Da Å fs⁻¹ | 1000 |
//! | [`AtomicMomentum`] | ℏ a₀⁻¹ | 12.001234736055 |
//! | [`KilogramMeterPerSecond`] | kg m s⁻¹ | 6.0221407537e24 |

use crate::units::quantity::define_quantity;

/// Marker trait for momentum units.
///
/// Implement this on a zero-sized type to define a new momentum unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of Da Å ps⁻¹
/// per one unit of `Self`.
pub trait MomentumUnit {
    /// Da Å ps⁻¹ per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Da Å ps⁻¹"`, `"ℏ a₀⁻¹"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A momentum parameterised by scalar type `V` and unit marker `U`.
    Momentum,
    MomentumUnit
);

/// The dalton-ångström per picosecond (Da Å ps⁻¹) — canonical momentum unit.
///
/// 1 Da Å ps⁻¹ ≈ 1.66053906892e-25 kg m s⁻¹ (CODATA 2022).
pub struct DaltonAngstromPerPicosecond;

impl MomentumUnit for DaltonAngstromPerPicosecond {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Da Å ps⁻¹";
}

/// The dalton-nanometre per picosecond (Da nm ps⁻¹).
///
/// 1 Da nm ps⁻¹ = 10 Da Å ps⁻¹.
pub struct DaltonNanometerPerPicosecond;

impl MomentumUnit for DaltonNanometerPerPicosecond {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "Da nm ps⁻¹";
}

/// The dalton-ångström per femtosecond (Da Å fs⁻¹).
///
/// 1 Da Å fs⁻¹ = 1000 Da Å ps⁻¹.
pub struct DaltonAngstromPerFemtosecond;

impl MomentumUnit for DaltonAngstromPerFemtosecond {
    const TO_CANONICAL: f64 = 1000.0;
    const SYMBOL: &'static str = "Da Å fs⁻¹";
}

/// The atomic momentum unit (ℏ a₀⁻¹) — atomic unit of momentum (CODATA 2022, derived).
///
/// 1 ℏ a₀⁻¹ ≈ 12.001234736055 Da Å ps⁻¹.
pub struct AtomicMomentum;

impl MomentumUnit for AtomicMomentum {
    const TO_CANONICAL: f64 = 12.001_234_736_055_5;
    const SYMBOL: &'static str = "ℏ a₀⁻¹";
}

/// The kilogram-metre per second (kg m s⁻¹) — SI unit of momentum (CODATA 2022, derived).
///
/// 1 kg m s⁻¹ ≈ 6.0221407537e24 Da Å ps⁻¹.
pub struct KilogramMeterPerSecond;

impl MomentumUnit for KilogramMeterPerSecond {
    const TO_CANONICAL: f64 = 6.022_140_7537e24;
    const SYMBOL: &'static str = "kg m s⁻¹";
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::format;
    use alloc::string::ToString;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(
            Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.52).value(),
            1.52
        );
    }

    #[test]
    fn from_scalar() {
        let p: Momentum<f64, DaltonNanometerPerPicosecond> = Momentum::from(3.0);
        assert_eq!(p.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(
            Momentum::<f64, DaltonAngstromPerPicosecond>::default().value(),
            0.0_f64
        );
    }

    #[test]
    fn copy_and_clone() {
        let a = Momentum::<f64, DaltonAngstromPerPicosecond>::new(2.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn dalton_angstrom_per_ps_to_dalton_nanometer_per_ps() {
        let nm: Momentum<f64, DaltonNanometerPerPicosecond> =
            Momentum::<f64, DaltonAngstromPerPicosecond>::new(10.0).to();
        assert!((nm.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dalton_nanometer_per_ps_to_dalton_angstrom_per_ps() {
        let a: Momentum<f64, DaltonAngstromPerPicosecond> =
            Momentum::<f64, DaltonNanometerPerPicosecond>::new(1.0).to();
        assert!((a.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn dalton_angstrom_per_ps_to_dalton_angstrom_per_fs() {
        let afs: Momentum<f64, DaltonAngstromPerFemtosecond> =
            Momentum::<f64, DaltonAngstromPerPicosecond>::new(1000.0).to();
        assert!((afs.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dalton_angstrom_per_ps_to_kilogram_meter_per_second() {
        let si: Momentum<f64, KilogramMeterPerSecond> =
            Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0).to();
        assert!((si.value() - 1.660_539_068_92e-25).abs() < 1e-35);
    }

    #[test]
    fn atomic_momentum_to_dalton_angstrom_per_ps() {
        let p: Momentum<f64, DaltonAngstromPerPicosecond> =
            Momentum::<f64, AtomicMomentum>::new(1.0).to();
        assert!((p.value() - 12.001_234_736_055_5).abs() < 1e-9);
    }

    #[test]
    fn roundtrip_dalton_nanometer_per_ps_atomic_momentum_dalton_nanometer_per_ps() {
        let orig = Momentum::<f64, DaltonNanometerPerPicosecond>::new(0.5);
        let back: Momentum<f64, DaltonNanometerPerPicosecond> = orig.to::<AtomicMomentum>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0) + Momentum::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut p = Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0);
        p += Momentum::new(0.5);
        assert_eq!(p.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Momentum::<f64, DaltonAngstromPerPicosecond>::new(3.0) - Momentum::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut p = Momentum::<f64, DaltonAngstromPerPicosecond>::new(3.0);
        p -= Momentum::new(1.0);
        assert_eq!(p.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!(
            (-Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.5)).value(),
            -1.5
        );
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(
            (Momentum::<f64, DaltonAngstromPerPicosecond>::new(2.0) * 3.0).value(),
            6.0
        );
    }

    #[test]
    fn mul_assign_scalar() {
        let mut p = Momentum::<f64, DaltonAngstromPerPicosecond>::new(2.0);
        p *= 3.0;
        assert_eq!(p.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!(
            (Momentum::<f64, DaltonAngstromPerPicosecond>::new(6.0) / 2.0).value(),
            3.0
        );
    }

    #[test]
    fn div_assign_scalar() {
        let mut p = Momentum::<f64, DaltonAngstromPerPicosecond>::new(6.0);
        p /= 2.0;
        assert_eq!(p.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Momentum::<f64, DaltonAngstromPerPicosecond>::new(6.0) / Momentum::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0);
        assert_eq!(a, Momentum::new(1.0));
        assert_ne!(a, Momentum::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0);
        let b = Momentum::<f64, DaltonAngstromPerPicosecond>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(
            Momentum::<f64, DaltonAngstromPerPicosecond>::new(-3.0)
                .abs()
                .value(),
            3.0
        );
        assert_eq!(
            Momentum::<f64, DaltonAngstromPerPicosecond>::new(3.0)
                .abs()
                .value(),
            3.0
        );
    }

    #[test]
    fn min_ignores_nan() {
        let p = Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0);
        let nan = Momentum::<f64, DaltonAngstromPerPicosecond>::new(f64::NAN);
        assert_eq!(p.min(nan).value(), 1.0);
        assert_eq!(nan.min(p).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let p = Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0);
        let nan = Momentum::<f64, DaltonAngstromPerPicosecond>::new(f64::NAN);
        assert_eq!(p.max(nan).value(), 1.0);
        assert_eq!(nan.max(p).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0);
        let hi = Momentum::<f64, DaltonAngstromPerPicosecond>::new(2.0);
        assert_eq!(Momentum::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Momentum::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Momentum::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Momentum::<f64, DaltonAngstromPerPicosecond>::new(2.0);
        let hi = Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0);
        Momentum::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0),
            Momentum::new(2.0),
            Momentum::new(3.0),
        ];
        let total: Momentum<f64, DaltonAngstromPerPicosecond> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0),
            Momentum::new(2.0),
            Momentum::new(3.0),
        ];
        let total: Momentum<f64, DaltonAngstromPerPicosecond> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Momentum<f64, DaltonAngstromPerPicosecond> =
            iter::empty::<Momentum<f64, DaltonAngstromPerPicosecond>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            Momentum::<f64, DaltonNanometerPerPicosecond>::new(1.5).to_string(),
            "1.5 Da nm ps⁻¹"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!(
                "{:?}",
                Momentum::<f64, DaltonAngstromPerPicosecond>::new(1.0)
            ),
            "Momentum(1.0)"
        );
    }

    #[test]
    fn f32_dalton_angstrom_per_ps_to_dalton_nanometer_per_ps() {
        let nm: Momentum<f32, DaltonNanometerPerPicosecond> =
            Momentum::<f32, DaltonAngstromPerPicosecond>::new(10.0_f32).to();
        assert!((nm.value() - 1.0_f32).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum =
            Momentum::<f32, DaltonAngstromPerPicosecond>::new(1.0_f32) + Momentum::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

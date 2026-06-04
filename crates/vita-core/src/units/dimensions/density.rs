//! Density quantities and unit markers.
//!
//! The canonical unit is the **gram per cubic centimeter** (g cm⁻³).
//!
//! | Type | Symbol | g cm⁻³ per unit |
//! |---|---|---|
//! | [`GramPerCubicCentimeter`] | g cm⁻³ | 1 |
//! | [`KilogramPerCubicMeter`] | kg m⁻³ | 0.001 |
//! | [`DaltonPerCubicAngstrom`] | Da Å⁻³ | 1.66053906892 |

use crate::units::quantity::define_quantity;

/// Marker trait for density units.
///
/// Implement this on a zero-sized type to define a new density unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of grams per
/// cubic centimeter per one unit of `Self`.
pub trait DensityUnit {
    /// Grams per cubic centimeter per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"g cm⁻³"`, `"kg m⁻³"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A density parameterized by scalar type `V` and unit marker `U`.
    Density,
    DensityUnit
);

/// The gram per cubic centimeter (g cm⁻³) — canonical density unit.
///
/// 1 g cm⁻³ = 1000 kg m⁻³.
pub struct GramPerCubicCentimeter;

impl DensityUnit for GramPerCubicCentimeter {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "g cm⁻³";
}

/// The kilogram per cubic meter (kg m⁻³) — SI unit of density.
///
/// 1 kg m⁻³ = 0.001 g cm⁻³.
pub struct KilogramPerCubicMeter;

impl DensityUnit for KilogramPerCubicMeter {
    const TO_CANONICAL: f64 = 0.001;
    const SYMBOL: &'static str = "kg m⁻³";
}

/// The dalton per cubic ångström (Da Å⁻³) — atomic unit of density (CODATA 2022, derived).
///
/// 1 Da Å⁻³ ≈ 1.66053906892 g cm⁻³.
pub struct DaltonPerCubicAngstrom;

impl DensityUnit for DaltonPerCubicAngstrom {
    const TO_CANONICAL: f64 = 1.660_539_068_92;
    const SYMBOL: &'static str = "Da Å⁻³";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(1.0).value(),
            1.0
        );
    }

    #[test]
    fn from_scalar() {
        let d: Density<f64, KilogramPerCubicMeter> = Density::from(1000.0);
        assert_eq!(d.value(), 1000.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::default().value(),
            0.0_f64
        );
    }

    #[test]
    fn copy_and_clone() {
        let a = Density::<f64, GramPerCubicCentimeter>::new(1.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn g_per_cm3_to_kg_per_m3() {
        let kg: Density<f64, KilogramPerCubicMeter> =
            Density::<f64, GramPerCubicCentimeter>::new(1.0).to();
        assert!((kg.value() - 1000.0).abs() < 1e-12);
    }

    #[test]
    fn kg_per_m3_to_g_per_cm3() {
        let g: Density<f64, GramPerCubicCentimeter> =
            Density::<f64, KilogramPerCubicMeter>::new(1000.0).to();
        assert!((g.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn da_per_a3_to_g_per_cm3() {
        let g: Density<f64, GramPerCubicCentimeter> =
            Density::<f64, DaltonPerCubicAngstrom>::new(1.0).to();
        assert!((g.value() - 1.660_539_068_92).abs() < 1e-12);
    }

    #[test]
    fn g_per_cm3_to_da_per_a3() {
        let da: Density<f64, DaltonPerCubicAngstrom> =
            Density::<f64, GramPerCubicCentimeter>::new(1.660_539_068_92).to();
        assert!((da.value() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn roundtrip_da_per_a3_kg_per_m3_da_per_a3() {
        let orig = Density::<f64, DaltonPerCubicAngstrom>::new(1.0);
        let back: Density<f64, DaltonPerCubicAngstrom> = orig.to::<KilogramPerCubicMeter>().to();
        assert!((back.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Density::<f64, GramPerCubicCentimeter>::new(1.0) + Density::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut d = Density::<f64, GramPerCubicCentimeter>::new(1.0);
        d += Density::new(0.5);
        assert_eq!(d.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Density::<f64, GramPerCubicCentimeter>::new(3.0) - Density::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut d = Density::<f64, GramPerCubicCentimeter>::new(3.0);
        d -= Density::new(1.0);
        assert_eq!(d.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = Density::<f64, GramPerCubicCentimeter>::new(7.0) % Density::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut d = Density::<f64, GramPerCubicCentimeter>::new(7.0);
        d %= Density::new(3.0);
        assert_eq!(d.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!(
            (-Density::<f64, GramPerCubicCentimeter>::new(1.5)).value(),
            -1.5
        );
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(
            (Density::<f64, GramPerCubicCentimeter>::new(2.0) * 3.0).value(),
            6.0
        );
    }

    #[test]
    fn mul_assign_scalar() {
        let mut d = Density::<f64, GramPerCubicCentimeter>::new(2.0);
        d *= 3.0;
        assert_eq!(d.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!(
            (Density::<f64, GramPerCubicCentimeter>::new(6.0) / 2.0).value(),
            3.0
        );
    }

    #[test]
    fn div_assign_scalar() {
        let mut d = Density::<f64, GramPerCubicCentimeter>::new(6.0);
        d /= 2.0;
        assert_eq!(d.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Density::<f64, GramPerCubicCentimeter>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut d = Density::<f64, GramPerCubicCentimeter>::new(7.0);
        d %= 3.0;
        assert_eq!(d.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Density::<f64, GramPerCubicCentimeter>::new(6.0) / Density::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Density::<f64, GramPerCubicCentimeter>::new(1.0);
        assert_eq!(a, Density::new(1.0));
        assert_ne!(a, Density::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Density::<f64, GramPerCubicCentimeter>::new(1.0);
        let b = Density::<f64, GramPerCubicCentimeter>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(-3.0)
                .abs()
                .value(),
            3.0
        );
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(3.0)
                .abs()
                .value(),
            3.0
        );
    }

    #[test]
    fn min_ignores_nan() {
        let d = Density::<f64, GramPerCubicCentimeter>::new(1.0);
        let nan = Density::<f64, GramPerCubicCentimeter>::new(f64::NAN);
        assert_eq!(d.min(nan).value(), 1.0);
        assert_eq!(nan.min(d).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let d = Density::<f64, GramPerCubicCentimeter>::new(1.0);
        let nan = Density::<f64, GramPerCubicCentimeter>::new(f64::NAN);
        assert_eq!(d.max(nan).value(), 1.0);
        assert_eq!(nan.max(d).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Density::<f64, GramPerCubicCentimeter>::new(1.0);
        let hi = Density::<f64, GramPerCubicCentimeter>::new(2.0);
        assert_eq!(Density::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Density::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Density::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Density::<f64, GramPerCubicCentimeter>::new(2.0);
        let hi = Density::<f64, GramPerCubicCentimeter>::new(1.0);
        Density::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(3.0).signum(),
            1.0
        );
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(-3.0).signum(),
            -1.0
        );
    }

    #[test]
    fn copysign() {
        let d = Density::<f64, GramPerCubicCentimeter>::new(3.0);
        let sign = Density::<f64, GramPerCubicCentimeter>::new(-1.0);
        assert_eq!(d.copysign(sign).value(), -3.0);
        assert_eq!((-d).copysign(d).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(2.7)
                .floor()
                .value(),
            2.0
        );
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(-2.3)
                .floor()
                .value(),
            -3.0
        );
    }

    #[test]
    fn ceil() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(2.3)
                .ceil()
                .value(),
            3.0
        );
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(-2.7)
                .ceil()
                .value(),
            -2.0
        );
    }

    #[test]
    fn round() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(2.5)
                .round()
                .value(),
            3.0
        );
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(-2.5)
                .round()
                .value(),
            -3.0
        );
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(2.5)
                .round_ties_even()
                .value(),
            2.0
        );
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(3.5)
                .round_ties_even()
                .value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(2.7)
                .trunc()
                .value(),
            2.0
        );
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(-2.7)
                .trunc()
                .value(),
            -2.0
        );
    }

    #[test]
    fn fract() {
        assert!(
            (Density::<f64, GramPerCubicCentimeter>::new(2.75)
                .fract()
                .value()
                - 0.75)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn div_euclid() {
        let q = Density::<f64, GramPerCubicCentimeter>::new(7.0).div_euclid(Density::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Density::<f64, GramPerCubicCentimeter>::new(-7.0).rem_euclid(Density::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Density::<f64, GramPerCubicCentimeter>::new(2.0).mul_add(3.0, Density::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Density::<f64, GramPerCubicCentimeter>::new(3.0).hypot(Density::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Density::<f64, GramPerCubicCentimeter>::new(f64::NAN).is_nan());
        assert!(!Density::<f64, GramPerCubicCentimeter>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Density::<f64, GramPerCubicCentimeter>::new(f64::INFINITY).is_infinite());
        assert!(!Density::<f64, GramPerCubicCentimeter>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Density::<f64, GramPerCubicCentimeter>::new(1.0).is_finite());
        assert!(!Density::<f64, GramPerCubicCentimeter>::new(f64::INFINITY).is_finite());
        assert!(!Density::<f64, GramPerCubicCentimeter>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Density::<f64, GramPerCubicCentimeter>::new(1.0).is_sign_positive());
        assert!(!Density::<f64, GramPerCubicCentimeter>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Density::<f64, GramPerCubicCentimeter>::new(-1.0).is_sign_negative());
        assert!(!Density::<f64, GramPerCubicCentimeter>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let v = [
            Density::<f64, GramPerCubicCentimeter>::new(1.0),
            Density::new(2.0),
            Density::new(3.0),
        ];
        let total: Density<f64, GramPerCubicCentimeter> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Density::<f64, GramPerCubicCentimeter>::new(1.0),
            Density::new(2.0),
            Density::new(3.0),
        ];
        let total: Density<f64, GramPerCubicCentimeter> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Density<f64, GramPerCubicCentimeter> =
            iter::empty::<Density<f64, GramPerCubicCentimeter>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            Density::<f64, GramPerCubicCentimeter>::new(1.5).to_string(),
            "1.5 g cm⁻³"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Density::<f64, GramPerCubicCentimeter>::new(1.0)),
            "Density(1.0)"
        );
    }

    #[test]
    fn f32_g_per_cm3_to_kg_per_m3() {
        let kg: Density<f32, KilogramPerCubicMeter> =
            Density::<f32, GramPerCubicCentimeter>::new(1.0_f32).to();
        assert!((kg.value() - 1000.0_f32).abs() < 1e-3_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Density::<f32, GramPerCubicCentimeter>::new(1.0_f32) + Density::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

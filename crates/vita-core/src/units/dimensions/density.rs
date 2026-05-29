//! Density quantities and unit markers.
//!
//! The canonical unit is the **gram per cubic centimetre** (g cm⁻³).
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
/// cubic centimetre per one unit of `Self`.
pub trait DensityUnit {
    /// Grams per cubic centimetre per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"g cm⁻³"`, `"kg m⁻³"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A density parameterised by scalar type `V` and unit marker `U`.
    Density,
    DensityUnit
);

/// The gram per cubic centimetre (g cm⁻³) — canonical density unit.
///
/// 1 g cm⁻³ = 1000 kg m⁻³.
pub struct GramPerCubicCentimeter;

impl DensityUnit for GramPerCubicCentimeter {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "g cm⁻³";
}

/// The kilogram per cubic metre (kg m⁻³) — SI unit of density.
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
        let c = a.clone();
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

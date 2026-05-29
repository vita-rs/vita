//! Concentration quantities and unit markers.
//!
//! The canonical unit is the **molar** (M).
//!
//! | Type | Symbol | M per unit |
//! |---|---|---|
//! | [`Molar`] | M | 1 |
//! | [`Millimolar`] | mM | 0.001 |
//! | [`Micromolar`] | μM | 1e-6 |
//! | [`Nanomolar`] | nM | 1e-9 |
//! | [`Picomolar`] | pM | 1e-12 |
//! | [`Femtomolar`] | fM | 1e-15 |
//! | [`MolePerCubicMeter`] | mol/m³ | 0.001 |

use crate::units::quantity::define_quantity;

/// Marker trait for concentration units.
///
/// Implement this on a zero-sized type to define a new concentration unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of mol L⁻¹
/// per one unit of `Self`.
pub trait ConcentrationUnit {
    /// mol L⁻¹ per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"M"`, `"μM"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A concentration parameterised by scalar type `V` and unit marker `U`.
    Concentration,
    ConcentrationUnit
);

/// The molar (M) — canonical concentration unit.
///
/// 1 M = 1 mol L⁻¹.
pub struct Molar;

impl ConcentrationUnit for Molar {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "M";
}

/// The millimolar (mM).
///
/// 1 mM = 0.001 M.
pub struct Millimolar;

impl ConcentrationUnit for Millimolar {
    const TO_CANONICAL: f64 = 0.001;
    const SYMBOL: &'static str = "mM";
}

/// The micromolar (μM).
///
/// 1 μM = 1e-6 M.
pub struct Micromolar;

impl ConcentrationUnit for Micromolar {
    const TO_CANONICAL: f64 = 1e-6;
    const SYMBOL: &'static str = "μM";
}

/// The nanomolar (nM).
///
/// 1 nM = 1e-9 M.
pub struct Nanomolar;

impl ConcentrationUnit for Nanomolar {
    const TO_CANONICAL: f64 = 1e-9;
    const SYMBOL: &'static str = "nM";
}

/// The picomolar (pM).
///
/// 1 pM = 1e-12 M.
pub struct Picomolar;

impl ConcentrationUnit for Picomolar {
    const TO_CANONICAL: f64 = 1e-12;
    const SYMBOL: &'static str = "pM";
}

/// The femtomolar (fM).
///
/// 1 fM = 1e-15 M.
pub struct Femtomolar;

impl ConcentrationUnit for Femtomolar {
    const TO_CANONICAL: f64 = 1e-15;
    const SYMBOL: &'static str = "fM";
}

/// The mole per cubic metre (mol/m³) — SI unit of concentration.
///
/// 1 mol/m³ = 0.001 M (exact, since 1 m³ = 1000 L by definition).
pub struct MolePerCubicMeter;

impl ConcentrationUnit for MolePerCubicMeter {
    const TO_CANONICAL: f64 = 0.001;
    const SYMBOL: &'static str = "mol/m³";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Concentration::<f64, Molar>::new(1.52).value(), 1.52);
    }

    #[test]
    fn from_scalar() {
        let c: Concentration<f64, Millimolar> = Concentration::from(3.0);
        assert_eq!(c.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Concentration::<f64, Molar>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Concentration::<f64, Molar>::new(2.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn molar_to_millimolar() {
        let mm: Concentration<f64, Millimolar> = Concentration::<f64, Molar>::new(1.0).to();
        assert!((mm.value() - 1000.0).abs() < 1e-12);
    }

    #[test]
    fn millimolar_to_molar() {
        let m: Concentration<f64, Molar> = Concentration::<f64, Millimolar>::new(1000.0).to();
        assert!((m.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn molar_to_nanomolar() {
        let nm: Concentration<f64, Nanomolar> = Concentration::<f64, Molar>::new(1.0).to();
        assert!((nm.value() - 1e9).abs() < 1e-3);
    }

    #[test]
    fn molar_to_mole_per_cubic_meter() {
        let si: Concentration<f64, MolePerCubicMeter> = Concentration::<f64, Molar>::new(1.0).to();
        assert!((si.value() - 1000.0).abs() < 1e-12);
    }

    #[test]
    fn roundtrip_molar_picomolar_molar() {
        let orig = Concentration::<f64, Molar>::new(0.5);
        let back: Concentration<f64, Molar> = orig.to::<Picomolar>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Concentration::<f64, Molar>::new(1.0) + Concentration::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut c = Concentration::<f64, Molar>::new(1.0);
        c += Concentration::new(0.5);
        assert_eq!(c.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Concentration::<f64, Molar>::new(3.0) - Concentration::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut c = Concentration::<f64, Molar>::new(3.0);
        c -= Concentration::new(1.0);
        assert_eq!(c.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Concentration::<f64, Molar>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Concentration::<f64, Molar>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut c = Concentration::<f64, Molar>::new(2.0);
        c *= 3.0;
        assert_eq!(c.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Concentration::<f64, Molar>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut c = Concentration::<f64, Molar>::new(6.0);
        c /= 2.0;
        assert_eq!(c.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Concentration::<f64, Molar>::new(6.0) / Concentration::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Concentration::<f64, Molar>::new(1.0);
        assert_eq!(a, Concentration::new(1.0));
        assert_ne!(a, Concentration::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Concentration::<f64, Molar>::new(1.0);
        let b = Concentration::<f64, Molar>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Concentration::<f64, Molar>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Concentration::<f64, Molar>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let c = Concentration::<f64, Molar>::new(1.0);
        let nan = Concentration::<f64, Molar>::new(f64::NAN);
        assert_eq!(c.min(nan).value(), 1.0);
        assert_eq!(nan.min(c).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let c = Concentration::<f64, Molar>::new(1.0);
        let nan = Concentration::<f64, Molar>::new(f64::NAN);
        assert_eq!(c.max(nan).value(), 1.0);
        assert_eq!(nan.max(c).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Concentration::<f64, Molar>::new(1.0);
        let hi = Concentration::<f64, Molar>::new(2.0);
        assert_eq!(Concentration::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Concentration::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Concentration::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Concentration::<f64, Molar>::new(2.0);
        let hi = Concentration::<f64, Molar>::new(1.0);
        Concentration::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            Concentration::<f64, Molar>::new(1.0),
            Concentration::new(2.0),
            Concentration::new(3.0),
        ];
        let total: Concentration<f64, Molar> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Concentration::<f64, Molar>::new(1.0),
            Concentration::new(2.0),
            Concentration::new(3.0),
        ];
        let total: Concentration<f64, Molar> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Concentration<f64, Molar> = iter::empty::<Concentration<f64, Molar>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            Concentration::<f64, Millimolar>::new(1.5).to_string(),
            "1.5 mM"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Concentration::<f64, Molar>::new(1.0)),
            "Concentration(1.0)"
        );
    }

    #[test]
    fn f32_molar_to_millimolar() {
        let mm: Concentration<f32, Millimolar> = Concentration::<f32, Molar>::new(1.0_f32).to();
        assert!((mm.value() - 1000.0_f32).abs() < 1e-3_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Concentration::<f32, Molar>::new(1.0_f32) + Concentration::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

//! Amount-of-substance quantities and unit markers.
//!
//! The canonical unit is the **mole** (mol).
//!
//! | Type | Symbol | mol per unit |
//! |---|---|---|
//! | [`Mole`] | mol | 1 |
//! | [`Millimole`] | mmol | 0.001 |
//! | [`Micromole`] | μmol | 1e-6 |
//! | [`Nanomole`] | nmol | 1e-9 |
//! | [`Picomole`] | pmol | 1e-12 |
//! | [`Femtomole`] | fmol | 1e-15 |
//! | [`Molecule`] | molecule | 1.660539067173847e-24 |

use crate::units::quantity::define_quantity;

/// Marker trait for amount-of-substance units.
///
/// Implement this on a zero-sized type to define a new amount-of-substance unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of moles
/// per one unit of `Self`.
pub trait AmountOfSubstanceUnit {
    /// Moles per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"mol"`, `"mmol"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An amount of substance parameterised by scalar type `V` and unit marker `U`.
    AmountOfSubstance,
    AmountOfSubstanceUnit
);

/// The mole (mol) — canonical amount-of-substance unit and SI base unit.
pub struct Mole;

impl AmountOfSubstanceUnit for Mole {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "mol";
}

/// The millimole (mmol).
///
/// 1 mmol = 0.001 mol.
pub struct Millimole;

impl AmountOfSubstanceUnit for Millimole {
    const TO_CANONICAL: f64 = 0.001;
    const SYMBOL: &'static str = "mmol";
}

/// The micromole (μmol).
///
/// 1 μmol = 1e-6 mol.
pub struct Micromole;

impl AmountOfSubstanceUnit for Micromole {
    const TO_CANONICAL: f64 = 1e-6;
    const SYMBOL: &'static str = "μmol";
}

/// The nanomole (nmol).
///
/// 1 nmol = 1e-9 mol.
pub struct Nanomole;

impl AmountOfSubstanceUnit for Nanomole {
    const TO_CANONICAL: f64 = 1e-9;
    const SYMBOL: &'static str = "nmol";
}

/// The picomole (pmol).
///
/// 1 pmol = 1e-12 mol.
pub struct Picomole;

impl AmountOfSubstanceUnit for Picomole {
    const TO_CANONICAL: f64 = 1e-12;
    const SYMBOL: &'static str = "pmol";
}

/// The femtomole (fmol).
///
/// 1 fmol = 1e-15 mol.
pub struct Femtomole;

impl AmountOfSubstanceUnit for Femtomole {
    const TO_CANONICAL: f64 = 1e-15;
    const SYMBOL: &'static str = "fmol";
}

/// The molecule — one indivisible entity (CODATA 2022).
///
/// 1 molecule = 1/Nₐ mol ≈ 1.660539067173847e-24 mol,
/// where Nₐ = 6.022 140 76 × 10²³ mol⁻¹ (exact, defined).
pub struct Molecule;

impl AmountOfSubstanceUnit for Molecule {
    const TO_CANONICAL: f64 = 1.660_539_067_173_846_6e-24;
    const SYMBOL: &'static str = "molecule";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(AmountOfSubstance::<f64, Mole>::new(1.52).value(), 1.52);
    }

    #[test]
    fn from_scalar() {
        let a: AmountOfSubstance<f64, Millimole> = AmountOfSubstance::from(3.0);
        assert_eq!(a.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(AmountOfSubstance::<f64, Mole>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = AmountOfSubstance::<f64, Mole>::new(2.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn mole_to_millimole() {
        let mmol: AmountOfSubstance<f64, Millimole> = AmountOfSubstance::<f64, Mole>::new(1.0).to();
        assert!((mmol.value() - 1000.0).abs() < 1e-12);
    }

    #[test]
    fn millimole_to_mole() {
        let mol: AmountOfSubstance<f64, Mole> =
            AmountOfSubstance::<f64, Millimole>::new(1000.0).to();
        assert!((mol.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mole_to_nanomole() {
        let nmol: AmountOfSubstance<f64, Nanomole> = AmountOfSubstance::<f64, Mole>::new(1.0).to();
        assert!((nmol.value() - 1e9).abs() < 1e-3);
    }

    #[test]
    fn mole_to_molecule() {
        let molecules: AmountOfSubstance<f64, Molecule> =
            AmountOfSubstance::<f64, Mole>::new(1.0).to();
        assert!((molecules.value() - 6.022_140_76e23).abs() < 1e9);
    }

    #[test]
    fn roundtrip_mole_molecule_mole() {
        let orig = AmountOfSubstance::<f64, Mole>::new(0.5);
        let back: AmountOfSubstance<f64, Mole> = orig.to::<Molecule>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = AmountOfSubstance::<f64, Mole>::new(1.0) + AmountOfSubstance::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut a = AmountOfSubstance::<f64, Mole>::new(1.0);
        a += AmountOfSubstance::new(0.5);
        assert_eq!(a.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = AmountOfSubstance::<f64, Mole>::new(3.0) - AmountOfSubstance::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut a = AmountOfSubstance::<f64, Mole>::new(3.0);
        a -= AmountOfSubstance::new(1.0);
        assert_eq!(a.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-AmountOfSubstance::<f64, Mole>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(
            (AmountOfSubstance::<f64, Mole>::new(2.0) * 3.0).value(),
            6.0
        );
    }

    #[test]
    fn mul_assign_scalar() {
        let mut a = AmountOfSubstance::<f64, Mole>::new(2.0);
        a *= 3.0;
        assert_eq!(a.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!(
            (AmountOfSubstance::<f64, Mole>::new(6.0) / 2.0).value(),
            3.0
        );
    }

    #[test]
    fn div_assign_scalar() {
        let mut a = AmountOfSubstance::<f64, Mole>::new(6.0);
        a /= 2.0;
        assert_eq!(a.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = AmountOfSubstance::<f64, Mole>::new(6.0) / AmountOfSubstance::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = AmountOfSubstance::<f64, Mole>::new(1.0);
        assert_eq!(a, AmountOfSubstance::new(1.0));
        assert_ne!(a, AmountOfSubstance::new(2.0));
    }

    #[test]
    fn ord() {
        let a = AmountOfSubstance::<f64, Mole>::new(1.0);
        let b = AmountOfSubstance::<f64, Mole>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(AmountOfSubstance::<f64, Mole>::new(-3.0).abs().value(), 3.0);
        assert_eq!(AmountOfSubstance::<f64, Mole>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let a = AmountOfSubstance::<f64, Mole>::new(1.0);
        let nan = AmountOfSubstance::<f64, Mole>::new(f64::NAN);
        assert_eq!(a.min(nan).value(), 1.0);
        assert_eq!(nan.min(a).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let a = AmountOfSubstance::<f64, Mole>::new(1.0);
        let nan = AmountOfSubstance::<f64, Mole>::new(f64::NAN);
        assert_eq!(a.max(nan).value(), 1.0);
        assert_eq!(nan.max(a).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = AmountOfSubstance::<f64, Mole>::new(1.0);
        let hi = AmountOfSubstance::<f64, Mole>::new(2.0);
        assert_eq!(AmountOfSubstance::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(AmountOfSubstance::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(AmountOfSubstance::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = AmountOfSubstance::<f64, Mole>::new(2.0);
        let hi = AmountOfSubstance::<f64, Mole>::new(1.0);
        AmountOfSubstance::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            AmountOfSubstance::<f64, Mole>::new(1.0),
            AmountOfSubstance::new(2.0),
            AmountOfSubstance::new(3.0),
        ];
        let total: AmountOfSubstance<f64, Mole> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            AmountOfSubstance::<f64, Mole>::new(1.0),
            AmountOfSubstance::new(2.0),
            AmountOfSubstance::new(3.0),
        ];
        let total: AmountOfSubstance<f64, Mole> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: AmountOfSubstance<f64, Mole> =
            iter::empty::<AmountOfSubstance<f64, Mole>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            AmountOfSubstance::<f64, Millimole>::new(1.5).to_string(),
            "1.5 mmol"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", AmountOfSubstance::<f64, Mole>::new(1.0)),
            "AmountOfSubstance(1.0)"
        );
    }

    #[test]
    fn f32_mole_to_millimole() {
        let mmol: AmountOfSubstance<f32, Millimole> =
            AmountOfSubstance::<f32, Mole>::new(1.0_f32).to();
        assert!((mmol.value() - 1000.0_f32).abs() < 1e-3_f32);
    }

    #[test]
    fn f32_add() {
        let sum = AmountOfSubstance::<f32, Mole>::new(1.0_f32) + AmountOfSubstance::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

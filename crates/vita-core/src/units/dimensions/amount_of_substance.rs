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

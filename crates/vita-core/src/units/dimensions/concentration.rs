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

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

//! Force-constant quantities and unit markers.
//!
//! The canonical unit is the **hartree per bohr squared** (Eₕ a₀⁻²).
//!
//! | Type | Symbol | Eₕ a₀⁻² per unit |
//! |---|---|---|
//! | [`HartreePerBohrSquared`] | Eₕ a₀⁻² | 1 |
//! | [`KcalPerMolPerAngstromSquared`] | kcal mol⁻¹ Å⁻² | 0.529177210544² / 627.509474063 |
//! | [`KjPerMolPerNanometerSquared`] | kJ mol⁻¹ nm⁻² | 0.0529177210544² / 2625.49963948 |
//! | [`NewtonPerMeter`] | N m⁻¹ | 5.29177210544e-11 / 8.2387235038e-8 |

use crate::units::quantity::define_quantity;

/// Marker trait for force-constant units.
///
/// Implement this on a zero-sized type to define a new force-constant unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of hartrees per
/// bohr squared (Eₕ a₀⁻²) per one unit of `Self`.
pub trait ForceConstantUnit {
    /// Eₕ a₀⁻² per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Eₕ a₀⁻²"`, `"N m⁻¹"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A force constant parameterised by scalar type `V` and unit marker `U`.
    ForceConstant,
    ForceConstantUnit
);

/// The hartree per bohr squared (Eₕ a₀⁻²) — canonical force-constant unit
/// (atomic unit of force constant, CODATA 2022).
///
/// 1 Eₕ a₀⁻² ≈ 1556.893105 N m⁻¹.
pub struct HartreePerBohrSquared;

impl ForceConstantUnit for HartreePerBohrSquared {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Eₕ a₀⁻²";
}

/// The kilocalorie per mole per ångström squared (kcal mol⁻¹ Å⁻²).
///
/// 1 kcal mol⁻¹ Å⁻² ≈ 0.529177210544² / 627.509474063 Eₕ a₀⁻² (CODATA 2022, derived).
pub struct KcalPerMolPerAngstromSquared;

impl ForceConstantUnit for KcalPerMolPerAngstromSquared {
    const TO_CANONICAL: f64 = 0.529_177_210_544 * 0.529_177_210_544 / 627.509_474_063;
    const SYMBOL: &'static str = "kcal mol⁻¹ Å⁻²";
}

/// The kilojoule per mole per nanometre squared (kJ mol⁻¹ nm⁻²).
///
/// 1 kJ mol⁻¹ nm⁻² ≈ 0.0529177210544² / 2625.49963948 Eₕ a₀⁻² (CODATA 2022, derived).
pub struct KjPerMolPerNanometerSquared;

impl ForceConstantUnit for KjPerMolPerNanometerSquared {
    const TO_CANONICAL: f64 = 0.052_917_721_054_4 * 0.052_917_721_054_4 / 2625.499_639_48;
    const SYMBOL: &'static str = "kJ mol⁻¹ nm⁻²";
}

/// The newton per metre (N m⁻¹) — SI unit of force constant (CODATA 2022).
///
/// 1 N m⁻¹ ≈ 5.29177210544e-11 / 8.2387235038e-8 Eₕ a₀⁻².
pub struct NewtonPerMeter;

impl ForceConstantUnit for NewtonPerMeter {
    const TO_CANONICAL: f64 = 5.291_772_105_44e-11 / 8.238_723_503_8e-8;
    const SYMBOL: &'static str = "N m⁻¹";
}

//! Force quantities and unit markers.
//!
//! The canonical unit is the **hartree per bohr** (Eₕ a₀⁻¹).
//!
//! | Type | Symbol | Eₕ a₀⁻¹ per unit |
//! |---|---|---|
//! | [`HartreePerBohr`] | Eₕ a₀⁻¹ | 1 |
//! | [`KcalPerMolPerAngstrom`] | kcal mol⁻¹ Å⁻¹ | 0.529177210544 / 627.509474063 |
//! | [`KjPerMolPerNanometer`] | kJ mol⁻¹ nm⁻¹ | 0.0529177210544 / 2625.49963948 |
//! | [`Newton`] | N | 1 / 8.2387235038e-8 |
//! | [`Piconewton`] | pN | 1e-12 / 8.2387235038e-8 |

use crate::units::quantity::define_quantity;

/// Marker trait for force units.
///
/// Implement this on a zero-sized type to define a new force unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of hartrees per
/// bohr (Eₕ a₀⁻¹) per one unit of `Self`.
pub trait ForceUnit {
    /// Eₕ a₀⁻¹ per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Eₕ a₀⁻¹"`, `"N"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A force parameterised by scalar type `V` and unit marker `U`.
    Force,
    ForceUnit
);

/// The hartree per bohr (Eₕ a₀⁻¹) — canonical force unit (atomic unit of force, CODATA 2022).
///
/// 1 Eₕ a₀⁻¹ ≈ 8.2387235038e-8 N.
pub struct HartreePerBohr;

impl ForceUnit for HartreePerBohr {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Eₕ a₀⁻¹";
}

/// The kilocalorie per mole per ångström (kcal mol⁻¹ Å⁻¹).
///
/// 1 kcal mol⁻¹ Å⁻¹ ≈ 0.529177210544 / 627.509474063 Eₕ a₀⁻¹ (CODATA 2022, derived).
pub struct KcalPerMolPerAngstrom;

impl ForceUnit for KcalPerMolPerAngstrom {
    const TO_CANONICAL: f64 = 0.529_177_210_544 / 627.509_474_063;
    const SYMBOL: &'static str = "kcal mol⁻¹ Å⁻¹";
}

/// The kilojoule per mole per nanometre (kJ mol⁻¹ nm⁻¹).
///
/// 1 kJ mol⁻¹ nm⁻¹ ≈ 0.0529177210544 / 2625.49963948 Eₕ a₀⁻¹ (CODATA 2022, derived).
pub struct KjPerMolPerNanometer;

impl ForceUnit for KjPerMolPerNanometer {
    const TO_CANONICAL: f64 = 0.052_917_721_054_4 / 2625.499_639_48;
    const SYMBOL: &'static str = "kJ mol⁻¹ nm⁻¹";
}

/// The newton (N) — SI unit of force (CODATA 2022).
///
/// 1 N ≈ 1 / 8.2387235038e-8 Eₕ a₀⁻¹.
pub struct Newton;

impl ForceUnit for Newton {
    const TO_CANONICAL: f64 = 1.0 / 8.238_723_503_8e-8;
    const SYMBOL: &'static str = "N";
}

/// The piconewton (pN) — used in single-molecule force spectroscopy (AFM, optical tweezers).
///
/// 1 pN ≈ 1e-12 / 8.2387235038e-8 Eₕ a₀⁻¹ (CODATA 2022).
pub struct Piconewton;

impl ForceUnit for Piconewton {
    const TO_CANONICAL: f64 = 1e-12 / 8.238_723_503_8e-8;
    const SYMBOL: &'static str = "pN";
}

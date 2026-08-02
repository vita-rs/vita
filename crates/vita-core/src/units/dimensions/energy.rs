//! Energy quantities and unit markers.
//!
//! The canonical unit is the **hartree** (Eₕ).
//!
//! | Type | Symbol | Eₕ per unit |
//! |---|---|---|
//! | [`Hartree`] | Eₕ | 1 |
//! | [`KilocaloriePerMole`] | kcal mol⁻¹ | 1 / 627.509474063 |
//! | [`KilojoulePerMole`] | kJ mol⁻¹ | 1 / 2625.49963948 |
//! | [`ElectronVolt`] | eV | 1 / 27.211386245981 |
//! | [`Wavenumber`] | cm⁻¹ | 1 / 219474.63136314 |
//! | [`MilliElectronVolt`] | meV | 1 / 27211.386245981 |
//! | [`Joule`] | J | 1 / 4.3597447222060e-18 |

use crate::units::quantity::define_quantity;

/// Marker trait for energy units.
///
/// Implement this on a zero-sized type to define a new energy unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of hartrees
/// per one unit of `Self`.
pub trait EnergyUnit {
    /// Hartrees per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Eₕ"`, `"eV"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An energy parameterized by scalar type `V` and unit marker `U`.
    Energy,
    EnergyUnit
);

/// The hartree (Eₕ) — canonical energy unit.
///
/// 1 Eₕ ≈ 4.3597447222060e-18 J (CODATA 2022).
pub struct Hartree;

impl EnergyUnit for Hartree {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Eₕ";
}

/// The kilocalorie per mole (kcal mol⁻¹).
///
/// 1 kcal mol⁻¹ ≈ 1 / 627.509474063 Eₕ (CODATA 2022, computed).
pub struct KilocaloriePerMole;

impl EnergyUnit for KilocaloriePerMole {
    const TO_CANONICAL: f64 = 1.0 / 627.509_474_063;
    const SYMBOL: &'static str = "kcal mol⁻¹";
}

/// The kilojoule per mole (kJ mol⁻¹).
///
/// 1 kJ mol⁻¹ ≈ 1 / 2625.49963948 Eₕ (CODATA 2022, computed).
pub struct KilojoulePerMole;

impl EnergyUnit for KilojoulePerMole {
    const TO_CANONICAL: f64 = 1.0 / 2_625.499_639_48;
    const SYMBOL: &'static str = "kJ mol⁻¹";
}

/// The electron volt (eV) (CODATA 2022).
///
/// 1 eV ≈ 1 / 27.211386245981 Eₕ.
pub struct ElectronVolt;

impl EnergyUnit for ElectronVolt {
    const TO_CANONICAL: f64 = 1.0 / 27.211_386_245_981;
    const SYMBOL: &'static str = "eV";
}

/// The wavenumber (cm⁻¹) — inverse centimeters (CODATA 2022).
///
/// 1 cm⁻¹ ≈ 1 / 219474.63136314 Eₕ.
pub struct Wavenumber;

impl EnergyUnit for Wavenumber {
    const TO_CANONICAL: f64 = 1.0 / 219_474.631_363_14;
    const SYMBOL: &'static str = "cm⁻¹";
}

/// The milli-electron volt (meV) (CODATA 2022).
///
/// 1 meV ≈ 1 / 27211.386245981 Eₕ.
pub struct MilliElectronVolt;

impl EnergyUnit for MilliElectronVolt {
    const TO_CANONICAL: f64 = 1.0 / 27_211.386_245_981;
    const SYMBOL: &'static str = "meV";
}

/// The joule (J) — SI derived unit of energy (CODATA 2022).
///
/// 1 J ≈ 1 / 4.3597447222060e-18 Eₕ.
pub struct Joule;

impl EnergyUnit for Joule {
    const TO_CANONICAL: f64 = 1.0 / 4.359_744_722_206_0e-18;
    const SYMBOL: &'static str = "J";
}

//! Mass quantities and unit markers.
//!
//! The canonical unit is the **dalton** (Da).
//!
//! | Type | Symbol | Da per unit |
//! |---|---|---|
//! | [`Dalton`] | Da | 1 |
//! | [`ElectronMass`] | mₑ | 5.485799090441e-4 |
//! | [`ProtonMass`] | mₚ | 1.0072764665789 |
//! | [`Gram`] | g | 6.0221407537e23 |
//! | [`Kilogram`] | kg | 6.0221407537e26 |

use crate::units::quantity::define_quantity;

/// Marker trait for mass units.
///
/// Implement this on a zero-sized type to define a new mass unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of daltons
/// per one unit of `Self`.
pub trait MassUnit {
    /// Daltons per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Da"`, `"kg"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A mass parameterized by scalar type `V` and unit marker `U`.
    Mass,
    MassUnit
);

/// The dalton (Da) — canonical mass unit.
///
/// 1 Da ≈ 1.66053906892e-27 kg (CODATA 2022).
pub struct Dalton;

impl MassUnit for Dalton {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Da";
}

/// The electron mass (mₑ) — atomic unit of mass (CODATA 2022).
///
/// 1 mₑ ≈ 5.485799090441e-4 Da.
pub struct ElectronMass;

impl MassUnit for ElectronMass {
    const TO_CANONICAL: f64 = 5.485_799_090_441e-4;
    const SYMBOL: &'static str = "mₑ";
}

/// The proton mass (mₚ) (CODATA 2022).
///
/// 1 mₚ ≈ 1.0072764665789 Da.
pub struct ProtonMass;

impl MassUnit for ProtonMass {
    const TO_CANONICAL: f64 = 1.007_276_466_578_9;
    const SYMBOL: &'static str = "mₚ";
}

/// The gram (g) (CODATA 2022).
///
/// 1 g ≈ 6.0221407537e23 Da.
pub struct Gram;

impl MassUnit for Gram {
    const TO_CANONICAL: f64 = 6.022_140_753_7e23;
    const SYMBOL: &'static str = "g";
}

/// The kilogram (kg) — SI base unit of mass (CODATA 2022).
///
/// 1 kg ≈ 6.0221407537e26 Da.
pub struct Kilogram;

impl MassUnit for Kilogram {
    const TO_CANONICAL: f64 = 6.022_140_753_7e26;
    const SYMBOL: &'static str = "kg";
}

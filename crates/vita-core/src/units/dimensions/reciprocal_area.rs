//! Reciprocal area quantities and unit markers.
//!
//! The canonical unit is the **reciprocal square ångström** (Å⁻²).
//!
//! | Type | Symbol | Å⁻² per unit |
//! |---|---|---|
//! | [`ReciprocalSquareAngstrom`] | Å⁻² | 1 |
//! | [`ReciprocalSquareBohr`] | a₀⁻² | 3.571064830938 |
//! | [`ReciprocalSquareNanometer`] | nm⁻² | 0.01 |
//! | [`ReciprocalSquarePicometer`] | pm⁻² | 1e4 |
//! | [`ReciprocalSquareMeter`] | m⁻² | 1e-20 |

use crate::units::quantity::define_quantity;

/// Marker trait for reciprocal area units.
///
/// Implement this on a zero-sized type to define a new reciprocal area unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of reciprocal
/// square ångströms per one unit of `Self`.
pub trait ReciprocalAreaUnit {
    /// Reciprocal square ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å⁻²"`, `"nm⁻²"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A reciprocal area parameterized by scalar type `V` and unit marker `U`.
    ReciprocalArea,
    ReciprocalAreaUnit
);

/// The reciprocal square ångström (Å⁻²) — canonical reciprocal area unit.
///
/// 1 Å⁻² = 1e20 m⁻².
pub struct ReciprocalSquareAngstrom;

impl ReciprocalAreaUnit for ReciprocalSquareAngstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å⁻²";
}

/// The reciprocal square bohr (a₀⁻²) — atomic unit of reciprocal area
/// (CODATA 2022, derived).
///
/// 1 a₀⁻² ≈ 3.571064830938 Å⁻².
pub struct ReciprocalSquareBohr;

impl ReciprocalAreaUnit for ReciprocalSquareBohr {
    const TO_CANONICAL: f64 = 3.571_064_830_938;
    const SYMBOL: &'static str = "a₀⁻²";
}

/// The reciprocal square nanometer (nm⁻²).
///
/// 1 nm⁻² = 0.01 Å⁻².
pub struct ReciprocalSquareNanometer;

impl ReciprocalAreaUnit for ReciprocalSquareNanometer {
    const TO_CANONICAL: f64 = 0.01;
    const SYMBOL: &'static str = "nm⁻²";
}

/// The reciprocal square picometer (pm⁻²).
///
/// 1 pm⁻² = 1e4 Å⁻².
pub struct ReciprocalSquarePicometer;

impl ReciprocalAreaUnit for ReciprocalSquarePicometer {
    const TO_CANONICAL: f64 = 1e4;
    const SYMBOL: &'static str = "pm⁻²";
}

/// The reciprocal square meter (m⁻²) — SI base unit of reciprocal area.
///
/// 1 m⁻² = 1e-20 Å⁻².
pub struct ReciprocalSquareMeter;

impl ReciprocalAreaUnit for ReciprocalSquareMeter {
    const TO_CANONICAL: f64 = 1e-20;
    const SYMBOL: &'static str = "m⁻²";
}

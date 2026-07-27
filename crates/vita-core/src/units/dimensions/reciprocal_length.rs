//! Reciprocal length quantities and unit markers.
//!
//! The canonical unit is the **reciprocal ångström** (Å⁻¹).
//!
//! | Type | Symbol | Å⁻¹ per unit |
//! |---|---|---|
//! | [`ReciprocalAngstrom`] | Å⁻¹ | 1 |
//! | [`ReciprocalBohr`] | a₀⁻¹ | 1.889726125908 |
//! | [`ReciprocalNanometer`] | nm⁻¹ | 0.1 |
//! | [`ReciprocalPicometer`] | pm⁻¹ | 100 |
//! | [`ReciprocalMeter`] | m⁻¹ | 1e-10 |

use crate::units::quantity::define_quantity;

/// Marker trait for reciprocal length units.
///
/// Implement this on a zero-sized type to define a new reciprocal length unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of reciprocal
/// ångströms per one unit of `Self`.
pub trait ReciprocalLengthUnit {
    /// Reciprocal ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å⁻¹"`, `"nm⁻¹"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A reciprocal length parameterized by scalar type `V` and unit marker
    /// `U`.
    ReciprocalLength,
    ReciprocalLengthUnit
);

/// The reciprocal ångström (Å⁻¹) — canonical reciprocal length unit.
///
/// 1 Å⁻¹ = 1e10 m⁻¹.
pub struct ReciprocalAngstrom;

impl ReciprocalLengthUnit for ReciprocalAngstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å⁻¹";
}

/// The reciprocal bohr (a₀⁻¹) — atomic unit of reciprocal length
/// (CODATA 2022, derived).
///
/// 1 a₀⁻¹ ≈ 1.889726125908 Å⁻¹.
pub struct ReciprocalBohr;

impl ReciprocalLengthUnit for ReciprocalBohr {
    const TO_CANONICAL: f64 = 1.889_726_125_908;
    const SYMBOL: &'static str = "a₀⁻¹";
}

/// The reciprocal nanometer (nm⁻¹).
///
/// 1 nm⁻¹ = 0.1 Å⁻¹.
pub struct ReciprocalNanometer;

impl ReciprocalLengthUnit for ReciprocalNanometer {
    const TO_CANONICAL: f64 = 0.1;
    const SYMBOL: &'static str = "nm⁻¹";
}

/// The reciprocal picometer (pm⁻¹).
///
/// 1 pm⁻¹ = 100 Å⁻¹.
pub struct ReciprocalPicometer;

impl ReciprocalLengthUnit for ReciprocalPicometer {
    const TO_CANONICAL: f64 = 100.0;
    const SYMBOL: &'static str = "pm⁻¹";
}

/// The reciprocal meter (m⁻¹) — SI base unit of reciprocal length.
///
/// 1 m⁻¹ = 1e-10 Å⁻¹.
pub struct ReciprocalMeter;

impl ReciprocalLengthUnit for ReciprocalMeter {
    const TO_CANONICAL: f64 = 1e-10;
    const SYMBOL: &'static str = "m⁻¹";
}

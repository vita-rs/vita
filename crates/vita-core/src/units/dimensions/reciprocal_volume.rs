//! Reciprocal volume quantities and unit markers.
//!
//! The canonical unit is the **reciprocal cubic ångström** (Å⁻³).
//!
//! | Type | Symbol | Å⁻³ per unit |
//! |---|---|---|
//! | [`ReciprocalCubicAngstrom`] | Å⁻³ | 1 |
//! | [`ReciprocalCubicBohr`] | a₀⁻³ | 6.748334508335 |
//! | [`ReciprocalCubicNanometer`] | nm⁻³ | 0.001 |
//! | [`ReciprocalCubicPicometer`] | pm⁻³ | 1e6 |
//! | [`ReciprocalCubicMeter`] | m⁻³ | 1e-30 |
//! | [`ReciprocalLiter`] | L⁻¹ | 1e-27 |
//! | [`ReciprocalMilliliter`] | mL⁻¹ | 1e-24 |

use crate::units::quantity::define_quantity;

/// Marker trait for reciprocal volume units.
///
/// Implement this on a zero-sized type to define a new reciprocal volume unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of reciprocal
/// cubic ångströms per one unit of `Self`.
pub trait ReciprocalVolumeUnit {
    /// Reciprocal cubic ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å⁻³"`, `"nm⁻³"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A reciprocal volume parameterized by scalar type `V` and unit marker
    /// `U`.
    ReciprocalVolume,
    ReciprocalVolumeUnit
);

/// The reciprocal cubic ångström (Å⁻³) — canonical reciprocal volume unit.
///
/// 1 Å⁻³ = 1e30 m⁻³.
pub struct ReciprocalCubicAngstrom;

impl ReciprocalVolumeUnit for ReciprocalCubicAngstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å⁻³";
}

/// The reciprocal cubic bohr (a₀⁻³) — atomic unit of reciprocal volume
/// (CODATA 2022, derived).
///
/// 1 a₀⁻³ ≈ 6.748334508335 Å⁻³.
pub struct ReciprocalCubicBohr;

impl ReciprocalVolumeUnit for ReciprocalCubicBohr {
    const TO_CANONICAL: f64 = 6.748_334_508_335;
    const SYMBOL: &'static str = "a₀⁻³";
}

/// The reciprocal cubic nanometer (nm⁻³).
///
/// 1 nm⁻³ = 0.001 Å⁻³.
pub struct ReciprocalCubicNanometer;

impl ReciprocalVolumeUnit for ReciprocalCubicNanometer {
    const TO_CANONICAL: f64 = 0.001;
    const SYMBOL: &'static str = "nm⁻³";
}

/// The reciprocal cubic picometer (pm⁻³).
///
/// 1 pm⁻³ = 1e6 Å⁻³.
pub struct ReciprocalCubicPicometer;

impl ReciprocalVolumeUnit for ReciprocalCubicPicometer {
    const TO_CANONICAL: f64 = 1e6;
    const SYMBOL: &'static str = "pm⁻³";
}

/// The reciprocal cubic meter (m⁻³) — SI derived unit of reciprocal volume.
///
/// 1 m⁻³ = 1e-30 Å⁻³.
pub struct ReciprocalCubicMeter;

impl ReciprocalVolumeUnit for ReciprocalCubicMeter {
    const TO_CANONICAL: f64 = 1e-30;
    const SYMBOL: &'static str = "m⁻³";
}

/// The reciprocal liter (L⁻¹).
///
/// 1 L⁻¹ = 1e-27 Å⁻³.
pub struct ReciprocalLiter;

impl ReciprocalVolumeUnit for ReciprocalLiter {
    const TO_CANONICAL: f64 = 1e-27;
    const SYMBOL: &'static str = "L⁻¹";
}

/// The reciprocal milliliter (mL⁻¹).
///
/// 1 mL⁻¹ = 1e-24 Å⁻³.
pub struct ReciprocalMilliliter;

impl ReciprocalVolumeUnit for ReciprocalMilliliter {
    const TO_CANONICAL: f64 = 1e-24;
    const SYMBOL: &'static str = "mL⁻¹";
}

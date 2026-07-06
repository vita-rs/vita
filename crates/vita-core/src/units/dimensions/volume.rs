//! Volume quantities and unit markers.
//!
//! The canonical unit is the **cubic ångström** (Å³).
//!
//! | Type | Symbol | Å³ per unit |
//! |---|---|---|
//! | [`CubicAngstrom`] | Å³ | 1 |
//! | [`CubicBohr`] | a₀³ | 0.148184711171 |
//! | [`CubicNanometer`] | nm³ | 1000 |
//! | [`CubicPicometer`] | pm³ | 1e-6 |
//! | [`CubicMeter`] | m³ | 1e30 |
//! | [`Liter`] | L | 1e27 |
//! | [`Milliliter`] | mL | 1e24 |

use crate::units::quantity::define_quantity;

/// Marker trait for volume units.
///
/// Implement this on a zero-sized type to define a new volume unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of cubic ångströms
/// per one unit of `Self`.
pub trait VolumeUnit {
    /// Cubic ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å³"`, `"nm³"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A volume parameterized by scalar type `V` and unit marker `U`.
    Volume,
    VolumeUnit
);

/// The cubic ångström (Å³) — canonical volume unit.
///
/// 1 Å³ = 1e-30 m³.
pub struct CubicAngstrom;

impl VolumeUnit for CubicAngstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å³";
}

/// The cubic bohr (a₀³) — atomic unit of volume (CODATA 2022, derived).
///
/// 1 a₀³ ≈ 0.148184711171 Å³.
pub struct CubicBohr;

impl VolumeUnit for CubicBohr {
    const TO_CANONICAL: f64 = 0.148_184_711_171;
    const SYMBOL: &'static str = "a₀³";
}

/// The cubic nanometer (nm³).
///
/// 1 nm³ = 1000 Å³.
pub struct CubicNanometer;

impl VolumeUnit for CubicNanometer {
    const TO_CANONICAL: f64 = 1_000.0;
    const SYMBOL: &'static str = "nm³";
}

/// The cubic picometer (pm³).
///
/// 1 pm³ = 1e-6 Å³.
pub struct CubicPicometer;

impl VolumeUnit for CubicPicometer {
    const TO_CANONICAL: f64 = 1e-6;
    const SYMBOL: &'static str = "pm³";
}

/// The cubic meter (m³) — SI base unit of volume.
///
/// 1 m³ = 1e30 Å³.
pub struct CubicMeter;

impl VolumeUnit for CubicMeter {
    const TO_CANONICAL: f64 = 1e30;
    const SYMBOL: &'static str = "m³";
}

/// The liter (L).
///
/// 1 L = 1e27 Å³.
pub struct Liter;

impl VolumeUnit for Liter {
    const TO_CANONICAL: f64 = 1e27;
    const SYMBOL: &'static str = "L";
}

/// The milliliter (mL).
///
/// 1 mL = 1e24 Å³.
pub struct Milliliter;

impl VolumeUnit for Milliliter {
    const TO_CANONICAL: f64 = 1e24;
    const SYMBOL: &'static str = "mL";
}

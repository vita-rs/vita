//! Area quantities and unit markers.
//!
//! The canonical unit is the **square ångström** (Å²).
//!
//! | Type | Symbol | Å² per unit |
//! |---|---|---|
//! | [`SquareAngstrom`] | Å² | 1 |
//! | [`SquareBohr`] | a₀² | 0.280028520159 |
//! | [`SquareNanometer`] | nm² | 100 |
//! | [`SquarePicometer`] | pm² | 1e-4 |
//! | [`SquareMeter`] | m² | 1e20 |

use crate::units::quantity::define_quantity;

/// Marker trait for area units.
///
/// Implement this on a zero-sized type to define a new area unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of square
/// ångströms per one unit of `Self`.
pub trait AreaUnit {
    /// Square ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å²"`, `"nm²"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An area parameterised by scalar type `V` and unit marker `U`.
    Area,
    AreaUnit
);

/// The square ångström (Å²) — canonical area unit.
///
/// 1 Å² = 1e-20 m².
pub struct SquareAngstrom;

impl AreaUnit for SquareAngstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å²";
}

/// The square bohr (a₀²) — atomic unit of area (CODATA 2022, derived).
///
/// 1 a₀² ≈ 0.280028520159 Å².
pub struct SquareBohr;

impl AreaUnit for SquareBohr {
    const TO_CANONICAL: f64 = 0.280_028_520_159;
    const SYMBOL: &'static str = "a₀²";
}

/// The square nanometre (nm²).
///
/// 1 nm² = 100 Å².
pub struct SquareNanometer;

impl AreaUnit for SquareNanometer {
    const TO_CANONICAL: f64 = 100.0;
    const SYMBOL: &'static str = "nm²";
}

/// The square picometre (pm²).
///
/// 1 pm² = 1e-4 Å².
pub struct SquarePicometer;

impl AreaUnit for SquarePicometer {
    const TO_CANONICAL: f64 = 1e-4;
    const SYMBOL: &'static str = "pm²";
}

/// The square metre (m²) — SI base unit of area.
///
/// 1 m² = 1e20 Å².
pub struct SquareMeter;

impl AreaUnit for SquareMeter {
    const TO_CANONICAL: f64 = 1e20;
    const SYMBOL: &'static str = "m²";
}

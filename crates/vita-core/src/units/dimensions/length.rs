//! Length quantities and unit markers.
//!
//! The canonical unit is the **ångström** (Å).
//!
//! | Type | Symbol | Å per unit |
//! |---|---|---|
//! | [`Angstrom`] | Å | 1 |
//! | [`Bohr`] | a₀ | 0.529177210544 |
//! | [`Nanometer`] | nm | 10 |
//! | [`Picometer`] | pm | 0.01 |
//! | [`Meter`] | m | 1e10 |

use crate::units::quantity::define_quantity;

/// Marker trait for length units.
///
/// Implement this on a zero-sized type to define a new length unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of ångströms
/// per one unit of `Self`.
pub trait LengthUnit {
    /// Ångströms per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å"`, `"nm"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A length parameterised by scalar type `V` and unit marker `U`.
    Length,
    LengthUnit
);

/// The ångström (Å) — canonical length unit.
///
/// 1 Å = 1e-10 m.
pub struct Angstrom;

impl LengthUnit for Angstrom {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å";
}

/// The bohr (a₀) — atomic unit of length (CODATA 2022).
///
/// 1 a₀ ≈ 0.529177210544 Å.
pub struct Bohr;

impl LengthUnit for Bohr {
    const TO_CANONICAL: f64 = 0.529_177_210_544;
    const SYMBOL: &'static str = "a₀";
}

/// The nanometre (nm).
///
/// 1 nm = 10 Å.
pub struct Nanometer;

impl LengthUnit for Nanometer {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "nm";
}

/// The picometre (pm).
///
/// 1 pm = 0.01 Å.
pub struct Picometer;

impl LengthUnit for Picometer {
    const TO_CANONICAL: f64 = 0.01;
    const SYMBOL: &'static str = "pm";
}

/// The metre (m) — SI base unit of length.
///
/// 1 m = 1e10 Å.
pub struct Meter;

impl LengthUnit for Meter {
    const TO_CANONICAL: f64 = 1e10;
    const SYMBOL: &'static str = "m";
}

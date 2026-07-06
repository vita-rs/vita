//! Pressure quantities and unit markers.
//!
//! The canonical unit is the **bar**.
//!
//! | Type | Symbol | bar per unit |
//! |---|---|---|
//! | [`Bar`] | bar | 1 |
//! | [`Atmosphere`] | atm | 1.01325 |
//! | [`Pascal`] | Pa | 1e-5 |
//! | [`Kilopascal`] | kPa | 0.01 |
//! | [`Megapascal`] | MPa | 10 |
//! | [`Gigapascal`] | GPa | 1e4 |

use crate::units::quantity::define_quantity;

/// Marker trait for pressure units.
///
/// Implement this on a zero-sized type to define a new pressure unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of bar
/// per one unit of `Self`.
pub trait PressureUnit {
    /// Bar per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"bar"`, `"Pa"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A pressure parameterized by scalar type `V` and unit marker `U`.
    Pressure,
    PressureUnit
);

/// The bar — canonical pressure unit.
///
/// 1 bar = 1e5 Pa (exact).
pub struct Bar;

impl PressureUnit for Bar {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "bar";
}

/// The standard atmosphere (atm).
///
/// 1 atm = 1.01325 bar (exact).
pub struct Atmosphere;

impl PressureUnit for Atmosphere {
    const TO_CANONICAL: f64 = 1.01325;
    const SYMBOL: &'static str = "atm";
}

/// The pascal (Pa) — SI base unit of pressure.
///
/// 1 Pa = 1e-5 bar (exact).
pub struct Pascal;

impl PressureUnit for Pascal {
    const TO_CANONICAL: f64 = 1e-5;
    const SYMBOL: &'static str = "Pa";
}

/// The kilopascal (kPa).
///
/// 1 kPa = 0.01 bar (exact).
pub struct Kilopascal;

impl PressureUnit for Kilopascal {
    const TO_CANONICAL: f64 = 0.01;
    const SYMBOL: &'static str = "kPa";
}

/// The megapascal (MPa).
///
/// 1 MPa = 10 bar (exact).
pub struct Megapascal;

impl PressureUnit for Megapascal {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "MPa";
}

/// The gigapascal (GPa) — common in high-pressure materials science.
///
/// 1 GPa = 1e4 bar (exact).
pub struct Gigapascal;

impl PressureUnit for Gigapascal {
    const TO_CANONICAL: f64 = 1e4;
    const SYMBOL: &'static str = "GPa";
}

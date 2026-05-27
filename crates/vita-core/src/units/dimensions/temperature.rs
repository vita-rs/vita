//! Temperature quantities and unit markers.
//!
//! The canonical unit is the **kelvin** (K).
//!
//! | Type | Symbol | K per unit |
//! |---|---|---|
//! | [`Kelvin`] | K | 1 |
//! | [`Rankine`] | °R | 5/9 |

use crate::units::quantity::define_quantity;

/// Marker trait for temperature units.
///
/// Implement this on a zero-sized type to define a new temperature unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of kelvins
/// per one unit of `Self`.
pub trait TemperatureUnit {
    /// Kelvins per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"K"`, `"°R"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A temperature parameterised by scalar type `V` and unit marker `U`.
    Temperature,
    TemperatureUnit
);

/// The kelvin (K) — canonical temperature unit.
///
/// SI base unit of thermodynamic temperature.
pub struct Kelvin;

impl TemperatureUnit for Kelvin {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "K";
}

/// The rankine (°R) — absolute temperature scale with Fahrenheit-sized degrees.
///
/// 1 °R = 5/9 K exactly.
pub struct Rankine;

impl TemperatureUnit for Rankine {
    const TO_CANONICAL: f64 = 5.0 / 9.0;
    const SYMBOL: &'static str = "°R";
}

//! Angle quantities and unit markers.
//!
//! The canonical unit is the **radian** (rad).
//!
//! | Type | Symbol | rad per unit |
//! |---|---|---|
//! | [`Radian`] | rad | 1 |
//! | [`Degree`] | ° | π / 180 |
//! | [`Milliradian`] | mrad | 1e-3 |
//! | [`Revolution`] | rev | 2π |

use crate::units::quantity::define_quantity;

/// Marker trait for angle units.
///
/// Implement this on a zero-sized type to define a new angle unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of radians
/// per one unit of `Self`.
pub trait AngleUnit {
    /// Radians per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"rad"`, `"°"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An angle parameterized by scalar type `V` and unit marker `U`.
    Angle,
    AngleUnit
);

/// The radian (rad) — canonical angle unit.
///
/// SI derived unit of plane angle.
pub struct Radian;

impl AngleUnit for Radian {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "rad";
}

/// The degree (°).
///
/// 1° = π / 180 rad.
pub struct Degree;

impl AngleUnit for Degree {
    const TO_CANONICAL: f64 = core::f64::consts::PI / 180.0;
    const SYMBOL: &'static str = "°";
}

/// The milliradian (mrad).
///
/// 1 mrad = 1e-3 rad.
pub struct Milliradian;

impl AngleUnit for Milliradian {
    const TO_CANONICAL: f64 = 1e-3;
    const SYMBOL: &'static str = "mrad";
}

/// The revolution (rev) — full turn.
///
/// 1 rev = 2π rad.
pub struct Revolution;

impl AngleUnit for Revolution {
    const TO_CANONICAL: f64 = 2.0 * core::f64::consts::PI;
    const SYMBOL: &'static str = "rev";
}

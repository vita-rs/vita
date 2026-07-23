//! Time quantities and unit markers.
//!
//! The canonical unit is the **picosecond** (ps).
//!
//! | Type | Symbol | ps per unit |
//! |---|---|---|
//! | [`Picosecond`] | ps | 1 |
//! | [`Femtosecond`] | fs | 0.001 |
//! | [`Nanosecond`] | ns | 1000 |
//! | [`Second`] | s | 1e12 |
//! | [`AtomicTime`] | atu | 2.4188843265864e-5 |

use crate::units::quantity::define_quantity;

/// Marker trait for time units.
///
/// Implement this on a zero-sized type to define a new time unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of picoseconds
/// per one unit of `Self`.
pub trait TimeUnit {
    /// Picoseconds per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"ps"`, `"ns"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A time parameterized by scalar type `V` and unit marker `U`.
    Time,
    TimeUnit
);

/// The picosecond (ps) — canonical time unit.
///
/// 1 ps = 1e-12 s.
pub struct Picosecond;

impl TimeUnit for Picosecond {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "ps";
}

/// The femtosecond (fs).
///
/// 1 fs = 0.001 ps.
pub struct Femtosecond;

impl TimeUnit for Femtosecond {
    const TO_CANONICAL: f64 = 0.001;
    const SYMBOL: &'static str = "fs";
}

/// The nanosecond (ns).
///
/// 1 ns = 1000 ps.
pub struct Nanosecond;

impl TimeUnit for Nanosecond {
    const TO_CANONICAL: f64 = 1000.0;
    const SYMBOL: &'static str = "ns";
}

/// The second (s) — SI base unit of time.
///
/// 1 s = 1e12 ps.
pub struct Second;

impl TimeUnit for Second {
    const TO_CANONICAL: f64 = 1e12;
    const SYMBOL: &'static str = "s";
}

/// The atomic time unit (atu) — atomic unit of time (CODATA 2022).
///
/// 1 atu ≈ 2.4188843265864e-5 ps.
pub struct AtomicTime;

impl TimeUnit for AtomicTime {
    const TO_CANONICAL: f64 = 2.418_884_326_586_4e-5;
    const SYMBOL: &'static str = "atu";
}

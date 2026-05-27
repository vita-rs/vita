//! Acceleration quantities and unit markers.
//!
//! The canonical unit is the **ångström per picosecond squared** (Å ps⁻²).
//!
//! | Type | Symbol | Å ps⁻² per unit |
//! |---|---|---|
//! | [`AngstromPerPicosecondSquared`] | Å ps⁻² | 1 |
//! | [`NanometerPerPicosecondSquared`] | nm ps⁻² | 10 |
//! | [`AngstromPerFemtosecondSquared`] | Å fs⁻² | 1e6 |
//! | [`MeterPerSecondSquared`] | m s⁻² | 1e-14 |
//! | [`AtomicAcceleration`] | a₀ atu⁻² | 9.04421612109e8 |

use crate::units::quantity::define_quantity;

/// Marker trait for acceleration units.
///
/// Implement this on a zero-sized type to define a new acceleration unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of Å ps⁻²
/// per one unit of `Self`.
pub trait AccelerationUnit {
    /// Å ps⁻² per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å ps⁻²"`, `"nm ps⁻²"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An acceleration parameterised by scalar type `V` and unit marker `U`.
    Acceleration,
    AccelerationUnit
);

/// The ångström per picosecond squared (Å ps⁻²) — canonical acceleration unit.
///
/// 1 Å ps⁻² = 1e14 m s⁻².
pub struct AngstromPerPicosecondSquared;

impl AccelerationUnit for AngstromPerPicosecondSquared {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å ps⁻²";
}

/// The nanometre per picosecond squared (nm ps⁻²).
///
/// 1 nm ps⁻² = 10 Å ps⁻².
pub struct NanometerPerPicosecondSquared;

impl AccelerationUnit for NanometerPerPicosecondSquared {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "nm ps⁻²";
}

/// The ångström per femtosecond squared (Å fs⁻²).
///
/// 1 Å fs⁻² = 1e6 Å ps⁻².
pub struct AngstromPerFemtosecondSquared;

impl AccelerationUnit for AngstromPerFemtosecondSquared {
    const TO_CANONICAL: f64 = 1e6;
    const SYMBOL: &'static str = "Å fs⁻²";
}

/// The metre per second squared (m s⁻²) — SI unit of acceleration.
///
/// 1 m s⁻² = 1e-14 Å ps⁻².
pub struct MeterPerSecondSquared;

impl AccelerationUnit for MeterPerSecondSquared {
    const TO_CANONICAL: f64 = 1e-14;
    const SYMBOL: &'static str = "m s⁻²";
}

/// The atomic acceleration unit (a₀ atu⁻²) — atomic unit of acceleration (CODATA 2022, derived).
///
/// 1 a₀ atu⁻² ≈ 9.04421612109e8 Å ps⁻².
pub struct AtomicAcceleration;

impl AccelerationUnit for AtomicAcceleration {
    const TO_CANONICAL: f64 = 9.044_216_121_09e8;
    const SYMBOL: &'static str = "a₀ atu⁻²";
}

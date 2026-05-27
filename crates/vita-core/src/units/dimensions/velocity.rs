//! Velocity quantities and unit markers.
//!
//! The canonical unit is the **ångström per picosecond** (Å ps⁻¹).
//!
//! | Type | Symbol | Å ps⁻¹ per unit |
//! |---|---|---|
//! | [`AngstromPerPicosecond`] | Å ps⁻¹ | 1 |
//! | [`NanometerPerPicosecond`] | nm ps⁻¹ | 10 |
//! | [`AngstromPerFemtosecond`] | Å fs⁻¹ | 1000 |
//! | [`MeterPerSecond`] | m s⁻¹ | 0.01 |
//! | [`AtomicVelocity`] | a₀ atu⁻¹ | 2.18769126216e4 |

use crate::units::quantity::define_quantity;

/// Marker trait for velocity units.
///
/// Implement this on a zero-sized type to define a new velocity unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of Å ps⁻¹
/// per one unit of `Self`.
pub trait VelocityUnit {
    /// Å ps⁻¹ per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Å ps⁻¹"`, `"nm ps⁻¹"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A velocity parameterised by scalar type `V` and unit marker `U`.
    Velocity,
    VelocityUnit
);

/// The ångström per picosecond (Å ps⁻¹) — canonical velocity unit.
///
/// 1 Å ps⁻¹ = 100 m s⁻¹.
pub struct AngstromPerPicosecond;

impl VelocityUnit for AngstromPerPicosecond {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Å ps⁻¹";
}

/// The nanometre per picosecond (nm ps⁻¹).
///
/// 1 nm ps⁻¹ = 10 Å ps⁻¹.
pub struct NanometerPerPicosecond;

impl VelocityUnit for NanometerPerPicosecond {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "nm ps⁻¹";
}

/// The ångström per femtosecond (Å fs⁻¹).
///
/// 1 Å fs⁻¹ = 1000 Å ps⁻¹.
pub struct AngstromPerFemtosecond;

impl VelocityUnit for AngstromPerFemtosecond {
    const TO_CANONICAL: f64 = 1000.0;
    const SYMBOL: &'static str = "Å fs⁻¹";
}

/// The metre per second (m s⁻¹) — SI unit of velocity.
///
/// 1 m s⁻¹ = 0.01 Å ps⁻¹.
pub struct MeterPerSecond;

impl VelocityUnit for MeterPerSecond {
    const TO_CANONICAL: f64 = 0.01;
    const SYMBOL: &'static str = "m s⁻¹";
}

/// The atomic velocity unit (a₀ atu⁻¹) — atomic unit of velocity (CODATA 2022).
///
/// 1 a₀ atu⁻¹ ≈ 2.18769126216e4 Å ps⁻¹.
pub struct AtomicVelocity;

impl VelocityUnit for AtomicVelocity {
    const TO_CANONICAL: f64 = 2.187_691_262_16e4;
    const SYMBOL: &'static str = "a₀ atu⁻¹";
}

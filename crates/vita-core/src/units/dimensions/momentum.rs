//! Momentum quantities and unit markers.
//!
//! The canonical unit is the **dalton-ångström per picosecond** (Da Å ps⁻¹).
//!
//! | Type | Symbol | Da Å ps⁻¹ per unit |
//! |---|---|---|
//! | [`DaltonAngstromPerPicosecond`] | Da Å ps⁻¹ | 1 |
//! | [`DaltonNanometerPerPicosecond`] | Da nm ps⁻¹ | 10 |
//! | [`DaltonAngstromPerFemtosecond`] | Da Å fs⁻¹ | 1000 |
//! | [`AtomicMomentum`] | ℏ a₀⁻¹ | 12.001234736055 |
//! | [`KilogramMeterPerSecond`] | kg m s⁻¹ | 6.0221407537e24 |

use crate::units::quantity::define_quantity;

/// Marker trait for momentum units.
///
/// Implement this on a zero-sized type to define a new momentum unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of Da Å ps⁻¹
/// per one unit of `Self`.
pub trait MomentumUnit {
    /// Da Å ps⁻¹ per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Da Å ps⁻¹"`, `"ℏ a₀⁻¹"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A momentum parameterized by scalar type `V` and unit marker `U`.
    Momentum,
    MomentumUnit
);

/// The dalton-ångström per picosecond (Da Å ps⁻¹) — canonical momentum unit.
///
/// 1 Da Å ps⁻¹ ≈ 1.66053906892e-25 kg m s⁻¹ (CODATA 2022).
pub struct DaltonAngstromPerPicosecond;

impl MomentumUnit for DaltonAngstromPerPicosecond {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Da Å ps⁻¹";
}

/// The dalton-nanometer per picosecond (Da nm ps⁻¹).
///
/// 1 Da nm ps⁻¹ = 10 Da Å ps⁻¹.
pub struct DaltonNanometerPerPicosecond;

impl MomentumUnit for DaltonNanometerPerPicosecond {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "Da nm ps⁻¹";
}

/// The dalton-ångström per femtosecond (Da Å fs⁻¹).
///
/// 1 Da Å fs⁻¹ = 1000 Da Å ps⁻¹.
pub struct DaltonAngstromPerFemtosecond;

impl MomentumUnit for DaltonAngstromPerFemtosecond {
    const TO_CANONICAL: f64 = 1000.0;
    const SYMBOL: &'static str = "Da Å fs⁻¹";
}

/// The atomic momentum unit (ℏ a₀⁻¹) — atomic unit of momentum (CODATA 2022, computed).
///
/// 1 ℏ a₀⁻¹ ≈ 12.001234736055 Da Å ps⁻¹.
pub struct AtomicMomentum;

impl MomentumUnit for AtomicMomentum {
    const TO_CANONICAL: f64 = 12.001_234_736_055_5;
    const SYMBOL: &'static str = "ℏ a₀⁻¹";
}

/// The kilogram-meter per second (kg m s⁻¹) — SI derived unit of momentum (CODATA 2022, computed).
///
/// 1 kg m s⁻¹ ≈ 6.0221407537e24 Da Å ps⁻¹.
pub struct KilogramMeterPerSecond;

impl MomentumUnit for KilogramMeterPerSecond {
    const TO_CANONICAL: f64 = 6.022_140_753_7e24;
    const SYMBOL: &'static str = "kg m s⁻¹";
}

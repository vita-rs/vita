//! Frequency quantities and unit markers.
//!
//! The canonical unit is the **wavenumber** (cm⁻¹).
//!
//! | Type | Symbol | cm⁻¹ per unit |
//! |---|---|---|
//! | [`Wavenumber`] | cm⁻¹ | 1 |
//! | [`AtomicFrequency`] | a.u. | 219474.63136314 |
//! | [`Terahertz`] | THz | 33.3564095198152 |
//! | [`Hertz`] | Hz | 1 / 29979245800 |

use crate::units::quantity::define_quantity;

/// Marker trait for frequency units.
///
/// Implement this on a zero-sized type to define a new frequency unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of wavenumbers
/// (cm⁻¹) per one unit of `Self`.
pub trait FrequencyUnit {
    /// Wavenumbers (cm⁻¹) per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"cm⁻¹"`, `"THz"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A frequency parameterized by scalar type `V` and unit marker `U`.
    Frequency,
    FrequencyUnit
);

/// The wavenumber (cm⁻¹) — canonical frequency unit.
///
/// 1 cm⁻¹ = 29979245800 Hz (exact).
pub struct Wavenumber;

impl FrequencyUnit for Wavenumber {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "cm⁻¹";
}

/// The atomic unit of frequency (a.u.) (CODATA 2022).
///
/// 1 a.u. ≈ 219474.63136314 cm⁻¹.
pub struct AtomicFrequency;

impl FrequencyUnit for AtomicFrequency {
    const TO_CANONICAL: f64 = 219_474.631_363_14;
    const SYMBOL: &'static str = "a.u.";
}

/// The terahertz (THz).
///
/// 1 THz ≈ 33.3564095198152 cm⁻¹ (exact, computed from c = 299 792 458 m s⁻¹).
pub struct Terahertz;

impl FrequencyUnit for Terahertz {
    const TO_CANONICAL: f64 = 33.356_409_519_815_2;
    const SYMBOL: &'static str = "THz";
}

/// The hertz (Hz) — SI derived unit of frequency.
///
/// 1 Hz = 1 / 29979245800 cm⁻¹ (exact).
pub struct Hertz;

impl FrequencyUnit for Hertz {
    const TO_CANONICAL: f64 = 1.0 / 29_979_245_800.0;
    const SYMBOL: &'static str = "Hz";
}

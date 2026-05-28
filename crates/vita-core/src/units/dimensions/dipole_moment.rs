//! Dipole-moment quantities and unit markers.
//!
//! The canonical unit is the **debye** (D).
//!
//! | Type | Symbol | D per unit |
//! |---|---|---|
//! | [`Debye`] | D | 1 |
//! | [`ElectronBohr`] | ea₀ | 8.4783536198e-30 / 3.335640951981521e-30 |
//! | [`CoulombMeter`] | C·m | 1 / 3.335640951981521e-30 |

use crate::units::quantity::define_quantity;

/// Marker trait for dipole-moment units.
///
/// Implement this on a zero-sized type to define a new dipole-moment unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of debyes
/// per one unit of `Self`.
pub trait DipoleMomentUnit {
    /// Debyes per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"D"`, `"ea₀"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A dipole moment parameterised by scalar type `V` and unit marker `U`.
    DipoleMoment,
    DipoleMomentUnit
);

/// The debye (D) — canonical dipole-moment unit.
///
/// 1 D = 3.335640951981521e-30 C·m (exact, derived from c = 299 792 458 m s⁻¹).
pub struct Debye;

impl DipoleMomentUnit for Debye {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "D";
}

/// The atomic unit of electric dipole moment (ea₀) (CODATA 2022).
///
/// 1 ea₀ = 8.4783536198e-30 C·m ≈ 2.5417464715 D.
pub struct ElectronBohr;

impl DipoleMomentUnit for ElectronBohr {
    const TO_CANONICAL: f64 = 8.478_353_619_8e-30 / 3.335_640_951_981_521e-30;
    const SYMBOL: &'static str = "ea₀";
}

/// The coulomb metre (C·m) — SI unit of electric dipole moment.
///
/// 1 C·m = 2.997924580e29 D (exact, derived from c = 299 792 458 m s⁻¹).
pub struct CoulombMeter;

impl DipoleMomentUnit for CoulombMeter {
    const TO_CANONICAL: f64 = 1.0 / 3.335_640_951_981_521e-30;
    const SYMBOL: &'static str = "C·m";
}

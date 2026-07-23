//! Charge quantities and unit markers.
//!
//! The canonical unit is the **elementary charge** (e).
//!
//! | Type | Symbol | e per unit |
//! |---|---|---|
//! | [`ElementaryCharge`] | e | 1 |
//! | [`Coulomb`] | C | 1 / 1.602176634e-19 |
//! | [`Nanocoulomb`] | nC | 1 / 1.602176634e-10 |
//! | [`Picocoulomb`] | pC | 1 / 1.602176634e-7 |

use crate::units::quantity::define_quantity;

/// Marker trait for charge units.
///
/// Implement this on a zero-sized type to define a new charge unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of elementary
/// charges per one unit of `Self`.
pub trait ChargeUnit {
    /// Elementary charges per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"e"`, `"C"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A charge parameterized by scalar type `V` and unit marker `U`.
    Charge,
    ChargeUnit
);

/// The elementary charge (e) — canonical charge unit.
///
/// 1 e = 1.602176634e-19 C (exact, CODATA 2022).
pub struct ElementaryCharge;

impl ChargeUnit for ElementaryCharge {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "e";
}

/// The coulomb (C) — SI base unit of electric charge.
///
/// 1 C ≈ 1 / 1.602176634e-19 e (CODATA 2022).
pub struct Coulomb;

impl ChargeUnit for Coulomb {
    const TO_CANONICAL: f64 = 1.0 / 1.602_176_634e-19;
    const SYMBOL: &'static str = "C";
}

/// The nanocoulomb (nC).
///
/// 1 nC = 1e-9 C.
pub struct Nanocoulomb;

impl ChargeUnit for Nanocoulomb {
    const TO_CANONICAL: f64 = 1.0 / 1.602_176_634e-10;
    const SYMBOL: &'static str = "nC";
}

/// The picocoulomb (pC).
///
/// 1 pC = 1e-12 C.
pub struct Picocoulomb;

impl ChargeUnit for Picocoulomb {
    const TO_CANONICAL: f64 = 1.0 / 1.602_176_634e-7;
    const SYMBOL: &'static str = "pC";
}

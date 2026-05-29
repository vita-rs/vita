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
    /// A charge parameterised by scalar type `V` and unit marker `U`.
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

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Charge::<f64, ElementaryCharge>::new(1.52).value(), 1.52);
    }

    #[test]
    fn from_scalar() {
        let q: Charge<f64, Coulomb> = Charge::from(3.0);
        assert_eq!(q.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Charge::<f64, ElementaryCharge>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Charge::<f64, ElementaryCharge>::new(2.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn elementary_charge_to_coulomb() {
        let c: Charge<f64, Coulomb> = Charge::<f64, ElementaryCharge>::new(1.0).to();
        assert!((c.value() - 1.602_176_634e-19).abs() < 1e-28);
    }

    #[test]
    fn coulomb_to_elementary_charge() {
        let q: Charge<f64, ElementaryCharge> = Charge::<f64, Coulomb>::new(1.602_176_634e-19).to();
        assert!((q.value() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn elementary_charge_to_picocoulomb() {
        let pc: Charge<f64, Picocoulomb> = Charge::<f64, ElementaryCharge>::new(1.0).to();
        assert!((pc.value() - 1.602_176_634e-7).abs() < 1e-17);
    }

    #[test]
    fn nanocoulomb_to_elementary_charge() {
        let q: Charge<f64, ElementaryCharge> =
            Charge::<f64, Nanocoulomb>::new(1.602_176_634e-10).to();
        assert!((q.value() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn roundtrip_coulomb_picocoulomb_coulomb() {
        let orig = Charge::<f64, Coulomb>::new(1.0);
        let back: Charge<f64, Coulomb> = orig.to::<Picocoulomb>().to();
        assert!((back.value() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn add() {
        let sum = Charge::<f64, ElementaryCharge>::new(1.0) + Charge::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut q = Charge::<f64, ElementaryCharge>::new(1.0);
        q += Charge::new(0.5);
        assert_eq!(q.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Charge::<f64, ElementaryCharge>::new(3.0) - Charge::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut q = Charge::<f64, ElementaryCharge>::new(3.0);
        q -= Charge::new(1.0);
        assert_eq!(q.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Charge::<f64, ElementaryCharge>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(
            (Charge::<f64, ElementaryCharge>::new(2.0) * 3.0).value(),
            6.0
        );
    }

    #[test]
    fn mul_assign_scalar() {
        let mut q = Charge::<f64, ElementaryCharge>::new(2.0);
        q *= 3.0;
        assert_eq!(q.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!(
            (Charge::<f64, ElementaryCharge>::new(6.0) / 2.0).value(),
            3.0
        );
    }

    #[test]
    fn div_assign_scalar() {
        let mut q = Charge::<f64, ElementaryCharge>::new(6.0);
        q /= 2.0;
        assert_eq!(q.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Charge::<f64, ElementaryCharge>::new(6.0) / Charge::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Charge::<f64, ElementaryCharge>::new(1.0);
        assert_eq!(a, Charge::new(1.0));
        assert_ne!(a, Charge::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Charge::<f64, ElementaryCharge>::new(1.0);
        let b = Charge::<f64, ElementaryCharge>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(
            Charge::<f64, ElementaryCharge>::new(-3.0).abs().value(),
            3.0
        );
        assert_eq!(Charge::<f64, ElementaryCharge>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let q = Charge::<f64, ElementaryCharge>::new(1.0);
        let nan = Charge::<f64, ElementaryCharge>::new(f64::NAN);
        assert_eq!(q.min(nan).value(), 1.0);
        assert_eq!(nan.min(q).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let q = Charge::<f64, ElementaryCharge>::new(1.0);
        let nan = Charge::<f64, ElementaryCharge>::new(f64::NAN);
        assert_eq!(q.max(nan).value(), 1.0);
        assert_eq!(nan.max(q).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Charge::<f64, ElementaryCharge>::new(1.0);
        let hi = Charge::<f64, ElementaryCharge>::new(2.0);
        assert_eq!(Charge::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Charge::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Charge::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Charge::<f64, ElementaryCharge>::new(2.0);
        let hi = Charge::<f64, ElementaryCharge>::new(1.0);
        Charge::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            Charge::<f64, ElementaryCharge>::new(1.0),
            Charge::new(2.0),
            Charge::new(3.0),
        ];
        let total: Charge<f64, ElementaryCharge> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Charge::<f64, ElementaryCharge>::new(1.0),
            Charge::new(2.0),
            Charge::new(3.0),
        ];
        let total: Charge<f64, ElementaryCharge> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Charge<f64, ElementaryCharge> =
            iter::empty::<Charge<f64, ElementaryCharge>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            Charge::<f64, ElementaryCharge>::new(1.5).to_string(),
            "1.5 e"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Charge::<f64, ElementaryCharge>::new(1.0)),
            "Charge(1.0)"
        );
    }

    #[test]
    fn f32_coulomb_to_elementary_charge() {
        let q: Charge<f32, ElementaryCharge> =
            Charge::<f32, Coulomb>::new(1.602_176_634e-19_f32).to();
        assert!((q.value() - 1.0_f32).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Charge::<f32, ElementaryCharge>::new(1.0_f32) + Charge::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

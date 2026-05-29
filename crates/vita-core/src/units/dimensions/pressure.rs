//! Pressure quantities and unit markers.
//!
//! The canonical unit is the **bar**.
//!
//! | Type | Symbol | bar per unit |
//! |---|---|---|
//! | [`Bar`] | bar | 1 |
//! | [`Atmosphere`] | atm | 1.01325 |
//! | [`Pascal`] | Pa | 1e-5 |
//! | [`Kilopascal`] | kPa | 0.01 |
//! | [`Megapascal`] | MPa | 10 |
//! | [`Gigapascal`] | GPa | 1e4 |

use crate::units::quantity::define_quantity;

/// Marker trait for pressure units.
///
/// Implement this on a zero-sized type to define a new pressure unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of bar
/// per one unit of `Self`.
pub trait PressureUnit {
    /// Bar per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"bar"`, `"Pa"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A pressure parameterised by scalar type `V` and unit marker `U`.
    Pressure,
    PressureUnit
);

/// The bar — canonical pressure unit.
///
/// 1 bar = 1e5 Pa (exact).
pub struct Bar;

impl PressureUnit for Bar {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "bar";
}

/// The standard atmosphere (atm).
///
/// 1 atm = 1.01325 bar (exact).
pub struct Atmosphere;

impl PressureUnit for Atmosphere {
    const TO_CANONICAL: f64 = 1.01325;
    const SYMBOL: &'static str = "atm";
}

/// The pascal (Pa) — SI base unit of pressure.
///
/// 1 Pa = 1e-5 bar (exact).
pub struct Pascal;

impl PressureUnit for Pascal {
    const TO_CANONICAL: f64 = 1e-5;
    const SYMBOL: &'static str = "Pa";
}

/// The kilopascal (kPa).
///
/// 1 kPa = 0.01 bar (exact).
pub struct Kilopascal;

impl PressureUnit for Kilopascal {
    const TO_CANONICAL: f64 = 0.01;
    const SYMBOL: &'static str = "kPa";
}

/// The megapascal (MPa).
///
/// 1 MPa = 10 bar (exact).
pub struct Megapascal;

impl PressureUnit for Megapascal {
    const TO_CANONICAL: f64 = 10.0;
    const SYMBOL: &'static str = "MPa";
}

/// The gigapascal (GPa) — common in high-pressure materials science.
///
/// 1 GPa = 1e4 bar (exact).
pub struct Gigapascal;

impl PressureUnit for Gigapascal {
    const TO_CANONICAL: f64 = 1e4;
    const SYMBOL: &'static str = "GPa";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Pressure::<f64, Bar>::new(1.013_25).value(), 1.013_25);
    }

    #[test]
    fn from_scalar() {
        let p: Pressure<f64, Atmosphere> = Pressure::from(3.0);
        assert_eq!(p.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Pressure::<f64, Bar>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Pressure::<f64, Bar>::new(2.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn atmosphere_to_bar() {
        let p: Pressure<f64, Bar> = Pressure::<f64, Atmosphere>::new(1.0).to();
        assert!((p.value() - 1.01325).abs() < 1e-12);
    }

    #[test]
    fn bar_to_pascal() {
        let p: Pressure<f64, Pascal> = Pressure::<f64, Bar>::new(1.0).to();
        assert!((p.value() - 1e5).abs() < 1e-8);
    }

    #[test]
    fn pascal_to_bar() {
        let p: Pressure<f64, Bar> = Pressure::<f64, Pascal>::new(1e5).to();
        assert!((p.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn gigapascal_to_bar() {
        let p: Pressure<f64, Bar> = Pressure::<f64, Gigapascal>::new(1.0).to();
        assert!((p.value() - 1e4).abs() < 1e-8);
    }

    #[test]
    fn bar_to_kilopascal() {
        let p: Pressure<f64, Kilopascal> = Pressure::<f64, Bar>::new(1.0).to();
        assert!((p.value() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn roundtrip_bar_atmosphere_bar() {
        let orig = Pressure::<f64, Bar>::new(2.5);
        let back: Pressure<f64, Bar> = orig.to::<Atmosphere>().to();
        assert!((back.value() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Pressure::<f64, Bar>::new(1.0) + Pressure::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut p = Pressure::<f64, Bar>::new(1.0);
        p += Pressure::new(0.5);
        assert_eq!(p.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Pressure::<f64, Bar>::new(3.0) - Pressure::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut p = Pressure::<f64, Bar>::new(3.0);
        p -= Pressure::new(1.0);
        assert_eq!(p.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Pressure::<f64, Bar>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Pressure::<f64, Bar>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut p = Pressure::<f64, Bar>::new(2.0);
        p *= 3.0;
        assert_eq!(p.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Pressure::<f64, Bar>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut p = Pressure::<f64, Bar>::new(6.0);
        p /= 2.0;
        assert_eq!(p.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Pressure::<f64, Bar>::new(6.0) / Pressure::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let p = Pressure::<f64, Bar>::new(1.0);
        assert_eq!(p, Pressure::new(1.0));
        assert_ne!(p, Pressure::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Pressure::<f64, Bar>::new(1.0);
        let b = Pressure::<f64, Bar>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Pressure::<f64, Bar>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Pressure::<f64, Bar>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let p = Pressure::<f64, Bar>::new(1.0);
        let nan = Pressure::<f64, Bar>::new(f64::NAN);
        assert_eq!(p.min(nan).value(), 1.0);
        assert_eq!(nan.min(p).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let p = Pressure::<f64, Bar>::new(1.0);
        let nan = Pressure::<f64, Bar>::new(f64::NAN);
        assert_eq!(p.max(nan).value(), 1.0);
        assert_eq!(nan.max(p).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Pressure::<f64, Bar>::new(1.0);
        let hi = Pressure::<f64, Bar>::new(2.0);
        assert_eq!(Pressure::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Pressure::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Pressure::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Pressure::<f64, Bar>::new(2.0);
        let hi = Pressure::<f64, Bar>::new(1.0);
        Pressure::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            Pressure::<f64, Bar>::new(1.0),
            Pressure::new(2.0),
            Pressure::new(3.0),
        ];
        let total: Pressure<f64, Bar> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Pressure::<f64, Bar>::new(1.0),
            Pressure::new(2.0),
            Pressure::new(3.0),
        ];
        let total: Pressure<f64, Bar> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Pressure<f64, Bar> = iter::empty::<Pressure<f64, Bar>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Pressure::<f64, Bar>::new(1.5).to_string(), "1.5 bar");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Pressure::<f64, Bar>::new(1.0)),
            "Pressure(1.0)"
        );
    }

    #[test]
    fn f32_atmosphere_to_bar() {
        let p: Pressure<f32, Bar> = Pressure::<f32, Atmosphere>::new(1.0_f32).to();
        assert!((p.value() - 1.013_25_f32).abs() < 1e-5_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Pressure::<f32, Bar>::new(1.0_f32) + Pressure::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

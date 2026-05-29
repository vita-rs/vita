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

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Time::<f64, Picosecond>::new(1.52).value(), 1.52);
    }

    #[test]
    fn from_scalar() {
        let t: Time<f64, Femtosecond> = Time::from(3.0);
        assert_eq!(t.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Time::<f64, Picosecond>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Time::<f64, Picosecond>::new(2.0);
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn femtosecond_to_picosecond() {
        let ps: Time<f64, Picosecond> = Time::<f64, Femtosecond>::new(1000.0).to();
        assert!((ps.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn picosecond_to_femtosecond() {
        let fs: Time<f64, Femtosecond> = Time::<f64, Picosecond>::new(1.0).to();
        assert!((fs.value() - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn picosecond_to_nanosecond() {
        let ns: Time<f64, Nanosecond> = Time::<f64, Picosecond>::new(1000.0).to();
        assert!((ns.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn picosecond_to_second() {
        let s: Time<f64, Second> = Time::<f64, Picosecond>::new(1e12).to();
        assert!((s.value() - 1.0).abs() < 1e-4);
    }

    #[test]
    fn atomic_time_to_picosecond() {
        let ps: Time<f64, Picosecond> = Time::<f64, AtomicTime>::new(1.0).to();
        assert!((ps.value() - 2.418_884_326_586_4e-5).abs() < 1e-18);
    }

    #[test]
    fn roundtrip_nanosecond_femtosecond_nanosecond() {
        let orig = Time::<f64, Nanosecond>::new(0.5);
        let back: Time<f64, Nanosecond> = orig.to::<Femtosecond>().to();
        assert!((back.value() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Time::<f64, Picosecond>::new(1.0) + Time::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut t = Time::<f64, Picosecond>::new(1.0);
        t += Time::new(0.5);
        assert_eq!(t.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Time::<f64, Picosecond>::new(3.0) - Time::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut t = Time::<f64, Picosecond>::new(3.0);
        t -= Time::new(1.0);
        assert_eq!(t.value(), 2.0);
    }

    #[test]
    fn rem() {
        let r = Time::<f64, Picosecond>::new(7.0) % Time::new(3.0);
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign() {
        let mut t = Time::<f64, Picosecond>::new(7.0);
        t %= Time::new(3.0);
        assert_eq!(t.value(), 1.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Time::<f64, Picosecond>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Time::<f64, Picosecond>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut t = Time::<f64, Picosecond>::new(2.0);
        t *= 3.0;
        assert_eq!(t.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Time::<f64, Picosecond>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut t = Time::<f64, Picosecond>::new(6.0);
        t /= 2.0;
        assert_eq!(t.value(), 3.0);
    }

    #[test]
    fn rem_scalar() {
        let r = Time::<f64, Picosecond>::new(7.0) % 3.0;
        assert_eq!(r.value(), 1.0);
    }

    #[test]
    fn rem_assign_scalar() {
        let mut t = Time::<f64, Picosecond>::new(7.0);
        t %= 3.0;
        assert_eq!(t.value(), 1.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Time::<f64, Picosecond>::new(6.0) / Time::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Time::<f64, Picosecond>::new(1.0);
        assert_eq!(a, Time::new(1.0));
        assert_ne!(a, Time::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Time::<f64, Picosecond>::new(1.0);
        let b = Time::<f64, Picosecond>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Time::<f64, Picosecond>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Time::<f64, Picosecond>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let t = Time::<f64, Picosecond>::new(1.0);
        let nan = Time::<f64, Picosecond>::new(f64::NAN);
        assert_eq!(t.min(nan).value(), 1.0);
        assert_eq!(nan.min(t).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let t = Time::<f64, Picosecond>::new(1.0);
        let nan = Time::<f64, Picosecond>::new(f64::NAN);
        assert_eq!(t.max(nan).value(), 1.0);
        assert_eq!(nan.max(t).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Time::<f64, Picosecond>::new(1.0);
        let hi = Time::<f64, Picosecond>::new(2.0);
        assert_eq!(Time::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Time::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Time::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Time::<f64, Picosecond>::new(2.0);
        let hi = Time::<f64, Picosecond>::new(1.0);
        Time::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn signum() {
        assert_eq!(Time::<f64, Picosecond>::new(3.0).signum(), 1.0);
        assert_eq!(Time::<f64, Picosecond>::new(-3.0).signum(), -1.0);
    }

    #[test]
    fn copysign() {
        let t = Time::<f64, Picosecond>::new(3.0);
        let sign = Time::<f64, Picosecond>::new(-1.0);
        assert_eq!(t.copysign(sign).value(), -3.0);
        assert_eq!((-t).copysign(t).value(), 3.0);
    }

    #[test]
    fn floor() {
        assert_eq!(Time::<f64, Picosecond>::new(2.7).floor().value(), 2.0);
        assert_eq!(Time::<f64, Picosecond>::new(-2.3).floor().value(), -3.0);
    }

    #[test]
    fn ceil() {
        assert_eq!(Time::<f64, Picosecond>::new(2.3).ceil().value(), 3.0);
        assert_eq!(Time::<f64, Picosecond>::new(-2.7).ceil().value(), -2.0);
    }

    #[test]
    fn round() {
        assert_eq!(Time::<f64, Picosecond>::new(2.5).round().value(), 3.0);
        assert_eq!(Time::<f64, Picosecond>::new(-2.5).round().value(), -3.0);
    }

    #[test]
    fn round_ties_even() {
        assert_eq!(
            Time::<f64, Picosecond>::new(2.5).round_ties_even().value(),
            2.0
        );
        assert_eq!(
            Time::<f64, Picosecond>::new(3.5).round_ties_even().value(),
            4.0
        );
    }

    #[test]
    fn trunc() {
        assert_eq!(Time::<f64, Picosecond>::new(2.7).trunc().value(), 2.0);
        assert_eq!(Time::<f64, Picosecond>::new(-2.7).trunc().value(), -2.0);
    }

    #[test]
    fn fract() {
        assert!((Time::<f64, Picosecond>::new(2.75).fract().value() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn div_euclid() {
        let q = Time::<f64, Picosecond>::new(7.0).div_euclid(Time::new(3.0));
        assert_eq!(q, 2.0);
    }

    #[test]
    fn rem_euclid() {
        let r = Time::<f64, Picosecond>::new(-7.0).rem_euclid(Time::new(3.0));
        assert_eq!(r.value(), 2.0);
    }

    #[test]
    fn mul_add() {
        let r = Time::<f64, Picosecond>::new(2.0).mul_add(3.0, Time::new(1.0));
        assert_eq!(r.value(), 7.0);
    }

    #[test]
    fn hypot() {
        let h = Time::<f64, Picosecond>::new(3.0).hypot(Time::new(4.0));
        assert!((h.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn is_nan() {
        assert!(Time::<f64, Picosecond>::new(f64::NAN).is_nan());
        assert!(!Time::<f64, Picosecond>::new(1.0).is_nan());
    }

    #[test]
    fn is_infinite() {
        assert!(Time::<f64, Picosecond>::new(f64::INFINITY).is_infinite());
        assert!(!Time::<f64, Picosecond>::new(1.0).is_infinite());
    }

    #[test]
    fn is_finite() {
        assert!(Time::<f64, Picosecond>::new(1.0).is_finite());
        assert!(!Time::<f64, Picosecond>::new(f64::INFINITY).is_finite());
        assert!(!Time::<f64, Picosecond>::new(f64::NAN).is_finite());
    }

    #[test]
    fn is_sign_positive() {
        assert!(Time::<f64, Picosecond>::new(1.0).is_sign_positive());
        assert!(!Time::<f64, Picosecond>::new(-1.0).is_sign_positive());
    }

    #[test]
    fn is_sign_negative() {
        assert!(Time::<f64, Picosecond>::new(-1.0).is_sign_negative());
        assert!(!Time::<f64, Picosecond>::new(1.0).is_sign_negative());
    }

    #[test]
    fn sum_owned() {
        let v = [
            Time::<f64, Picosecond>::new(1.0),
            Time::new(2.0),
            Time::new(3.0),
        ];
        let total: Time<f64, Picosecond> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Time::<f64, Picosecond>::new(1.0),
            Time::new(2.0),
            Time::new(3.0),
        ];
        let total: Time<f64, Picosecond> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Time<f64, Picosecond> = iter::empty::<Time<f64, Picosecond>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Time::<f64, Femtosecond>::new(1.5).to_string(), "1.5 fs");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Time::<f64, Picosecond>::new(1.0)),
            "Time(1.0)"
        );
    }

    #[test]
    fn f32_femtosecond_to_picosecond() {
        let ps: Time<f32, Picosecond> = Time::<f32, Femtosecond>::new(1000.0_f32).to();
        assert!((ps.value() - 1.0_f32).abs() < 1e-6_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Time::<f32, Picosecond>::new(1.0_f32) + Time::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

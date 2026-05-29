//! Energy quantities and unit markers.
//!
//! The canonical unit is the **hartree** (Eₕ).
//!
//! | Type | Symbol | Eₕ per unit |
//! |---|---|---|
//! | [`Hartree`] | Eₕ | 1 |
//! | [`KilocaloriePerMole`] | kcal mol⁻¹ | 1 / 627.509474063 |
//! | [`KilojoulePerMole`] | kJ mol⁻¹ | 1 / 2625.49963948 |
//! | [`ElectronVolt`] | eV | 1 / 27.211386245981 |
//! | [`Wavenumber`] | cm⁻¹ | 1 / 219474.63136314 |
//! | [`MilliElectronVolt`] | meV | 1 / 27211.386245981 |
//! | [`Joule`] | J | 1 / 4.3597447222060e-18 |

use crate::units::quantity::define_quantity;

/// Marker trait for energy units.
///
/// Implement this on a zero-sized type to define a new energy unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of hartrees
/// per one unit of `Self`.
pub trait EnergyUnit {
    /// Hartrees per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Eₕ"`, `"eV"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// An energy parameterised by scalar type `V` and unit marker `U`.
    Energy,
    EnergyUnit
);

/// The hartree (Eₕ) — canonical energy unit.
///
/// 1 Eₕ ≈ 4.3597447222060e-18 J (CODATA 2022).
pub struct Hartree;

impl EnergyUnit for Hartree {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Eₕ";
}

/// The kilocalorie per mole (kcal mol⁻¹).
///
/// 1 kcal mol⁻¹ ≈ 1 / 627.509474063 Eₕ (CODATA 2022, derived).
pub struct KilocaloriePerMole;

impl EnergyUnit for KilocaloriePerMole {
    const TO_CANONICAL: f64 = 1.0 / 627.509_474_063;
    const SYMBOL: &'static str = "kcal mol⁻¹";
}

/// The kilojoule per mole (kJ mol⁻¹).
///
/// 1 kJ mol⁻¹ ≈ 1 / 2625.49963948 Eₕ (CODATA 2022, derived).
pub struct KilojoulePerMole;

impl EnergyUnit for KilojoulePerMole {
    const TO_CANONICAL: f64 = 1.0 / 2625.499_639_48;
    const SYMBOL: &'static str = "kJ mol⁻¹";
}

/// The electron volt (eV) (CODATA 2022).
///
/// 1 eV ≈ 1 / 27.211386245981 Eₕ.
pub struct ElectronVolt;

impl EnergyUnit for ElectronVolt {
    const TO_CANONICAL: f64 = 1.0 / 27.211_386_245_981;
    const SYMBOL: &'static str = "eV";
}

/// The wavenumber (cm⁻¹) — inverse centimetres (CODATA 2022).
///
/// 1 cm⁻¹ ≈ 1 / 219474.63136314 Eₕ.
pub struct Wavenumber;

impl EnergyUnit for Wavenumber {
    const TO_CANONICAL: f64 = 1.0 / 219_474.631_363_14;
    const SYMBOL: &'static str = "cm⁻¹";
}

/// The milli-electron volt (meV) (CODATA 2022).
///
/// 1 meV ≈ 1 / 27211.386245981 Eₕ.
pub struct MilliElectronVolt;

impl EnergyUnit for MilliElectronVolt {
    const TO_CANONICAL: f64 = 1.0 / 27_211.386_245_981;
    const SYMBOL: &'static str = "meV";
}

/// The joule (J) — SI base unit of energy (CODATA 2022).
///
/// 1 J ≈ 1 / 4.3597447222060e-18 Eₕ.
pub struct Joule;

impl EnergyUnit for Joule {
    const TO_CANONICAL: f64 = 1.0 / 4.359_744_722_206_0e-18;
    const SYMBOL: &'static str = "J";
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Energy::<f64, Hartree>::new(2.5).value(), 2.5);
    }

    #[test]
    fn from_scalar() {
        let e: Energy<f64, ElectronVolt> = Energy::from(3.0);
        assert_eq!(e.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Energy::<f64, Hartree>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Energy::<f64, Hartree>::new(1.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn hartree_to_kilocalorie_per_mole() {
        let e: Energy<f64, KilocaloriePerMole> = Energy::<f64, Hartree>::new(1.0).to();
        assert!((e.value() - 627.509_474_063).abs() < 1e-9);
    }

    #[test]
    fn kilocalorie_per_mole_to_hartree() {
        let h: Energy<f64, Hartree> = Energy::<f64, KilocaloriePerMole>::new(627.509_474_063).to();
        assert!((h.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn hartree_to_electron_volt() {
        let ev: Energy<f64, ElectronVolt> = Energy::<f64, Hartree>::new(1.0).to();
        assert!((ev.value() - 27.211_386_245_981).abs() < 1e-12);
    }

    #[test]
    fn joule_to_hartree() {
        let h: Energy<f64, Hartree> = Energy::<f64, Joule>::new(4.359_744_722_206_0e-18).to();
        assert!((h.value() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn hartree_to_wavenumber() {
        let wn: Energy<f64, Wavenumber> = Energy::<f64, Hartree>::new(1.0).to();
        assert!((wn.value() - 219_474.631_363_14).abs() < 1e-6);
    }

    #[test]
    fn roundtrip_kilocalorie_per_mole_kilojoule_per_mole_kilocalorie_per_mole() {
        let orig = Energy::<f64, KilocaloriePerMole>::new(10.0);
        let back: Energy<f64, KilocaloriePerMole> = orig.to::<KilojoulePerMole>().to();
        assert!((back.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Energy::<f64, Hartree>::new(1.0) + Energy::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut e = Energy::<f64, Hartree>::new(1.0);
        e += Energy::new(0.5);
        assert_eq!(e.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Energy::<f64, Hartree>::new(3.0) - Energy::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut e = Energy::<f64, Hartree>::new(3.0);
        e -= Energy::new(1.0);
        assert_eq!(e.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Energy::<f64, Hartree>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Energy::<f64, Hartree>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut e = Energy::<f64, Hartree>::new(2.0);
        e *= 3.0;
        assert_eq!(e.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Energy::<f64, Hartree>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut e = Energy::<f64, Hartree>::new(6.0);
        e /= 2.0;
        assert_eq!(e.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Energy::<f64, Hartree>::new(6.0) / Energy::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Energy::<f64, Hartree>::new(1.0);
        assert_eq!(a, Energy::new(1.0));
        assert_ne!(a, Energy::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Energy::<f64, Hartree>::new(1.0);
        let b = Energy::<f64, Hartree>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Energy::<f64, Hartree>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Energy::<f64, Hartree>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let e = Energy::<f64, Hartree>::new(1.0);
        let nan = Energy::<f64, Hartree>::new(f64::NAN);
        assert_eq!(e.min(nan).value(), 1.0);
        assert_eq!(nan.min(e).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let e = Energy::<f64, Hartree>::new(1.0);
        let nan = Energy::<f64, Hartree>::new(f64::NAN);
        assert_eq!(e.max(nan).value(), 1.0);
        assert_eq!(nan.max(e).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Energy::<f64, Hartree>::new(1.0);
        let hi = Energy::<f64, Hartree>::new(2.0);
        assert_eq!(Energy::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Energy::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Energy::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Energy::<f64, Hartree>::new(2.0);
        let hi = Energy::<f64, Hartree>::new(1.0);
        Energy::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            Energy::<f64, Hartree>::new(1.0),
            Energy::new(2.0),
            Energy::new(3.0),
        ];
        let total: Energy<f64, Hartree> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Energy::<f64, Hartree>::new(1.0),
            Energy::new(2.0),
            Energy::new(3.0),
        ];
        let total: Energy<f64, Hartree> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Energy<f64, Hartree> = iter::empty::<Energy<f64, Hartree>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Energy::<f64, ElectronVolt>::new(1.5).to_string(), "1.5 eV");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Energy::<f64, Hartree>::new(1.0)),
            "Energy(1.0)"
        );
    }

    #[test]
    fn f32_hartree_to_electron_volt() {
        let ev: Energy<f32, ElectronVolt> = Energy::<f32, Hartree>::new(1.0_f32).to();
        assert!((ev.value() - 27.211_386_f32).abs() < 1e-4_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Energy::<f32, Hartree>::new(1.0_f32) + Energy::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

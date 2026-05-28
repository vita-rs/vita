//! Force quantities and unit markers.
//!
//! The canonical unit is the **hartree per bohr** (Eₕ a₀⁻¹).
//!
//! | Type | Symbol | Eₕ a₀⁻¹ per unit |
//! |---|---|---|
//! | [`HartreePerBohr`] | Eₕ a₀⁻¹ | 1 |
//! | [`KilocaloriePerMolePerAngstrom`] | kcal mol⁻¹ Å⁻¹ | 0.529177210544 / 627.509474063 |
//! | [`KilojoulePerMolePerNanometer`] | kJ mol⁻¹ nm⁻¹ | 0.0529177210544 / 2625.49963948 |
//! | [`ElectronVoltPerAngstrom`] | eV Å⁻¹ | 0.529177210544 / 27.211386245981 |
//! | [`Newton`] | N | 1 / 8.2387235038e-8 |
//! | [`Piconewton`] | pN | 1e-12 / 8.2387235038e-8 |

use crate::units::quantity::define_quantity;

/// Marker trait for force units.
///
/// Implement this on a zero-sized type to define a new force unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of hartrees per
/// bohr (Eₕ a₀⁻¹) per one unit of `Self`.
pub trait ForceUnit {
    /// Eₕ a₀⁻¹ per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Eₕ a₀⁻¹"`, `"N"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A force parameterised by scalar type `V` and unit marker `U`.
    Force,
    ForceUnit
);

/// The hartree per bohr (Eₕ a₀⁻¹) — canonical force unit (atomic unit of force, CODATA 2022).
///
/// 1 Eₕ a₀⁻¹ ≈ 8.2387235038e-8 N.
pub struct HartreePerBohr;

impl ForceUnit for HartreePerBohr {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Eₕ a₀⁻¹";
}

/// The kilocalorie per mole per ångström (kcal mol⁻¹ Å⁻¹).
///
/// 1 kcal mol⁻¹ Å⁻¹ ≈ 0.529177210544 / 627.509474063 Eₕ a₀⁻¹ (CODATA 2022, derived).
pub struct KilocaloriePerMolePerAngstrom;

impl ForceUnit for KilocaloriePerMolePerAngstrom {
    const TO_CANONICAL: f64 = 0.529_177_210_544 / 627.509_474_063;
    const SYMBOL: &'static str = "kcal mol⁻¹ Å⁻¹";
}

/// The kilojoule per mole per nanometre (kJ mol⁻¹ nm⁻¹).
///
/// 1 kJ mol⁻¹ nm⁻¹ ≈ 0.0529177210544 / 2625.49963948 Eₕ a₀⁻¹ (CODATA 2022, derived).
pub struct KilojoulePerMolePerNanometer;

impl ForceUnit for KilojoulePerMolePerNanometer {
    const TO_CANONICAL: f64 = 0.052_917_721_054_4 / 2625.499_639_48;
    const SYMBOL: &'static str = "kJ mol⁻¹ nm⁻¹";
}

/// The electronvolt per ångström (eV Å⁻¹) — standard DFT force unit.
///
/// 1 eV Å⁻¹ ≈ 0.529177210544 / 27.211386245981 Eₕ a₀⁻¹ (CODATA 2022, derived).
pub struct ElectronVoltPerAngstrom;

impl ForceUnit for ElectronVoltPerAngstrom {
    const TO_CANONICAL: f64 = 0.529_177_210_544 / 27.211_386_245_981;
    const SYMBOL: &'static str = "eV Å⁻¹";
}

/// The newton (N) — SI unit of force (CODATA 2022).
///
/// 1 N ≈ 1 / 8.2387235038e-8 Eₕ a₀⁻¹.
pub struct Newton;

impl ForceUnit for Newton {
    const TO_CANONICAL: f64 = 1.0 / 8.238_723_503_8e-8;
    const SYMBOL: &'static str = "N";
}

/// The piconewton (pN) — used in single-molecule force spectroscopy (AFM, optical tweezers).
///
/// 1 pN ≈ 1e-12 / 8.2387235038e-8 Eₕ a₀⁻¹ (CODATA 2022).
pub struct Piconewton;

impl ForceUnit for Piconewton {
    const TO_CANONICAL: f64 = 1e-12 / 8.238_723_503_8e-8;
    const SYMBOL: &'static str = "pN";
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::format;
    use alloc::string::ToString;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(Force::<f64, HartreePerBohr>::new(1.52).value(), 1.52);
    }

    #[test]
    fn from_scalar() {
        let f: Force<f64, Newton> = Force::from(3.0);
        assert_eq!(f.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Force::<f64, HartreePerBohr>::default().value(), 0.0_f64);
    }

    #[test]
    fn copy_and_clone() {
        let a = Force::<f64, HartreePerBohr>::new(2.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn hartree_per_bohr_to_kilocalorie_per_mole_per_angstrom() {
        let f: Force<f64, KilocaloriePerMolePerAngstrom> =
            Force::<f64, HartreePerBohr>::new(1.0).to();
        assert!((f.value() - 1185.821_047).abs() < 1e-6);
    }

    #[test]
    fn kilocalorie_per_mole_per_angstrom_to_hartree_per_bohr() {
        let f: Force<f64, HartreePerBohr> =
            Force::<f64, KilocaloriePerMolePerAngstrom>::new(1185.821_047_391_503).to();
        assert!((f.value() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn hartree_per_bohr_to_kilojoule_per_mole_per_nanometer() {
        let f: Force<f64, KilojoulePerMolePerNanometer> =
            Force::<f64, HartreePerBohr>::new(1.0).to();
        assert!((f.value() - 49_614.752_622_87).abs() < 1e-3);
    }

    #[test]
    fn hartree_per_bohr_to_electron_volt_per_angstrom() {
        let f: Force<f64, ElectronVoltPerAngstrom> = Force::<f64, HartreePerBohr>::new(1.0).to();
        assert!((f.value() - 51.422_067_511).abs() < 1e-6);
    }

    #[test]
    fn hartree_per_bohr_to_newton() {
        let f: Force<f64, Newton> = Force::<f64, HartreePerBohr>::new(1.0).to();
        assert!((f.value() - 8.238_723_503_8e-8).abs() < 1e-18);
    }

    #[test]
    fn hartree_per_bohr_to_piconewton() {
        let f: Force<f64, Piconewton> = Force::<f64, HartreePerBohr>::new(1.0).to();
        assert!((f.value() - 82_387.235_038).abs() < 1e-3);
    }

    #[test]
    fn roundtrip_kilocalorie_per_mole_per_angstrom_newton_kilocalorie_per_mole_per_angstrom() {
        let orig = Force::<f64, KilocaloriePerMolePerAngstrom>::new(10.0);
        let back: Force<f64, KilocaloriePerMolePerAngstrom> = orig.to::<Newton>().to();
        assert!((back.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = Force::<f64, HartreePerBohr>::new(1.0) + Force::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut f = Force::<f64, HartreePerBohr>::new(1.0);
        f += Force::new(0.5);
        assert_eq!(f.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = Force::<f64, HartreePerBohr>::new(3.0) - Force::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut f = Force::<f64, HartreePerBohr>::new(3.0);
        f -= Force::new(1.0);
        assert_eq!(f.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!((-Force::<f64, HartreePerBohr>::new(1.5)).value(), -1.5);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!((Force::<f64, HartreePerBohr>::new(2.0) * 3.0).value(), 6.0);
    }

    #[test]
    fn mul_assign_scalar() {
        let mut f = Force::<f64, HartreePerBohr>::new(2.0);
        f *= 3.0;
        assert_eq!(f.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!((Force::<f64, HartreePerBohr>::new(6.0) / 2.0).value(), 3.0);
    }

    #[test]
    fn div_assign_scalar() {
        let mut f = Force::<f64, HartreePerBohr>::new(6.0);
        f /= 2.0;
        assert_eq!(f.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = Force::<f64, HartreePerBohr>::new(6.0) / Force::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = Force::<f64, HartreePerBohr>::new(1.0);
        assert_eq!(a, Force::new(1.0));
        assert_ne!(a, Force::new(2.0));
    }

    #[test]
    fn ord() {
        let a = Force::<f64, HartreePerBohr>::new(1.0);
        let b = Force::<f64, HartreePerBohr>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(Force::<f64, HartreePerBohr>::new(-3.0).abs().value(), 3.0);
        assert_eq!(Force::<f64, HartreePerBohr>::new(3.0).abs().value(), 3.0);
    }

    #[test]
    fn min_ignores_nan() {
        let f = Force::<f64, HartreePerBohr>::new(1.0);
        let nan = Force::<f64, HartreePerBohr>::new(f64::NAN);
        assert_eq!(f.min(nan).value(), 1.0);
        assert_eq!(nan.min(f).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let f = Force::<f64, HartreePerBohr>::new(1.0);
        let nan = Force::<f64, HartreePerBohr>::new(f64::NAN);
        assert_eq!(f.max(nan).value(), 1.0);
        assert_eq!(nan.max(f).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = Force::<f64, HartreePerBohr>::new(1.0);
        let hi = Force::<f64, HartreePerBohr>::new(2.0);
        assert_eq!(Force::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(Force::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(Force::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = Force::<f64, HartreePerBohr>::new(2.0);
        let hi = Force::<f64, HartreePerBohr>::new(1.0);
        Force::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            Force::<f64, HartreePerBohr>::new(1.0),
            Force::new(2.0),
            Force::new(3.0),
        ];
        let total: Force<f64, HartreePerBohr> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            Force::<f64, HartreePerBohr>::new(1.0),
            Force::new(2.0),
            Force::new(3.0),
        ];
        let total: Force<f64, HartreePerBohr> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: Force<f64, HartreePerBohr> = iter::empty::<Force<f64, HartreePerBohr>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(Force::<f64, Piconewton>::new(1.5).to_string(), "1.5 pN");
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Force::<f64, HartreePerBohr>::new(1.0)),
            "Force(1.0)"
        );
    }

    #[test]
    fn f32_hartree_per_bohr_to_kilocalorie_per_mole_per_angstrom() {
        let f: Force<f32, KilocaloriePerMolePerAngstrom> =
            Force::<f32, HartreePerBohr>::new(1.0_f32).to();
        assert!((f.value() - 1185.821_f32).abs() < 0.01_f32);
    }

    #[test]
    fn f32_add() {
        let sum = Force::<f32, HartreePerBohr>::new(1.0_f32) + Force::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

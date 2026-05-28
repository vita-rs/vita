//! Force-constant quantities and unit markers.
//!
//! The canonical unit is the **hartree per bohr squared** (Eₕ a₀⁻²).
//!
//! | Type | Symbol | Eₕ a₀⁻² per unit |
//! |---|---|---|
//! | [`HartreePerBohrSquared`] | Eₕ a₀⁻² | 1 |
//! | [`KcalPerMolPerAngstromSquared`] | kcal mol⁻¹ Å⁻² | 0.529177210544² / 627.509474063 |
//! | [`KjPerMolPerNanometerSquared`] | kJ mol⁻¹ nm⁻² | 0.0529177210544² / 2625.49963948 |
//! | [`NewtonPerMeter`] | N m⁻¹ | 5.29177210544e-11 / 8.2387235038e-8 |

use crate::units::quantity::define_quantity;

/// Marker trait for force-constant units.
///
/// Implement this on a zero-sized type to define a new force-constant unit.
/// [`TO_CANONICAL`][Self::TO_CANONICAL] must give the number of hartrees per
/// bohr squared (Eₕ a₀⁻²) per one unit of `Self`.
pub trait ForceConstantUnit {
    /// Eₕ a₀⁻² per one unit of `Self`.
    const TO_CANONICAL: f64;
    /// Display symbol (e.g. `"Eₕ a₀⁻²"`, `"N m⁻¹"`).
    const SYMBOL: &'static str;
}

define_quantity!(
    /// A force constant parameterised by scalar type `V` and unit marker `U`.
    ForceConstant,
    ForceConstantUnit
);

/// The hartree per bohr squared (Eₕ a₀⁻²) — canonical force-constant unit
/// (atomic unit of force constant, CODATA 2022).
///
/// 1 Eₕ a₀⁻² ≈ 1556.893105 N m⁻¹.
pub struct HartreePerBohrSquared;

impl ForceConstantUnit for HartreePerBohrSquared {
    const TO_CANONICAL: f64 = 1.0;
    const SYMBOL: &'static str = "Eₕ a₀⁻²";
}

/// The kilocalorie per mole per ångström squared (kcal mol⁻¹ Å⁻²).
///
/// 1 kcal mol⁻¹ Å⁻² ≈ 0.529177210544² / 627.509474063 Eₕ a₀⁻² (CODATA 2022, derived).
pub struct KcalPerMolPerAngstromSquared;

impl ForceConstantUnit for KcalPerMolPerAngstromSquared {
    const TO_CANONICAL: f64 = 0.529_177_210_544 * 0.529_177_210_544 / 627.509_474_063;
    const SYMBOL: &'static str = "kcal mol⁻¹ Å⁻²";
}

/// The kilojoule per mole per nanometre squared (kJ mol⁻¹ nm⁻²).
///
/// 1 kJ mol⁻¹ nm⁻² ≈ 0.0529177210544² / 2625.49963948 Eₕ a₀⁻² (CODATA 2022, derived).
pub struct KjPerMolPerNanometerSquared;

impl ForceConstantUnit for KjPerMolPerNanometerSquared {
    const TO_CANONICAL: f64 = 0.052_917_721_054_4 * 0.052_917_721_054_4 / 2625.499_639_48;
    const SYMBOL: &'static str = "kJ mol⁻¹ nm⁻²";
}

/// The newton per metre (N m⁻¹) — SI unit of force constant (CODATA 2022).
///
/// 1 N m⁻¹ ≈ 5.29177210544e-11 / 8.2387235038e-8 Eₕ a₀⁻².
pub struct NewtonPerMeter;

impl ForceConstantUnit for NewtonPerMeter {
    const TO_CANONICAL: f64 = 5.291_772_105_44e-11 / 8.238_723_503_8e-8;
    const SYMBOL: &'static str = "N m⁻¹";
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::format;
    use alloc::string::ToString;
    use core::iter;

    #[test]
    fn new_value_roundtrip() {
        assert_eq!(
            ForceConstant::<f64, HartreePerBohrSquared>::new(1.52).value(),
            1.52
        );
    }

    #[test]
    fn from_scalar() {
        let k: ForceConstant<f64, NewtonPerMeter> = ForceConstant::from(3.0);
        assert_eq!(k.value(), 3.0);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(
            ForceConstant::<f64, HartreePerBohrSquared>::default().value(),
            0.0_f64
        );
    }

    #[test]
    fn copy_and_clone() {
        let a = ForceConstant::<f64, HartreePerBohrSquared>::new(2.0);
        let b = a;
        let c = a.clone();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn hartree_per_bohr_squared_to_kcal_per_mol_per_angstrom_squared() {
        let k: ForceConstant<f64, KcalPerMolPerAngstromSquared> =
            ForceConstant::<f64, HartreePerBohrSquared>::new(1.0).to();
        assert!((k.value() - 2240.877_014).abs() < 1e-3);
    }

    #[test]
    fn kcal_per_mol_per_angstrom_squared_to_hartree_per_bohr_squared() {
        let k: ForceConstant<f64, HartreePerBohrSquared> =
            ForceConstant::<f64, KcalPerMolPerAngstromSquared>::new(2240.877_013_907_053).to();
        assert!((k.value() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn hartree_per_bohr_squared_to_kj_per_mol_per_nanometer_squared() {
        let k: ForceConstant<f64, KjPerMolPerNanometerSquared> =
            ForceConstant::<f64, HartreePerBohrSquared>::new(1.0).to();
        assert!((k.value() - 937_582.942_62).abs() < 1e-2);
    }

    #[test]
    fn hartree_per_bohr_squared_to_newton_per_meter() {
        let k: ForceConstant<f64, NewtonPerMeter> =
            ForceConstant::<f64, HartreePerBohrSquared>::new(1.0).to();
        assert!((k.value() - 1556.893_105).abs() < 1e-6);
    }

    #[test]
    fn roundtrip_kcal_per_mol_per_angstrom_squared_newton_per_meter_kcal_per_mol_per_angstrom_squared()
     {
        let orig = ForceConstant::<f64, KcalPerMolPerAngstromSquared>::new(10.0);
        let back: ForceConstant<f64, KcalPerMolPerAngstromSquared> =
            orig.to::<NewtonPerMeter>().to();
        assert!((back.value() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn add() {
        let sum = ForceConstant::<f64, HartreePerBohrSquared>::new(1.0) + ForceConstant::new(2.5);
        assert_eq!(sum.value(), 3.5);
    }

    #[test]
    fn add_assign() {
        let mut k = ForceConstant::<f64, HartreePerBohrSquared>::new(1.0);
        k += ForceConstant::new(0.5);
        assert_eq!(k.value(), 1.5);
    }

    #[test]
    fn sub() {
        let diff = ForceConstant::<f64, HartreePerBohrSquared>::new(3.0) - ForceConstant::new(1.0);
        assert_eq!(diff.value(), 2.0);
    }

    #[test]
    fn sub_assign() {
        let mut k = ForceConstant::<f64, HartreePerBohrSquared>::new(3.0);
        k -= ForceConstant::new(1.0);
        assert_eq!(k.value(), 2.0);
    }

    #[test]
    fn neg() {
        assert_eq!(
            (-ForceConstant::<f64, HartreePerBohrSquared>::new(1.5)).value(),
            -1.5
        );
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(
            (ForceConstant::<f64, HartreePerBohrSquared>::new(2.0) * 3.0).value(),
            6.0
        );
    }

    #[test]
    fn mul_assign_scalar() {
        let mut k = ForceConstant::<f64, HartreePerBohrSquared>::new(2.0);
        k *= 3.0;
        assert_eq!(k.value(), 6.0);
    }

    #[test]
    fn div_scalar() {
        assert_eq!(
            (ForceConstant::<f64, HartreePerBohrSquared>::new(6.0) / 2.0).value(),
            3.0
        );
    }

    #[test]
    fn div_assign_scalar() {
        let mut k = ForceConstant::<f64, HartreePerBohrSquared>::new(6.0);
        k /= 2.0;
        assert_eq!(k.value(), 3.0);
    }

    #[test]
    fn div_same_unit_yields_ratio() {
        let ratio = ForceConstant::<f64, HartreePerBohrSquared>::new(6.0) / ForceConstant::new(2.0);
        assert_eq!(ratio, 3.0);
    }

    #[test]
    fn eq() {
        let a = ForceConstant::<f64, HartreePerBohrSquared>::new(1.0);
        assert_eq!(a, ForceConstant::new(1.0));
        assert_ne!(a, ForceConstant::new(2.0));
    }

    #[test]
    fn ord() {
        let a = ForceConstant::<f64, HartreePerBohrSquared>::new(1.0);
        let b = ForceConstant::<f64, HartreePerBohrSquared>::new(2.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn abs() {
        assert_eq!(
            ForceConstant::<f64, HartreePerBohrSquared>::new(-3.0)
                .abs()
                .value(),
            3.0
        );
        assert_eq!(
            ForceConstant::<f64, HartreePerBohrSquared>::new(3.0)
                .abs()
                .value(),
            3.0
        );
    }

    #[test]
    fn min_ignores_nan() {
        let k = ForceConstant::<f64, HartreePerBohrSquared>::new(1.0);
        let nan = ForceConstant::<f64, HartreePerBohrSquared>::new(f64::NAN);
        assert_eq!(k.min(nan).value(), 1.0);
        assert_eq!(nan.min(k).value(), 1.0);
    }

    #[test]
    fn max_ignores_nan() {
        let k = ForceConstant::<f64, HartreePerBohrSquared>::new(1.0);
        let nan = ForceConstant::<f64, HartreePerBohrSquared>::new(f64::NAN);
        assert_eq!(k.max(nan).value(), 1.0);
        assert_eq!(nan.max(k).value(), 1.0);
    }

    #[test]
    fn clamp() {
        let lo = ForceConstant::<f64, HartreePerBohrSquared>::new(1.0);
        let hi = ForceConstant::<f64, HartreePerBohrSquared>::new(2.0);
        assert_eq!(ForceConstant::new(1.5_f64).clamp(lo, hi).value(), 1.5);
        assert_eq!(ForceConstant::new(0.5_f64).clamp(lo, hi).value(), 1.0);
        assert_eq!(ForceConstant::new(3.0_f64).clamp(lo, hi).value(), 2.0);
    }

    #[test]
    #[should_panic]
    fn clamp_panics_when_lo_gt_hi() {
        let lo = ForceConstant::<f64, HartreePerBohrSquared>::new(2.0);
        let hi = ForceConstant::<f64, HartreePerBohrSquared>::new(1.0);
        ForceConstant::new(1.5_f64).clamp(lo, hi);
    }

    #[test]
    fn sum_owned() {
        let v = [
            ForceConstant::<f64, HartreePerBohrSquared>::new(1.0),
            ForceConstant::new(2.0),
            ForceConstant::new(3.0),
        ];
        let total: ForceConstant<f64, HartreePerBohrSquared> = v.iter().copied().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_borrowed() {
        let v = [
            ForceConstant::<f64, HartreePerBohrSquared>::new(1.0),
            ForceConstant::new(2.0),
            ForceConstant::new(3.0),
        ];
        let total: ForceConstant<f64, HartreePerBohrSquared> = v.iter().sum();
        assert_eq!(total.value(), 6.0);
    }

    #[test]
    fn sum_empty() {
        let total: ForceConstant<f64, HartreePerBohrSquared> =
            iter::empty::<ForceConstant<f64, HartreePerBohrSquared>>().sum();
        assert_eq!(total.value(), 0.0);
    }

    #[test]
    fn display() {
        assert_eq!(
            ForceConstant::<f64, NewtonPerMeter>::new(1.5).to_string(),
            "1.5 N m⁻¹"
        );
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!(
                "{:?}",
                ForceConstant::<f64, HartreePerBohrSquared>::new(1.0)
            ),
            "ForceConstant(1.0)"
        );
    }

    #[test]
    fn f32_hartree_per_bohr_squared_to_kcal_per_mol_per_angstrom_squared() {
        let k: ForceConstant<f32, KcalPerMolPerAngstromSquared> =
            ForceConstant::<f32, HartreePerBohrSquared>::new(1.0_f32).to();
        assert!((k.value() - 2240.877_f32).abs() < 0.1_f32);
    }

    #[test]
    fn f32_add() {
        let sum =
            ForceConstant::<f32, HartreePerBohrSquared>::new(1.0_f32) + ForceConstant::new(2.0_f32);
        assert_eq!(sum.value(), 3.0_f32);
    }
}

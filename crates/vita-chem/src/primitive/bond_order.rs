/// The order of a chemical bond between two sites.
///
/// `BondOrder` encodes bond multiplicity as a per-bond enum variant, not as
/// multiple parallel edges. `C=C` is one [`Double`](BondOrder::Double) bond,
/// not two [`Single`](BondOrder::Single) bonds.
///
/// ## Aromatic vs. Kekulé representation
///
/// The same aromatic ring may be stored in two equivalent forms:
///
/// | Form | benzene C–C |
/// |------|-------------|
/// | Aromatic | six [`Aromatic`](BondOrder::Aromatic) bonds |
/// | Kekulé | alternating [`Single`](BondOrder::Single)/[`Double`](BondOrder::Double) |
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondOrder {
    /// A single covalent bond (σ only).
    Single,
    /// A double bond (σ + π).
    Double,
    /// A triple bond (σ + 2π).
    Triple,
    /// A quadruple bond (σ + 2π + δ).
    Quadruple,
    /// A quintuple bond (σ + 2π + 2δ).
    Quintuple,
    /// A hextuple bond (σ + 2π + 2δ + 2φ).
    Hextuple,
    /// A bond in a π-aromatic ring system (Hückel 4*n*+2).
    Aromatic,
}

impl BondOrder {
    /// Returns the bond order as a floating-point number.
    #[inline]
    pub fn as_f64(self) -> f64 {
        match self {
            Self::Single => 1.0,
            Self::Aromatic => 1.5,
            Self::Double => 2.0,
            Self::Triple => 3.0,
            Self::Quadruple => 4.0,
            Self::Quintuple => 5.0,
            Self::Hextuple => 6.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn copy_and_clone() {
        let a = BondOrder::Triple;
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn eq() {
        assert_eq!(BondOrder::Single, BondOrder::Single);
        assert_ne!(BondOrder::Single, BondOrder::Double);
        assert_ne!(BondOrder::Double, BondOrder::Triple);
        assert_ne!(BondOrder::Triple, BondOrder::Aromatic);
    }

    #[test]
    fn debug() {
        assert_eq!(format!("{:?}", BondOrder::Single), "Single");
        assert_eq!(format!("{:?}", BondOrder::Double), "Double");
        assert_eq!(format!("{:?}", BondOrder::Triple), "Triple");
        assert_eq!(format!("{:?}", BondOrder::Quadruple), "Quadruple");
        assert_eq!(format!("{:?}", BondOrder::Quintuple), "Quintuple");
        assert_eq!(format!("{:?}", BondOrder::Hextuple), "Hextuple");
        assert_eq!(format!("{:?}", BondOrder::Aromatic), "Aromatic");
    }

    #[test]
    fn all_variants_distinct() {
        use std::collections::HashSet;
        let variants = [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Quadruple,
            BondOrder::Quintuple,
            BondOrder::Hextuple,
            BondOrder::Aromatic,
        ];
        let set: HashSet<_> = variants.into_iter().collect();
        assert_eq!(set.len(), 7);
    }

    #[test]
    fn as_f64() {
        assert_eq!(BondOrder::Single.as_f64(), 1.0);
        assert_eq!(BondOrder::Double.as_f64(), 2.0);
        assert_eq!(BondOrder::Triple.as_f64(), 3.0);
        assert_eq!(BondOrder::Quadruple.as_f64(), 4.0);
        assert_eq!(BondOrder::Quintuple.as_f64(), 5.0);
        assert_eq!(BondOrder::Hextuple.as_f64(), 6.0);
        assert_eq!(BondOrder::Aromatic.as_f64(), 1.5);
    }
}

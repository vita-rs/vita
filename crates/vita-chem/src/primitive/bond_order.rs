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
    pub fn to_f64(self) -> f64 {
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
    fn integer_orders_map_to_their_multiplicity() {
        assert_eq!(BondOrder::Single.to_f64(), 1.0);
        assert_eq!(BondOrder::Double.to_f64(), 2.0);
        assert_eq!(BondOrder::Triple.to_f64(), 3.0);
        assert_eq!(BondOrder::Quadruple.to_f64(), 4.0);
        assert_eq!(BondOrder::Quintuple.to_f64(), 5.0);
        assert_eq!(BondOrder::Hextuple.to_f64(), 6.0);
    }

    #[test]
    fn aromatic_maps_to_one_and_a_half() {
        assert_eq!(BondOrder::Aromatic.to_f64(), 1.5);
    }
}

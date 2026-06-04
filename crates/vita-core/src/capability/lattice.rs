use crate::{Lattice, Scalar};

/// The periodic [`Lattice`] of a system.
///
/// [`lattice`](HasLattice::lattice) returns the basis vectors generating the system's
/// translational symmetry. Implementing this trait asserts periodicity.
pub trait HasLattice<V: Scalar> {
    /// Returns the system's lattice.
    fn lattice(&self) -> Lattice<V>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::length::{Angstrom, Length};
    use crate::units::volume::CubicAngstrom;

    struct Crystal {
        lattice: Lattice<f64>,
    }
    impl HasLattice<f64> for Crystal {
        fn lattice(&self) -> Lattice<f64> {
            self.lattice
        }
    }

    #[test]
    fn lattice() {
        let crystal = Crystal {
            lattice: Lattice::cubic(Length::<f64, Angstrom>::new(2.0)).unwrap(),
        };
        assert_eq!(crystal.lattice().volume::<CubicAngstrom>().value(), 8.0);
    }
}

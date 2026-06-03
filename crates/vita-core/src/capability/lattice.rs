use crate::{Lattice, Scalar};

/// The periodic [`Lattice`] of a system.
///
/// [`lattice`](HasLattice::lattice) returns the basis vectors generating the system's
/// translational symmetry. Implementing this trait asserts periodicity.
pub trait HasLattice<V: Scalar> {
    /// Returns the system's lattice.
    fn lattice(&self) -> Lattice<V>;
}

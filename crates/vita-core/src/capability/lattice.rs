use crate::{Lattice, Scalar};

/// The periodic [`Lattice`] of a system.
///
/// [`lattice`](HasLattice::lattice) returns the group of translations the system
/// repeats under, and the quotient of space it induces. Implementing this trait
/// asserts periodicity.
pub trait HasLattice<V: Scalar> {
    /// Returns the system's lattice.
    fn lattice(&self) -> Lattice<V>;
}

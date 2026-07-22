use vita_core::HasIsotopes;

use super::{Composition, Constituent};
use crate::HasFormalCharges;

/// The isotopic composition of a molecule: every site counted under its
/// declared nuclide, with the net formal charge.
///
/// The finer of the two folds — [`HasIsotopes`] declares one nuclide per
/// site, and every count keeps it. For natural-mixture counting, use
/// [`elemental`](super::elemental).
///
/// # Complexity
///
/// O(V · log V) time and O(V) space, over the molecule's `V` sites; the log
/// factor orders the counts canonically.
pub fn isotopic<M: HasIsotopes + HasFormalCharges>(mol: &M) -> Composition {
    Composition::from_counts(
        mol.isotopes()
            .map(|isotope| (Constituent::Nuclide(isotope), 1)),
        mol.formal_charges().map(i32::from).sum(),
    )
}

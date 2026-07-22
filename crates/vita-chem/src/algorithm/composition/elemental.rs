use vita_core::HasElements;

use super::{Composition, Constituent};
use crate::HasFormalCharges;

/// The elemental composition of a molecule: every site counted under its
/// element at natural isotopic precision, with the net formal charge.
///
/// The coarser of the two folds — isotopic declarations, if any, are not
/// consulted. For a composition that keeps them, use
/// [`isotopic`](super::isotopic).
///
/// # Complexity
///
/// O(V · log V) time and O(V) space, over the molecule's `V` sites; the log
/// factor orders the counts canonically.
pub fn elemental<M: HasElements + HasFormalCharges>(mol: &M) -> Composition {
    Composition::from_counts(
        mol.elements()
            .map(|element| (Constituent::Element(element), 1)),
        mol.formal_charges().map(i32::from).sum(),
    )
}

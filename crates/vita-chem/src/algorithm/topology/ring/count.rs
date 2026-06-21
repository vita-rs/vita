use vita_core::HasSites;

use crate::HasBonds;
use crate::topology::connectivity::components;

/// Number of independent rings in a molecule.
///
/// Equals the cycle rank of the molecular graph, μ = E − V + C, where C is the
/// number of connected components. This is the size of any minimum cycle basis
/// — the value [`Rings::len`](super::Rings::len) reports — but is obtained
/// without enumerating the rings. Acyclic and empty molecules have zero.
///
/// # Complexity
///
/// O(V + E) time.
pub fn count<M: HasBonds + HasSites>(mol: &M) -> usize {
    let v = mol.sites().count();
    let e = mol.bonds().count();
    let c = components(mol).len();
    e + c - v
}

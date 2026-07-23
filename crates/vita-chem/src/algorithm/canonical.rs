//! Canonical labeling: the molecule's portable identity.
//!
//! [`canonicalize`] assigns every site its rank in a total order fixed by
//! the molecular graph and the caller's coloring — never by the order the
//! sites arrived — and records the relabeled graph, with its symmetry
//! classes, as a [`Canonical`] form. Molecules the coloring makes
//! isomorphic share the form exactly, so identity is `==`.

mod canonicalize;

pub use canonicalize::{Canonical, canonicalize};

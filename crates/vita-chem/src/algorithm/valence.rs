//! Per-site electron bookkeeping over the Lewis graph.
//!
//! Four accountings of one bond set: [`valence`] sums the integer bond orders;
//! [`lone_pairs`] pairs up what the covalent limit leaves, every bond split
//! evenly; [`oxidation_state`] takes the ionic limit, every bond awarded wholly
//! to its more electronegative end; [`steric_numbers`] counts the electron
//! domains those electrons fall into rather than the electrons. None guesses
//! where its count is not exact — above all across a delocalized
//! [`Aromatic`](crate::BondOrder::Aromatic) bond, which has no integer order
//! until the ring is kekulized — and each leaves such a site unanswered.

mod explicit;
mod lone_pairs;
mod oxidation_state;
mod steric_numbers;

pub use explicit::valence;
pub use lone_pairs::lone_pairs;
pub use oxidation_state::oxidation_state;
pub use steric_numbers::{StericNumbers, steric_numbers};

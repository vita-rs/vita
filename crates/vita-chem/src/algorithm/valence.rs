//! Per-site electron bookkeeping over the Lewis graph.
//!
//! Three accountings of one bond set: [`valence`] sums the integer bond
//! orders; [`lone_pairs`] pairs up what the covalent limit leaves, every bond
//! split evenly; [`oxidation_state`] takes the ionic limit, every bond
//! awarded wholly to its more electronegative end. Each is a point function
//! answering for a single site, and each returns `None` rather than guess
//! where its count is not exact — above all across a delocalized
//! [`Aromatic`](crate::BondOrder::Aromatic) bond, which has no integer order
//! until the ring is kekulized.

mod explicit;
mod lone_pairs;
mod oxidation_state;

pub use explicit::valence;
pub use lone_pairs::lone_pairs;
pub use oxidation_state::oxidation_state;

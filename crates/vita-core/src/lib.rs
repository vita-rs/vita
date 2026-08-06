//! Foundational vocabulary for the vita ecosystem: the questions a physical system can
//! answer, not the structures that store the answers.
//!
//! > What we observe is not nature itself, but nature exposed to our method of
//! > questioning.
//! >
//! > — Werner Heisenberg
//!
//! A system is never a concrete type here. Each kind of question it can answer is a
//! trait: [`HasElements`] what occupies a site, [`HasPositions`] where it sits,
//! [`HasNetCharge`] the total charge. Code bounds on exactly the capabilities it needs
//! and stays blind to the storage behind them.
//!
//! # Sites
//!
//! Per-site data is keyed on an opaque [`SiteId`]; [`HasSites`] enumerates those keys and
//! is the supertrait of every per-site capability. A capability is a keyed getter paired
//! with a value iterator in [`sites`](HasSites::sites) order — zip the two for keyed
//! access. System-wide quantities ([`HasLattice`], [`HasNetCharge`]) are standalone,
//! single-valued traits.
//!
//! # Quantities
//!
//! Physical values carry their dimension and unit in the type via [`units`]; spatial ones
//! use the three-dimensional primitives in [`tensor`]. [`Quantity`] names what the two share —
//! the operations that leave a dimension unchanged — so a tensor's element may be a bare
//! [`Scalar`] (`f32` or `f64`) or any dimensioned quantity.
//!
//! # Geometry
//!
//! [`geometry`] reads what a system's placement determines:
//! [`measure`](geometry::measure) answers for a fixed tuple of sites — a separation, an
//! angle, a dihedral, a handedness — [`moment`](geometry::moment) for the sites as a
//! whole, the center they are spread about and the spread itself, and
//! [`proximity`](geometry::proximity) for what lies near what.

mod capability;
mod primitive;

pub mod geometry;
pub mod tensor;
pub mod units;

pub mod prelude;

pub use primitive::{Element, Isotope, Lattice, Quantity, ReciprocalLattice, Scalar, SiteId};

pub use capability::{
    HasAccelerations, HasElements, HasIsotopes, HasLattice, HasMasses, HasNetCharge, HasPositions,
    HasSites, HasVelocities,
};

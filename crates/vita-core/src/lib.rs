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
//! with a `(SiteId, _)` iterator — storage order is never implied. System-wide quantities
//! ([`HasLattice`], [`HasNetCharge`]) are standalone, single-valued traits.
//!
//! # Quantities
//!
//! Physical values carry their dimension and unit in the type via [`units`]; geometry uses
//! the three-dimensional primitives in [`tensor`]. Both are generic over the [`Scalar`]
//! element type, `f32` or `f64`.

mod capability;
mod primitive;

pub mod tensor;
pub mod units;

pub mod prelude;

pub use primitive::{Element, Isotope, Lattice, Scalar, SiteId};

pub use capability::{
    HasAccelerations, HasElements, HasIsotopes, HasLattice, HasMasses, HasNetCharge, HasPositions,
    HasSites, HasVelocities,
};

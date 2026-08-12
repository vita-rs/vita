//! Discrete chemistry for the vita ecosystem: what a declared molecular
//! structure determines, never what a heuristic would guess.
//!
//! > Let us learn to dream, gentlemen, then perhaps we shall find the truth —
//! > but let us beware of publishing our dreams before they have been put to
//! > the proof by the waking understanding.
//! >
//! > — August Kekulé
//!
//! A molecule here is what it is in [`vita_core`]: the questions it can
//! answer. This crate asks the chemical ones — [`HasBondOrders`],
//! [`HasFormalCharges`], [`HasStereoConfigurations`], spoken in the bond
//! vocabulary of [`BondId`], [`BondOrder`], and [`StereoConfiguration`] —
//! and derives what the answers entail. Every algorithm bounds on exactly
//! the capabilities it consumes; every result is an owned value or an
//! overlay view binding the derived fact back onto the unchanged molecule
//! as a capability of its own. Perception widens what a molecule can
//! answer, never what it is.
//!
//! # Determination
//!
//! Each answer is fixed by the discrete structure alone: absolutely, as in
//! [`topology`] and [`canonical`], or within a model declared outright —
//! Hückel's rule for [`aromaticity`], p-orbital availability for
//! [`conjugation`]. Where determination runs out, the answer is honestly
//! absent: a delocalized bond has no integer order until kekulized, a
//! d-block site no fixed electron count. Where only convention could
//! choose, nothing is offered: a configuration is a coset, not a named
//! handedness, and no tautomer is chosen as dominant. Identity is computed;
//! values are never invented.
//!
//! # Modules
//!
//! [`topology`] reads the bare skeleton; [`valence`] keeps the electron
//! books; [`aromaticity`] and [`conjugation`] perceive within their models;
//! [`canonical`] and [`stereo`] settle identity; [`isomorphism`] finds
//! structure inside structure; [`composition`] takes the particle
//! inventory; [`fingerprint`] measures similarity; [`notation`] renders
//! facts as standard symbols.

mod algorithm;
mod capability;
mod primitive;

pub mod prelude;

pub use primitive::{
    BondId, BondOrder, CoordinationGeometry, StereoConfiguration, StereoDescriptor, StereoKind,
    StereoLocus,
};

pub use capability::{
    HasAromaticity, HasBondOrders, HasBonds, HasFormalCharges, HasPartialCharges,
    HasRadicalElectrons, HasStereoConfigurations,
};

pub use algorithm::{
    aromaticity, canonical, composition, conjugation, fingerprint, isomorphism, notation, stereo,
    topology, valence,
};

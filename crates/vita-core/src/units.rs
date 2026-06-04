//! Dimensionally-typed physical quantities.
//!
//! Each physical dimension is a submodule exposing a quantity newtype `Q<V, U>` — a
//! scalar `V` tagged with a zero-sized unit marker `U` — alongside the marker trait and
//! unit types implementing it. Converting between units is the explicit `.to()`.

mod dimensions;
mod quantity;

pub use dimensions::{
    acceleration, amount_of_substance, angle, area, charge, concentration, density, dipole_moment,
    energy, force, force_constant, frequency, length, mass, momentum, pressure, temperature, time,
    velocity, volume,
};

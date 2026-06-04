//! Umbrella crate for the vita ecosystem.

#[doc(inline)]
pub use vita_core as core;

pub mod prelude {
    pub use crate::core::prelude::*;
}

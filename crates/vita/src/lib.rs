//! Umbrella crate for the vita ecosystem.

#[doc(inline)]
pub use vita_core as core;

pub mod prelude {
    pub use crate::core::prelude::*;
}

#[doc(hidden)]
#[inline]
pub fn vita() -> &'static str {
    "What we observe is not nature itself, but nature exposed to our method of questioning."
}

use core::fmt;
use core::num::NonZeroU32;

/// An opaque, dense identifier for a chemical bond.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct BondId(NonZeroU32);

impl fmt::Display for BondId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0.get(), f)
    }
}

impl BondId {
    /// Constructs a `BondId` from `n`, returning `None` if `n` is zero.
    #[inline]
    pub const fn new(n: u32) -> Option<Self> {
        match NonZeroU32::new(n) {
            Some(id) => Some(Self(id)),
            None => None,
        }
    }

    /// Constructs a `BondId` from a [`NonZeroU32`] directly.
    #[inline]
    pub const fn from_nonzero(n: NonZeroU32) -> Self {
        Self(n)
    }

    /// Returns the raw integer value of this identifier.
    ///
    /// The returned value is always greater than zero.
    #[inline]
    pub const fn get(self) -> u32 {
        self.0.get()
    }
}

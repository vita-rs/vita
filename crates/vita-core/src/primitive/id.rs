use core::fmt;
use core::num::NonZeroU32;

/// An opaque, dense identifier for a simulation site (e.g., an atom).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct SiteId(NonZeroU32);

impl fmt::Display for SiteId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0.get(), f)
    }
}

impl SiteId {
    /// Constructs a `SiteId` from `n`, returning `None` if `n` is zero.
    #[inline]
    pub const fn new(n: u32) -> Option<Self> {
        match NonZeroU32::new(n) {
            Some(id) => Some(Self(id)),
            None => None,
        }
    }

    /// Constructs a `SiteId` from a [`NonZeroU32`] directly.
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

#[cfg(test)]
mod tests {
    use super::*;

    fn id(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    #[test]
    fn new_with_zero_is_none() {
        assert!(SiteId::new(0).is_none());
    }

    #[test]
    fn new_with_one_is_some() {
        assert!(SiteId::new(1).is_some());
    }

    #[test]
    fn get_returns_the_value_passed_to_new() {
        assert_eq!(id(7).get(), 7);
    }

    #[test]
    fn from_nonzero_preserves_the_value() {
        let n = NonZeroU32::new(7).unwrap();
        assert_eq!(SiteId::from_nonzero(n).get(), 7);
    }

    #[test]
    fn display_shows_the_integer_value() {
        assert_eq!(format!("{}", id(42)), "42");
    }

    #[test]
    fn new_accepts_the_maximum_value() {
        assert_eq!(id(u32::MAX).get(), u32::MAX);
    }

    #[test]
    fn ids_order_by_their_integer_value() {
        assert!(id(1) < id(2));
        assert!(id(2) < id(10));
        assert!(id(2) > id(1));
    }

    #[test]
    fn ids_compare_equal_when_their_values_match() {
        assert_eq!(id(5), id(5));
        assert_ne!(id(5), id(6));
    }

    #[test]
    fn option_is_the_same_size_as_the_raw_integer() {
        assert_eq!(
            core::mem::size_of::<Option<SiteId>>(),
            core::mem::size_of::<u32>(),
        );
    }
}

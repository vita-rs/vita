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

    #[test]
    fn new() {
        assert!(SiteId::new(0).is_none());
        assert!(SiteId::new(1).is_some());
    }

    #[test]
    fn from_nonzero() {
        let n = NonZeroU32::new(7).unwrap();
        assert_eq!(SiteId::from_nonzero(n).get(), 7);
    }

    #[test]
    fn new_get_roundtrip() {
        assert_eq!(SiteId::new(42).unwrap().get(), 42);
    }

    #[test]
    fn copy_and_clone() {
        let a = SiteId::new(1).unwrap();
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn eq() {
        let a = SiteId::new(3).unwrap();
        assert_eq!(a, SiteId::new(3).unwrap());
        assert_ne!(a, SiteId::new(4).unwrap());
    }

    #[test]
    fn ord() {
        let a = SiteId::new(1).unwrap();
        let b = SiteId::new(2).unwrap();
        assert!(a < b);
        assert!(b > a);
        assert_eq!(a.cmp(&a), ::core::cmp::Ordering::Equal);
    }

    #[test]
    fn debug() {
        assert_eq!(format!("{:?}", SiteId::new(1).unwrap()), "SiteId(1)");
    }

    #[test]
    fn display() {
        assert_eq!(format!("{}", SiteId::new(42).unwrap()), "42");
    }

    #[test]
    fn option_is_same_size_as_u32() {
        assert_eq!(
            ::core::mem::size_of::<Option<SiteId>>(),
            ::core::mem::size_of::<u32>()
        );
    }
}

use core::fmt;

use crate::Element;

/// A nuclide: an [`Element`] together with its mass number *A*.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Isotope {
    element: Element,
    mass_number: u16,
}

impl Isotope {
    /// Constructs a nuclide of `element` with mass number `mass_number`,
    /// returning `None` if `mass_number` is below the atomic number (a
    /// negative neutron count).
    #[inline]
    pub const fn new(element: Element, mass_number: u16) -> Option<Self> {
        if mass_number < element.atomic_number() as u16 {
            return None;
        }
        Some(Self {
            element,
            mass_number,
        })
    }

    /// Constructs a nuclide of `element` with `neutron_count` neutrons,
    /// returning `None` if the mass number would overflow [`u16`].
    #[inline]
    pub const fn from_neutron_count(element: Element, neutron_count: u16) -> Option<Self> {
        match (element.atomic_number() as u16).checked_add(neutron_count) {
            Some(mass_number) => Some(Self {
                element,
                mass_number,
            }),
            None => None,
        }
    }

    /// Returns the element (proton count) of this nuclide.
    #[inline]
    pub const fn element(self) -> Element {
        self.element
    }

    /// Returns the atomic number *Z* (the proton count).
    #[inline]
    pub const fn atomic_number(self) -> u8 {
        self.element.atomic_number()
    }

    /// Returns the mass number *A* (protons plus neutrons).
    #[inline]
    pub const fn mass_number(self) -> u16 {
        self.mass_number
    }

    /// Returns the neutron count *N = A − Z*.
    #[inline]
    pub const fn neutron_count(self) -> u16 {
        self.mass_number - self.element.atomic_number() as u16
    }
}

impl fmt::Display for Isotope {
    /// Formats as `symbol-A`, e.g. `C-12`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{}", self.element.symbol(), self.mass_number)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn carbon() -> Element {
        Element::new(6).unwrap()
    }

    #[test]
    fn new() {
        let c12 = Isotope::new(carbon(), 12).unwrap();
        assert_eq!(c12.atomic_number(), 6);
        assert_eq!(c12.mass_number(), 12);
    }

    #[test]
    fn new_rejects_mass_below_atomic_number() {
        assert!(Isotope::new(carbon(), 5).is_none());
        assert!(Isotope::new(carbon(), 6).is_some());
    }

    #[test]
    fn from_neutron_count() {
        let c13 = Isotope::from_neutron_count(carbon(), 7).unwrap();
        assert_eq!(c13.mass_number(), 13);
        assert_eq!(c13.neutron_count(), 7);
    }

    #[test]
    fn from_neutron_count_rejects_overflow() {
        assert!(Isotope::from_neutron_count(carbon(), u16::MAX).is_none());
    }

    #[test]
    fn element() {
        assert_eq!(Isotope::new(carbon(), 12).unwrap().element(), carbon());
    }

    #[test]
    fn atomic_number() {
        assert_eq!(Isotope::new(carbon(), 12).unwrap().atomic_number(), 6);
    }

    #[test]
    fn mass_number() {
        assert_eq!(Isotope::new(carbon(), 12).unwrap().mass_number(), 12);
    }

    #[test]
    fn neutron_count() {
        assert_eq!(Isotope::new(carbon(), 14).unwrap().neutron_count(), 8);
    }

    #[test]
    fn neutron_count_is_zero_for_hydrogen_one() {
        let h1 = Isotope::new(Element::new(1).unwrap(), 1).unwrap();
        assert_eq!(h1.neutron_count(), 0);
    }

    #[test]
    fn copy_and_clone() {
        let a = Isotope::new(carbon(), 12).unwrap();
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn eq() {
        let a = Isotope::new(carbon(), 12).unwrap();
        assert_eq!(a, Isotope::new(carbon(), 12).unwrap());
        assert_ne!(a, Isotope::new(carbon(), 13).unwrap());
    }

    #[test]
    fn ord() {
        let c12 = Isotope::new(carbon(), 12).unwrap();
        let c13 = Isotope::new(carbon(), 13).unwrap();
        let n12 = Isotope::new(Element::new(7).unwrap(), 12).unwrap();
        assert!(c12 < c13);
        assert!(c13 < n12);
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Isotope::new(carbon(), 12).unwrap()),
            "Isotope { element: Element(6), mass_number: 12 }"
        );
    }

    #[test]
    fn display() {
        assert_eq!(format!("{}", Isotope::new(carbon(), 12).unwrap()), "C-12");
        assert_eq!(
            format!("{}", Isotope::new(Element::new(92).unwrap(), 238).unwrap()),
            "U-238"
        );
    }

    #[test]
    fn option_is_same_size_as_isotope() {
        assert_eq!(
            ::core::mem::size_of::<Option<Isotope>>(),
            ::core::mem::size_of::<Isotope>()
        );
    }
}

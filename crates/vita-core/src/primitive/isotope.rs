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

    fn hydrogen() -> Element {
        Element::new(1).unwrap()
    }

    #[test]
    fn element_returns_the_underlying_element() {
        assert_eq!(Isotope::new(carbon(), 14).unwrap().element(), carbon());
    }

    #[test]
    fn atomic_number_returns_the_proton_count() {
        assert_eq!(Isotope::new(carbon(), 14).unwrap().atomic_number(), 6);
    }

    #[test]
    fn mass_number_returns_the_nucleon_count() {
        assert_eq!(Isotope::new(carbon(), 14).unwrap().mass_number(), 14);
    }

    #[test]
    fn neutron_count_is_mass_number_minus_atomic_number() {
        assert_eq!(Isotope::new(carbon(), 14).unwrap().neutron_count(), 8);
    }

    #[test]
    fn from_neutron_count_builds_the_mass_number() {
        assert_eq!(
            Isotope::from_neutron_count(carbon(), 8)
                .unwrap()
                .mass_number(),
            14,
        );
    }

    #[test]
    fn displays_as_symbol_and_mass_number() {
        assert_eq!(format!("{}", Isotope::new(carbon(), 12).unwrap()), "C-12");
    }

    #[test]
    fn new_rejects_a_mass_number_below_the_atomic_number() {
        assert_eq!(Isotope::new(carbon(), 5), None);
    }

    #[test]
    fn from_neutron_count_rejects_an_overflowing_mass_number() {
        assert_eq!(Isotope::from_neutron_count(carbon(), u16::MAX), None);
    }

    #[test]
    fn mass_number_equal_to_the_atomic_number_has_zero_neutrons() {
        assert_eq!(Isotope::new(carbon(), 6).unwrap().neutron_count(), 0);
    }

    #[test]
    fn the_largest_non_overflowing_neutron_count_is_valid() {
        let iso = Isotope::from_neutron_count(carbon(), u16::MAX - 6).unwrap();
        assert_eq!(iso.mass_number(), u16::MAX);
    }

    #[test]
    fn isotopes_order_by_element_then_mass_number() {
        assert!(Isotope::new(hydrogen(), 2).unwrap() < Isotope::new(carbon(), 12).unwrap());
        assert!(Isotope::new(carbon(), 12).unwrap() < Isotope::new(carbon(), 14).unwrap());
    }

    #[test]
    fn isotopes_are_equal_exactly_when_element_and_mass_number_match() {
        let carbon_12 = Isotope::new(carbon(), 12).unwrap();
        assert_eq!(carbon_12, Isotope::new(carbon(), 12).unwrap());
        assert_ne!(carbon_12, Isotope::new(carbon(), 14).unwrap());
        assert_ne!(carbon_12, Isotope::new(hydrogen(), 12).unwrap());
    }
}

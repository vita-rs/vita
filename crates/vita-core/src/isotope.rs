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

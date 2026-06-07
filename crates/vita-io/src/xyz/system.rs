use vita_core::tensor::Point3;
use vita_core::units::length::{Angstrom, Length, LengthUnit};
use vita_core::{Element, HasElements, HasPositions, HasSites, Scalar, SiteId};

/// A set of atoms with their Cartesian positions, as recorded by the XYZ format.
///
/// Provides [`HasSites`], [`HasElements`], and [`HasPositions`] — exactly the
/// capabilities the format records.
#[derive(Clone, Debug, PartialEq)]
pub struct System<V: Scalar = f64> {
    comment: Box<str>,
    elements: Box<[Element]>,
    positions: Box<[Point3<Length<V, Angstrom>>]>,
}

impl<V: Scalar> System<V> {
    /// Assembles a system from its parsed columns.
    ///
    /// `elements` and `positions` are parallel columns,
    /// they must have equal length.
    pub(super) fn from_parts(
        comment: Box<str>,
        elements: Box<[Element]>,
        positions: Box<[Point3<Length<V, Angstrom>>]>,
    ) -> Self {
        debug_assert_eq!(elements.len(), positions.len());
        debug_assert!(elements.len() <= u32::MAX as usize);
        Self {
            comment,
            elements,
            positions,
        }
    }

    /// Returns the comment line, kept verbatim from the file.
    #[inline]
    pub fn comment(&self) -> &str {
        &self.comment
    }
}

impl<V: Scalar> HasSites for System<V> {
    #[inline]
    fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        (1..=self.elements.len() as u32).map(|n| SiteId::new(n).unwrap())
    }

    #[inline]
    fn site_count(&self) -> usize {
        self.elements.len()
    }

    #[inline]
    fn contains_site(&self, site: SiteId) -> bool {
        (site.get() as usize) <= self.elements.len()
    }
}

impl<V: Scalar> HasElements for System<V> {
    #[inline]
    fn element(&self, site: SiteId) -> Element {
        self.elements[(site.get() - 1) as usize]
    }

    #[inline]
    fn elements(&self) -> impl Iterator<Item = (SiteId, Element)> + '_ {
        self.sites().zip(self.elements.iter().copied())
    }
}

impl<V: Scalar> HasPositions<V> for System<V> {
    #[inline]
    fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<V, U>> {
        self.positions[(site.get() - 1) as usize].map(|length| length.to())
    }

    #[inline]
    fn positions<U: LengthUnit>(
        &self,
    ) -> impl Iterator<Item = (SiteId, Point3<Length<V, U>>)> + '_ {
        self.sites().zip(
            self.positions
                .iter()
                .copied()
                .map(|position| position.map(|length| length.to::<U>())),
        )
    }
}

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
    fn elements(&self) -> impl Iterator<Item = Element> + '_ {
        self.elements.iter().copied()
    }
}

impl<V: Scalar> HasPositions<V> for System<V> {
    #[inline]
    fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<V, U>> {
        self.positions[(site.get() - 1) as usize].map(|length| length.to())
    }

    #[inline]
    fn positions<U: LengthUnit>(&self) -> impl Iterator<Item = Point3<Length<V, U>>> + '_ {
        self.positions
            .iter()
            .copied()
            .map(|position| position.map(|length| length.to::<U>()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::Quantity;
    use vita_core::units::length::Nanometer;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    fn system(comment: &str, atoms: &[(&str, f64, f64, f64)]) -> System {
        System::from_parts(
            comment.into(),
            atoms.iter().map(|&(symbol, ..)| elem(symbol)).collect(),
            atoms
                .iter()
                .map(|&(_, x, y, z)| Point3::new(Length::new(x), Length::new(y), Length::new(z)))
                .collect(),
        )
    }

    fn empty() -> System {
        system("", &[])
    }

    fn molecule() -> System {
        system(
            "test frame",
            &[
                ("C", 1.0, 2.0, 3.0),
                ("O", 4.0, 5.0, 6.0),
                ("H", 7.0, 8.0, 9.0),
            ],
        )
    }

    #[test]
    fn empty_system_has_no_sites() {
        let system = empty();
        assert_eq!(system.site_count(), 0);
        assert_eq!(system.sites().count(), 0);
    }

    #[test]
    fn empty_system_contains_no_site() {
        assert!(!empty().contains_site(s(1)));
    }

    #[test]
    fn empty_system_has_no_elements() {
        assert_eq!(empty().elements().count(), 0);
    }

    #[test]
    fn empty_system_has_no_positions() {
        assert_eq!(empty().positions::<Angstrom>().count(), 0);
    }

    #[test]
    fn site_count_matches_the_number_of_atoms() {
        assert_eq!(molecule().site_count(), 3);
    }

    #[test]
    fn sites_are_numbered_consecutively_from_one() {
        assert_eq!(
            molecule().sites().collect::<Vec<_>>(),
            vec![s(1), s(2), s(3)],
        );
    }

    #[test]
    fn contains_site_holds_for_present_sites() {
        let system = molecule();
        assert!(system.contains_site(s(1)));
        assert!(system.contains_site(s(3)));
    }

    #[test]
    fn element_returns_the_atom_at_each_site() {
        let system = molecule();
        assert_eq!(system.element(s(1)), elem("C"));
        assert_eq!(system.element(s(3)), elem("H"));
    }

    #[test]
    fn elements_are_listed_in_site_order() {
        assert_eq!(
            molecule().elements().collect::<Vec<_>>(),
            vec![elem("C"), elem("O"), elem("H")],
        );
    }

    #[test]
    fn position_returns_the_recorded_coordinates() {
        let system = molecule();
        assert_eq!(
            system.position::<Angstrom>(s(1)),
            Point3::new(Length::new(1.0), Length::new(2.0), Length::new(3.0)),
        );
        assert_eq!(
            system.position::<Angstrom>(s(3)),
            Point3::new(Length::new(7.0), Length::new(8.0), Length::new(9.0)),
        );
    }

    #[test]
    fn positions_are_listed_in_site_order() {
        assert_eq!(
            molecule().positions::<Angstrom>().collect::<Vec<_>>(),
            vec![
                Point3::new(Length::new(1.0), Length::new(2.0), Length::new(3.0)),
                Point3::new(Length::new(4.0), Length::new(5.0), Length::new(6.0)),
                Point3::new(Length::new(7.0), Length::new(8.0), Length::new(9.0)),
            ],
        );
    }

    #[test]
    fn comment_is_kept_verbatim() {
        assert_eq!(molecule().comment(), "test frame");
    }

    #[test]
    fn contains_site_fails_past_the_last_site() {
        assert!(!molecule().contains_site(s(4)));
    }

    #[test]
    fn position_converts_to_the_requested_unit() {
        let p = molecule().position::<Nanometer>(s(1));
        assert!((p.x.value() - 0.1).abs() < 1e-12);
        assert!((p.y.value() - 0.2).abs() < 1e-12);
        assert!((p.z.value() - 0.3).abs() < 1e-12);
    }

    #[test]
    fn supports_the_f32_scalar_type() {
        let system: System<f32> = System::from_parts(
            "".into(),
            Box::from([elem("H")]),
            Box::from([Point3::new(
                Length::new(1.5),
                Length::new(2.5),
                Length::new(3.5),
            )]),
        );
        assert_eq!(
            system.position::<Angstrom>(s(1)),
            Point3::new(Length::new(1.5_f32), Length::new(2.5), Length::new(3.5)),
        );
    }

    #[test]
    fn equality_accounts_for_the_comment() {
        let atoms: &[(&str, f64, f64, f64)] = &[("H", 0.0, 0.0, 0.0)];
        assert_eq!(system("a", atoms), system("a", atoms));
        assert_ne!(system("a", atoms), system("b", atoms));
    }
}

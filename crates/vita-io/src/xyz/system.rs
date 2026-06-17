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
    use vita_core::units::length::Nanometer;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn hydrogen() -> Element {
        Element::new(1).unwrap()
    }

    fn oxygen() -> Element {
        Element::new(8).unwrap()
    }

    fn angstrom(x: f64, y: f64, z: f64) -> Point3<Length<f64, Angstrom>> {
        Point3::new(Length::new(x), Length::new(y), Length::new(z))
    }

    fn water() -> System<f64> {
        System::from_parts(
            "water".into(),
            vec![oxygen(), hydrogen(), hydrogen()].into_boxed_slice(),
            vec![
                angstrom(0.0, 0.0, 0.0),
                angstrom(0.757, 0.586, 0.0),
                angstrom(-0.757, 0.586, 0.0),
            ]
            .into_boxed_slice(),
        )
    }

    fn empty() -> System<f64> {
        System::from_parts(
            "".into(),
            vec![].into_boxed_slice(),
            vec![].into_boxed_slice(),
        )
    }

    #[test]
    fn comment() {
        assert_eq!(water().comment(), "water");
    }

    #[test]
    fn comment_empty() {
        assert_eq!(empty().comment(), "");
    }

    #[test]
    fn site_count() {
        assert_eq!(water().site_count(), 3);
    }

    #[test]
    fn site_count_empty_is_zero() {
        assert_eq!(empty().site_count(), 0);
    }

    #[test]
    fn sites_are_one_based() {
        assert_eq!(water().sites().next(), Some(site(1)));
    }

    #[test]
    fn sites_are_sequential() {
        assert_eq!(
            water().sites().collect::<Vec<_>>(),
            vec![site(1), site(2), site(3)],
        );
    }

    #[test]
    fn sites_empty() {
        assert_eq!(empty().sites().count(), 0);
    }

    #[test]
    fn contains_site() {
        let sys = water();
        assert!(sys.contains_site(site(1)));
        assert!(sys.contains_site(site(3)));
        assert!(!sys.contains_site(site(4)));
    }

    #[test]
    fn element() {
        let sys = water();
        assert_eq!(sys.element(site(1)), oxygen());
        assert_eq!(sys.element(site(2)), hydrogen());
        assert_eq!(sys.element(site(3)), hydrogen());
    }

    #[test]
    fn elements() {
        let sys = water();
        assert_eq!(
            sys.elements().collect::<Vec<_>>(),
            vec![oxygen(), hydrogen(), hydrogen()],
        );
    }

    #[test]
    fn elements_empty() {
        assert_eq!(empty().elements().count(), 0);
    }

    #[test]
    fn position() {
        let sys = water();
        assert_eq!(sys.position::<Angstrom>(site(1)), angstrom(0.0, 0.0, 0.0));
        assert_eq!(
            sys.position::<Angstrom>(site(2)),
            angstrom(0.757, 0.586, 0.0),
        );
    }

    #[test]
    fn position_unit_conversion() {
        let sys = System::from_parts(
            "".into(),
            vec![hydrogen()].into_boxed_slice(),
            vec![angstrom(1.0, 0.0, 0.0)].into_boxed_slice(),
        );
        let p = sys.position::<Nanometer>(site(1));
        assert_eq!(p.x, Length::new(0.1));
        assert_eq!(p.y, Length::new(0.0));
        assert_eq!(p.z, Length::new(0.0));
    }

    #[test]
    fn positions() {
        let sys = water();
        assert_eq!(
            sys.positions::<Angstrom>().collect::<Vec<_>>(),
            vec![
                angstrom(0.0, 0.0, 0.0),
                angstrom(0.757, 0.586, 0.0),
                angstrom(-0.757, 0.586, 0.0),
            ],
        );
    }

    #[test]
    fn positions_empty() {
        assert_eq!(empty().positions::<Angstrom>().count(), 0);
    }

    #[test]
    fn f32_scalar() {
        let sys: System<f32> = System::from_parts(
            "c".into(),
            vec![hydrogen()].into_boxed_slice(),
            vec![Point3::new(
                Length::<f32, Angstrom>::new(1.5),
                Length::<f32, Angstrom>::new(0.0),
                Length::<f32, Angstrom>::new(0.0),
            )]
            .into_boxed_slice(),
        );
        assert_eq!(sys.site_count(), 1);
        assert_eq!(sys.element(site(1)), hydrogen());
        assert_eq!(
            sys.position::<Angstrom>(site(1)).x,
            Length::<f32, Angstrom>::new(1.5)
        );
    }

    #[test]
    fn clone_and_eq() {
        let a = water();
        let b = a.clone();
        assert_eq!(a, b);
        let c = System::from_parts(
            "other".into(),
            vec![hydrogen()].into_boxed_slice(),
            vec![angstrom(0.0, 0.0, 0.0)].into_boxed_slice(),
        );
        assert_ne!(a, c);
    }
}

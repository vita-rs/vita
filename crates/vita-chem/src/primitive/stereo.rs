use vita_core::SiteId;

use crate::BondId;

/// Where a stereogenic unit is anchored in a molecule.
///
/// Stereochemistry is carried by the arrangement of a site's substituents (a
/// coordination centre), of the substituents across a bond (a double bond's two
/// ends), or across an axis (an allene's two termini); `StereoLocus` names that
/// anchor. It labels a unit for reporting and reconciliation and carries no meaning
/// of its own — the arrangement that fixes which orderings are equivalent is the
/// unit's [`StereoKind`].
///
/// A site locus orders before a bond, a bond before an axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoLocus {
    /// A site whose substituents' arrangement is stereogenic.
    Site(SiteId),
    /// A bond whose ends' substituents' arrangement is stereogenic.
    Bond(BondId),
    /// An axis, named by its central site, whose termini's substituents'
    /// arrangement is stereogenic.
    Axis(SiteId),
}

impl StereoLocus {
    /// Whether this locus is the anchor `kind` lives on: a site for a coordination
    /// centre, a bond for a double bond, an axis for an allene.
    #[inline]
    pub const fn anchors(self, kind: StereoKind) -> bool {
        matches!(
            (self, kind),
            (
                Self::Site(_),
                StereoKind::Tetrahedral
                    | StereoKind::SquarePlanar
                    | StereoKind::TrigonalBipyramidal
                    | StereoKind::SquarePyramidal
                    | StereoKind::Octahedral
                    | StereoKind::TrigonalPrismatic,
            ) | (Self::Bond(_), StereoKind::CisTrans)
                | (Self::Axis(_), StereoKind::Allene)
        )
    }
}

/// The kind of a stereogenic unit: the idealised local geometry whose rotation
/// group fixes which of its neighbour orderings are equivalent.
///
/// A kind is a pure data key. It selects the permutation group under which a
/// configuration's neighbour ordering reduces, and the geometric reference against
/// which coordinates are perceived.
///
/// The kinds order by neighbour count, then configuration count, then locus —
/// a centre before a bond before an axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoKind {
    /// A tetrahedral centre (4 neighbours, 2 configurations).
    Tetrahedral,
    /// A double bond (4 neighbours, 2 configurations).
    CisTrans,
    /// An allene axis (4 neighbours, 2 configurations).
    Allene,
    /// A square-planar centre (4 neighbours, 3 configurations).
    SquarePlanar,
    /// A trigonal-bipyramidal centre (5 neighbours, 20 configurations).
    TrigonalBipyramidal,
    /// A square-pyramidal centre (5 neighbours, 30 configurations).
    SquarePyramidal,
    /// An octahedral centre (6 neighbours, 30 configurations).
    Octahedral,
    /// A trigonal-prismatic centre (6 neighbours, 120 configurations).
    TrigonalPrismatic,
}

impl StereoKind {
    /// Returns the number of neighbour slots the geometry arranges.
    #[inline]
    pub const fn slot_count(self) -> usize {
        match self {
            Self::Tetrahedral | Self::CisTrans | Self::Allene | Self::SquarePlanar => 4,
            Self::TrigonalBipyramidal | Self::SquarePyramidal => 5,
            Self::Octahedral | Self::TrigonalPrismatic => 6,
        }
    }

    /// Returns the number of distinct configurations the geometry admits — the
    /// stereoisomers its slots realise when every substituent differs.
    #[inline]
    pub const fn configuration_count(self) -> usize {
        match self {
            Self::Tetrahedral | Self::CisTrans | Self::Allene => 2,
            Self::SquarePlanar => 3,
            Self::TrigonalBipyramidal => 20,
            Self::SquarePyramidal | Self::Octahedral => 30,
            Self::TrigonalPrismatic => 120,
        }
    }

    /// Returns whether the geometry is chiral — whether a configuration and its
    /// mirror image are distinct, no rotation carrying one onto the other.
    #[inline]
    pub const fn is_chiral(self) -> bool {
        match self {
            Self::CisTrans | Self::SquarePlanar => false,
            Self::Tetrahedral
            | Self::Allene
            | Self::TrigonalBipyramidal
            | Self::SquarePyramidal
            | Self::Octahedral
            | Self::TrigonalPrismatic => true,
        }
    }
}

/// A stereo configuration at one [`StereoLocus`].
///
/// The *order* of the neighbours fixes the configuration: it is a reference
/// arrangement, and every reordering the unit's [`StereoKind`] treats as equivalent
/// denotes the same one. The library reads this order and never invents it — a
/// source (a wedge, coordinates, a SMILES `@`) already committed to it, and must
/// commit to the convention here.
///
/// Each kind fills its slots from the neighbour list; a chiral kind's reference order
/// is additionally the one whose positions span a positive signed volume, its mirror
/// the order the kind's reflection gives. An achiral kind — a double bond, a square
/// plane — equates an arrangement with its mirror, so its order fixes no handedness.
///
/// - [`Tetrahedral`](StereoKind::Tetrahedral) — the four substituents
///   `[n0, n1, n2, n3]`; reference `(n1 − n0) · ((n2 − n0) × (n3 − n0))`.
/// - [`CisTrans`](StereoKind::CisTrans) — one end's two, then the other's,
///   `[e1a, e1b, e2a, e2b]`, `e1a` on `e2a`'s side.
/// - [`Allene`](StereoKind::Allene) — one terminus's two, then the other's,
///   `[t1a, t1b, t2a, t2b]`; reference the twist across the axis, its mirror
///   swapping `t1a` and `t1b`.
/// - [`SquarePlanar`](StereoKind::SquarePlanar) — the four in cyclic order
///   `[n0, n1, n2, n3]`, `n0`/`n2` and `n1`/`n3` trans.
/// - [`TrigonalBipyramidal`](StereoKind::TrigonalBipyramidal) — the two axial, then
///   the three equatorial, `[a0, a1, e0, e1, e2]`; reference
///   `(a0 − a1) · ((e1 − e0) × (e2 − e0))`.
/// - [`SquarePyramidal`](StereoKind::SquarePyramidal) — the apical, then the four
///   basal in cyclic order, `[p, b0, b1, b2, b3]`; reference
///   `(p − b0) · ((b1 − b0) × (b2 − b0))`.
/// - [`Octahedral`](StereoKind::Octahedral) — three trans pairs
///   `[n0, n1, n2, n3, n4, n5]`, `n2i` opposite `n2i+1`; reference
///   `(n0 − n1) · ((n2 − n3) × (n4 − n5))`.
/// - [`TrigonalPrismatic`](StereoKind::TrigonalPrismatic) — one triangular face,
///   then the other, `[t0, t1, t2, b0, b1, b2]`, `bi` eclipsing `ti`; reference
///   `(t0 − b0) · ((t1 − t0) × (t2 − t0))`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoConfiguration {
    locus: StereoLocus,
    kind: StereoKind,
    neighbors: Vec<SiteId>,
}

impl StereoConfiguration {
    /// Builds a configuration at `locus` of `kind` from its `neighbors`, given in the
    /// order the kind's contract prescribes.
    ///
    /// Returns `None` unless `locus` [anchors](StereoLocus::anchors) `kind` and its
    /// neighbour count is the kind's [`StereoKind::slot_count`].
    #[inline]
    pub fn new(
        locus: StereoLocus,
        kind: StereoKind,
        neighbors: impl IntoIterator<Item = SiteId>,
    ) -> Option<Self> {
        let neighbors: Vec<SiteId> = neighbors.into_iter().collect();
        (locus.anchors(kind) && neighbors.len() == kind.slot_count()).then_some(Self {
            locus,
            kind,
            neighbors,
        })
    }

    /// The anchor of the stereogenic unit.
    #[inline]
    pub fn locus(&self) -> StereoLocus {
        self.locus
    }

    /// The idealised geometry fixing which orderings are equivalent.
    #[inline]
    pub fn kind(&self) -> StereoKind {
        self.kind
    }

    /// The neighbours in reference order (see the type's contract).
    #[inline]
    pub fn neighbors(&self) -> &[SiteId] {
        &self.neighbors
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::SiteId;

    use crate::BondId;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn neighbors(count: usize) -> Vec<SiteId> {
        (1..=count as u32).map(s).collect()
    }

    #[test]
    fn slot_count_is_the_neighbor_count_of_the_geometry() {
        assert_eq!(StereoKind::Tetrahedral.slot_count(), 4);
        assert_eq!(StereoKind::CisTrans.slot_count(), 4);
        assert_eq!(StereoKind::Allene.slot_count(), 4);
        assert_eq!(StereoKind::SquarePlanar.slot_count(), 4);
        assert_eq!(StereoKind::TrigonalBipyramidal.slot_count(), 5);
        assert_eq!(StereoKind::SquarePyramidal.slot_count(), 5);
        assert_eq!(StereoKind::Octahedral.slot_count(), 6);
        assert_eq!(StereoKind::TrigonalPrismatic.slot_count(), 6);
    }

    #[test]
    fn configuration_count_is_the_number_of_distinct_stereoisomers() {
        assert_eq!(StereoKind::Tetrahedral.configuration_count(), 2);
        assert_eq!(StereoKind::CisTrans.configuration_count(), 2);
        assert_eq!(StereoKind::Allene.configuration_count(), 2);
        assert_eq!(StereoKind::SquarePlanar.configuration_count(), 3);
        assert_eq!(StereoKind::TrigonalBipyramidal.configuration_count(), 20);
        assert_eq!(StereoKind::SquarePyramidal.configuration_count(), 30);
        assert_eq!(StereoKind::Octahedral.configuration_count(), 30);
        assert_eq!(StereoKind::TrigonalPrismatic.configuration_count(), 120);
    }

    #[test]
    fn chiral_geometries_are_chiral() {
        assert!(StereoKind::Tetrahedral.is_chiral());
        assert!(StereoKind::Allene.is_chiral());
        assert!(StereoKind::TrigonalBipyramidal.is_chiral());
        assert!(StereoKind::SquarePyramidal.is_chiral());
        assert!(StereoKind::Octahedral.is_chiral());
        assert!(StereoKind::TrigonalPrismatic.is_chiral());
    }

    #[test]
    fn a_site_anchors_every_coordination_center() {
        for kind in [
            StereoKind::Tetrahedral,
            StereoKind::SquarePlanar,
            StereoKind::TrigonalBipyramidal,
            StereoKind::SquarePyramidal,
            StereoKind::Octahedral,
            StereoKind::TrigonalPrismatic,
        ] {
            assert!(StereoLocus::Site(s(1)).anchors(kind), "{kind:?}");
        }
    }

    #[test]
    fn a_bond_anchors_a_double_bond() {
        assert!(StereoLocus::Bond(b(1)).anchors(StereoKind::CisTrans));
    }

    #[test]
    fn an_axis_anchors_an_allene() {
        assert!(StereoLocus::Axis(s(1)).anchors(StereoKind::Allene));
    }

    #[test]
    fn locus_returns_the_anchor() {
        let config = StereoConfiguration::new(
            StereoLocus::Site(s(7)),
            StereoKind::Tetrahedral,
            neighbors(4),
        )
        .unwrap();
        assert_eq!(config.locus(), StereoLocus::Site(s(7)));
    }

    #[test]
    fn kind_returns_the_geometry() {
        let config = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            neighbors(4),
        )
        .unwrap();
        assert_eq!(config.kind(), StereoKind::Tetrahedral);
    }

    #[test]
    fn neighbors_are_returned_in_the_given_order() {
        let order = [s(4), s(2), s(3), s(1)];
        let config =
            StereoConfiguration::new(StereoLocus::Site(s(9)), StereoKind::Tetrahedral, order)
                .unwrap();
        assert_eq!(config.neighbors(), order.as_slice());
    }

    #[test]
    fn cis_trans_and_square_planar_are_not_chiral() {
        assert!(!StereoKind::CisTrans.is_chiral());
        assert!(!StereoKind::SquarePlanar.is_chiral());
    }

    #[test]
    fn only_a_site_anchors_a_coordination_center() {
        assert!(!StereoLocus::Bond(b(1)).anchors(StereoKind::Tetrahedral));
        assert!(!StereoLocus::Axis(s(1)).anchors(StereoKind::Tetrahedral));
    }

    #[test]
    fn only_a_bond_anchors_a_double_bond() {
        assert!(!StereoLocus::Site(s(1)).anchors(StereoKind::CisTrans));
        assert!(!StereoLocus::Axis(s(1)).anchors(StereoKind::CisTrans));
    }

    #[test]
    fn only_an_axis_anchors_an_allene() {
        assert!(!StereoLocus::Site(s(1)).anchors(StereoKind::Allene));
        assert!(!StereoLocus::Bond(b(1)).anchors(StereoKind::Allene));
    }

    #[test]
    fn new_rejects_a_kind_off_its_anchor() {
        assert!(
            StereoConfiguration::new(
                StereoLocus::Bond(b(1)),
                StereoKind::Tetrahedral,
                neighbors(4)
            )
            .is_none()
        );
    }

    #[test]
    fn new_rejects_an_empty_neighbor_list() {
        assert!(
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::Tetrahedral,
                neighbors(0)
            )
            .is_none()
        );
    }

    #[test]
    fn new_rejects_more_neighbors_than_the_slot_count() {
        assert!(
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::Tetrahedral,
                neighbors(5)
            )
            .is_none()
        );
    }

    #[test]
    fn loci_order_by_anchor_then_by_identifier() {
        assert!(StereoLocus::Site(s(9)) < StereoLocus::Bond(b(1)));
        assert!(StereoLocus::Bond(b(9)) < StereoLocus::Axis(s(1)));
        assert!(StereoLocus::Site(s(1)) < StereoLocus::Site(s(2)));
    }

    #[test]
    fn kinds_order_by_neighbor_count_then_configuration_count_then_locus() {
        assert!(StereoKind::SquarePlanar < StereoKind::TrigonalBipyramidal);
        assert!(StereoKind::SquarePyramidal < StereoKind::Octahedral);
        assert!(StereoKind::TrigonalBipyramidal < StereoKind::SquarePyramidal);
        assert!(StereoKind::Octahedral < StereoKind::TrigonalPrismatic);
        assert!(StereoKind::Tetrahedral < StereoKind::CisTrans);
        assert!(StereoKind::CisTrans < StereoKind::Allene);
    }

    #[test]
    fn configurations_are_equal_exactly_when_their_parts_match() {
        let base = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            neighbors(4),
        )
        .unwrap();
        assert_eq!(
            base,
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::Tetrahedral,
                neighbors(4),
            )
            .unwrap(),
        );
        assert_ne!(
            base,
            StereoConfiguration::new(
                StereoLocus::Site(s(2)),
                StereoKind::Tetrahedral,
                neighbors(4),
            )
            .unwrap(),
        );
        assert_ne!(
            base,
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::SquarePlanar,
                neighbors(4),
            )
            .unwrap(),
        );
        assert_ne!(
            base,
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::Tetrahedral,
                [s(1), s(2), s(3), s(5)],
            )
            .unwrap(),
        );
    }

    #[test]
    fn reordering_the_neighbors_yields_a_different_configuration() {
        let forward = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(1), s(2), s(3), s(4)],
        )
        .unwrap();
        let swapped = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(2), s(1), s(3), s(4)],
        )
        .unwrap();
        assert_ne!(forward, swapped);
    }
}

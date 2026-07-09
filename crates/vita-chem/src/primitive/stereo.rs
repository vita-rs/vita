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
    /// Returns `None` unless `locus` is the anchor `kind` lives on — a site for a
    /// centre, a bond for a double bond, an axis for an allene — and the neighbour
    /// count is the kind's [`StereoKind::slot_count`].
    #[inline]
    pub fn new(
        locus: StereoLocus,
        kind: StereoKind,
        neighbors: impl IntoIterator<Item = SiteId>,
    ) -> Option<Self> {
        let neighbors: Vec<SiteId> = neighbors.into_iter().collect();
        let anchored = matches!(
            (locus, kind),
            (
                StereoLocus::Site(_),
                StereoKind::Tetrahedral
                    | StereoKind::SquarePlanar
                    | StereoKind::TrigonalBipyramidal
                    | StereoKind::SquarePyramidal
                    | StereoKind::Octahedral
                    | StereoKind::TrigonalPrismatic,
            ) | (StereoLocus::Bond(_), StereoKind::CisTrans)
                | (StereoLocus::Axis(_), StereoKind::Allene)
        );
        (anchored && neighbors.len() == kind.slot_count()).then_some(Self {
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

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

use vita_core::HasElements;

use super::CoordinationGeometries;
use crate::algorithm::valence::steric_numbers;
use crate::{CoordinationGeometry, HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// Perceives the coordination geometry of each of a molecule's sites.
///
/// A site's electron domains take the arrangement that separates them furthest — the
/// parent polyhedron of its [`steric_numbers`] count — and the substituents take
/// whichever of its vertices the lone pairs leave, the lone pairs claiming the roomiest
/// first. Naming what remains is VSEPR, and this is its whole content: a site is read
/// through its domain count and its coordination number, and nothing else about it
/// enters.
///
/// The model reaches the eleven geometries it names; the two a coordination sphere can
/// take without a domain count predicting them —
/// [`PyramidalizedSquare`](CoordinationGeometry::PyramidalizedSquare) and
/// [`TrigonalPrismatic`](CoordinationGeometry::TrigonalPrismatic) — never come out of it.
/// A site is left unanswered where the model says nothing: where [`steric_numbers`]
/// settles no count, where one substituent or none fixes no arrangement, or where the
/// domains outrun the vocabulary.
///
/// # Complexity
///
/// O((V + E) · log (V + E)) time and O(V + E) space, over the molecule's `V` sites and
/// `E` bonds, assuming [`degree`](crate::HasBonds::degree) runs in O(degree); counting the
/// domains dominates.
pub fn perceive<M>(mol: &M) -> CoordinationGeometries
where
    M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons,
{
    CoordinationGeometries::from_pairs(
        steric_numbers(mol)
            .iter()
            .filter_map(|(site, domains)| Some((site, vsepr(domains, mol.degree(site) as u32)?))),
    )
}

/// The geometry VSEPR gives a site holding `domains` electron domains of which
/// `coordination` are substituents, or `None` where the model names none.
///
/// Reading down a column is one parent arrangement losing a vertex to each further lone
/// pair: the trigonal bipyramid empties its equatorial belt before its poles, the
/// octahedron its trans pairs.
const fn vsepr(domains: u32, coordination: u32) -> Option<CoordinationGeometry> {
    match (domains, coordination) {
        (2, 2) => Some(CoordinationGeometry::Linear),
        (3, 3) => Some(CoordinationGeometry::TrigonalPlanar),
        (3, 2) => Some(CoordinationGeometry::Angular),
        (4, 4) => Some(CoordinationGeometry::Tetrahedral),
        (4, 3) => Some(CoordinationGeometry::TrigonalPyramidal),
        (4, 2) => Some(CoordinationGeometry::Angular),
        (5, 5) => Some(CoordinationGeometry::TrigonalBipyramidal),
        (5, 4) => Some(CoordinationGeometry::Seesaw),
        (5, 3) => Some(CoordinationGeometry::TShaped),
        (5, 2) => Some(CoordinationGeometry::Linear),
        (6, 6) => Some(CoordinationGeometry::Octahedral),
        (6, 5) => Some(CoordinationGeometry::SquarePyramidal),
        (6, 4) => Some(CoordinationGeometry::SquarePlanar),
        (6, 3) => Some(CoordinationGeometry::TShaped),
        (6, 2) => Some(CoordinationGeometry::Linear),
        _ => None,
    }
}

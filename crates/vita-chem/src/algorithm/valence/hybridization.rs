use vita_core::{HasElements, SiteId};

use super::lone_pairs::lone_pairs;
use crate::{HasBondOrders, HasFormalCharges, HasRadicalElectrons, Hybridization};

/// Hybridization of `site` from its electron-domain count.
///
/// Counts electron domains — bonded neighbours plus [`lone_pairs`] — and names
/// the geometry: one or none is [`S`](Hybridization::S), two
/// [`Sp`](Hybridization::Sp), three [`Sp2`](Hybridization::Sp2), four
/// [`Sp3`](Hybridization::Sp3), five [`Sp3d`](Hybridization::Sp3d), six
/// [`Sp3d2`](Hybridization::Sp3d2), seven [`Sp3d3`](Hybridization::Sp3d3), and
/// eight or more [`Other`](Hybridization::Other).
///
/// This is the steric (VSEPR) count, blind to conjugation and to coordinates:
/// an amide nitrogen is [`Sp3`](Hybridization::Sp3) rather than `Sp2`, and
/// square-planar [`Sp2d`](Hybridization::Sp2d) is not told apart from
/// [`Sp3`](Hybridization::Sp3).
///
/// Returns `None` exactly when [`lone_pairs`] does — a d-/f-block element, an
/// aromatic (delocalised) bond, or an impossible valence — leaving the domain
/// count undefined.
///
/// # Complexity
///
/// O(degree) time.
pub fn hybridization<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>(
    mol: &M,
    site: SiteId,
) -> Option<Hybridization> {
    let domains = mol.degree(site) as u32 + lone_pairs(mol, site)?;
    Some(match domains {
        0 | 1 => Hybridization::S,
        2 => Hybridization::Sp,
        3 => Hybridization::Sp2,
        4 => Hybridization::Sp3,
        5 => Hybridization::Sp3d,
        6 => Hybridization::Sp3d2,
        7 => Hybridization::Sp3d3,
        _ => Hybridization::Other,
    })
}

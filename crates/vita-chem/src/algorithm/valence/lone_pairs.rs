use vita_core::{Element, HasElements, SiteId};

use super::explicit::valence;
use crate::{HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// Number of lone (non-bonding) electron pairs on `site`.
///
/// What is left of the element's valence electrons once the bonding and
/// unpaired ones are removed: `(valence electrons − formal charge − bond-order
/// sum − radicals) / 2`. Oxygen in water has two; the ammonium nitrogen has
/// none.
///
/// Returns `None` when no exact count exists:
/// - `site` holds a d- or f-block element, whose valence-electron count is
///   ambiguous;
/// - an incident bond is aromatic, so the bonding electrons are not localised
///   (see [`valence`] — kekulise the ring first);
/// - the bookkeeping is negative, describing an impossible structure.
///
/// # Complexity
///
/// O(degree) time.
pub fn lone_pairs<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>(
    mol: &M,
    site: SiteId,
) -> Option<u32> {
    let electrons = valence_electrons(mol.element(site))? as i32;
    let bonding = valence(mol, site)? as i32;
    let charge = mol.formal_charge(site) as i32;
    let radicals = mol.radical_electron(site) as i32;
    let free = electrons - charge - bonding - radicals;
    if free < 0 {
        return None;
    }
    Some(free as u32 / 2)
}

/// Valence (outer-shell s and p) electron count of a main-group element.
///
/// Returns `None` for the d- and f-block, where the count is ambiguous.
fn valence_electrons(element: Element) -> Option<u8> {
    Some(match element.atomic_number() {
        1 | 3 | 11 | 19 | 37 | 55 | 87 => 1,
        2 | 4 | 12 | 20 | 38 | 56 | 88 => 2,
        5 | 13 | 31 | 49 | 81 | 113 => 3,
        6 | 14 | 32 | 50 | 82 | 114 => 4,
        7 | 15 | 33 | 51 | 83 | 115 => 5,
        8 | 16 | 34 | 52 | 84 | 116 => 6,
        9 | 17 | 35 | 53 | 85 | 117 => 7,
        10 | 18 | 36 | 54 | 86 | 118 => 8,
        _ => return None,
    })
}

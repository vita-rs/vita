use std::cmp::Ordering;

use vita_core::{HasElements, SiteId};

use crate::algorithm::utils::electronegativity_rank;
use crate::{BondOrder, HasBondOrders, HasFormalCharges};

/// Oxidation state of `site`: the charge left once every bond's electrons are
/// awarded wholly to the bond's more electronegative end.
///
/// The ionic approximation of the IUPAC 2016 recommendation: relative to the
/// formal charge, each bond adds its integer order when the partner outranks
/// `site` on the Allen electronegativity order (Pauling across the f-block),
/// subtracts it when `site` outranks the partner, and moves nothing between
/// equal ranks — homonuclear bonds above all. Methane's carbon is −4, the
/// carbon of CO₂ is +4, and the d-block answers where
/// [`lone_pairs`](super::lone_pairs) cannot. The recommendation's caveat — a
/// reversibly bonded Lewis-acid ligand keeps its donated pair — rests on a
/// reversibility the bond graph does not record, and is not applied.
///
/// Returns `None` when no exact state exists:
/// - an incident bond is [`Aromatic`](BondOrder::Aromatic), its electrons not
///   localised into an integer order (see [`valence`](super::valence) —
///   kekulise the ring first);
/// - an incident bond reaches an element beyond lawrencium, past every
///   measured electronegativity;
/// - the total leaves `i8`, describing an impossible structure.
///
/// # Complexity
///
/// O(d) time and O(1) space, where `d` is the degree of `site`, assuming
/// [`bonds_of`](crate::HasBonds::bonds_of) runs in O(degree).
pub fn oxidation_state<M: HasBondOrders + HasElements + HasFormalCharges>(
    mol: &M,
    site: SiteId,
) -> Option<i8> {
    let own = electronegativity_rank(mol.element(site));
    let mut state = i32::from(mol.formal_charge(site));
    for (bond, partner) in mol.bonds_of(site) {
        let order: i32 = match mol.bond_order(bond) {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Quintuple => 5,
            BondOrder::Hextuple => 6,
            BondOrder::Aromatic => return None,
        };
        state += order
            * match own?.cmp(&electronegativity_rank(mol.element(partner))?) {
                Ordering::Less => 1,
                Ordering::Equal => 0,
                Ordering::Greater => -1,
            };
    }
    i8::try_from(state).ok()
}

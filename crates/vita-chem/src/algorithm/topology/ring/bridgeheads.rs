use std::collections::{HashMap, HashSet};

use vita_core::{HasSites, SiteId};

use crate::{BondId, HasBonds};

/// The bridgehead atoms of a molecule.
///
/// A bridgehead is an atom where a bridge meets a ring: an endpoint of the
/// bonds two rings share, when they share two or more. Ortho-fused rings
/// (sharing one bond) and spiro rings (sharing one site) have none; a bridged
/// system such as bicyclo[2.2.2]octane has two.
///
/// # Complexity
///
/// O(V² · E) time, dominated by the minimum cycle basis it builds.
pub fn bridgeheads<M: HasBonds + HasSites>(mol: &M) -> impl Iterator<Item = SiteId> {
    let basis = super::rings(mol);
    let cycles: Vec<HashSet<BondId>> = basis
        .iter()
        .map(|r| r.bonds().iter().copied().collect())
        .collect();

    let mut heads: HashSet<SiteId> = HashSet::new();
    for (i, a) in cycles.iter().enumerate() {
        for b in &cycles[i + 1..] {
            let shared: Vec<BondId> = a.intersection(b).copied().collect();
            if shared.len() < 2 {
                continue;
            }
            let mut incident: HashMap<SiteId, usize> = HashMap::new();
            for bond in shared {
                let (u, v) = mol.bond_endpoints(bond);
                *incident.entry(u).or_default() += 1;
                *incident.entry(v).or_default() += 1;
            }
            heads.extend(
                incident
                    .into_iter()
                    .filter(|&(_, n)| n == 1)
                    .map(|(s, _)| s),
            );
        }
    }

    let mut heads: Vec<SiteId> = heads.into_iter().collect();
    heads.sort_unstable();
    heads.into_iter()
}

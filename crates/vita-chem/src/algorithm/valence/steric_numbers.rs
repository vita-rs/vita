use vita_core::{HasElements, SiteId};

use super::lone_pairs;
use crate::algorithm::conjugation::systems;
use crate::algorithm::utils::SortedMap;
use crate::{HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// How many electron domains lie about each of a molecule's sites.
///
/// The number counts a site's bonds, one apiece whatever their order, plus the
/// lone pairs it keeps localized. A site whose arithmetic settles no exact
/// count goes uncounted.
///
/// Obtain via [`steric_numbers`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StericNumbers {
    counts: SortedMap<SiteId, u32>,
}

impl StericNumbers {
    /// Number of counted sites.
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns `true` if no site is counted.
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Returns the number of electron domains about `site`.
    ///
    /// Returns `None` if `site` is absent from the molecule or goes uncounted.
    pub fn get(&self, site: SiteId) -> Option<u32> {
        self.counts.get(&site).copied()
    }

    /// Iterates the `(site, steric number)` pairs, ordered by site.
    pub fn iter(&self) -> impl Iterator<Item = (SiteId, u32)> + '_ {
        self.counts.iter().map(|(&site, &count)| (site, count))
    }
}

/// The steric number of each of a molecule's sites.
///
/// A site's domains are the regions its valence electrons occupy: one per
/// bond, whatever its order — a multiple bond's π component lies along the
/// σ-bond already counted — and one per lone pair the site keeps localized. A
/// pair donated into a conjugated π network is not localized; it spreads along
/// the network and adds no domain. Discounting it leaves the total independent
/// of the Lewis form drawn: an amide nitrogen has three domains whether the
/// neutral form donates from it or the zwitterionic form donates from the
/// oxygen.
///
/// A site goes uncounted where [`lone_pairs`] settles no exact count — a d- or
/// f-block element, an incident aromatic bond, or arithmetic describing an
/// impossible structure.
///
/// # Complexity
///
/// O((V + E) · log (V + E)) time and O(V + E) space, over the molecule's `V`
/// sites and `E` bonds, assuming [`bonds_of`](crate::HasBonds::bonds_of) and
/// [`degree`](crate::HasBonds::degree) run in O(degree); perceiving the
/// conjugated systems dominates.
pub fn steric_numbers<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>(
    mol: &M,
) -> StericNumbers {
    let conjugation = systems(mol);
    StericNumbers {
        counts: SortedMap::from_pairs(mol.sites().filter_map(|site| {
            let pairs = lone_pairs(mol, site)?;
            let donated: u32 = conjugation
                .of_site(site)
                .map(|system| system.donated_pairs(site))
                .sum();
            let kept = pairs
                .checked_sub(donated)
                .expect("a system draws its donations from the donor's own lone pairs");
            Some((site, mol.degree(site) as u32 + kept))
        })),
    }
}

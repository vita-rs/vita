mod bridgeheads;
mod count;
mod families;
mod membership;
mod rings;

pub use bridgeheads::{Bridgeheads, bridgeheads};
pub use count::count;
pub use families::{RingFamilies, RingFamily, families};
pub use membership::{RingMembership, membership};
pub use rings::{Ring, Rings, rings};

use vita_core::SiteId;

/// A ring system: a maximal group of rings connected through shared sites, as
/// its set of sites.
///
/// Fused, bridged, and spiro rings all coalesce into one system.
///
/// Obtain from [`RingSystems`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingSystem {
    sites: Vec<SiteId>,
}

impl RingSystem {
    /// Number of sites in the ring system.
    pub fn len(&self) -> usize {
        self.sites.len()
    }

    /// Returns `true` if the ring system contains no sites. Always `false` — a
    /// ring system is non-empty — but provided alongside [`len`](Self::len).
    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }

    /// Returns `true` if `site` lies in this ring system.
    pub fn contains(&self, site: SiteId) -> bool {
        self.sites.binary_search(&site).is_ok()
    }

    /// Iterates the ring system's sites in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.sites.iter().copied()
    }
}

/// The ring systems of a molecule.
///
/// Each is a [`RingSystem`]; the systems are ordered by their sites. An acyclic
/// molecule has none.
///
/// Obtain via [`Rings::systems`] or [`RingFamilies::systems`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingSystems {
    systems: Vec<RingSystem>,
}

impl RingSystems {
    /// Wraps ascending groups of sites, each an already-sorted ring system.
    fn new(groups: Vec<Vec<SiteId>>) -> Self {
        RingSystems {
            systems: groups
                .into_iter()
                .map(|sites| RingSystem { sites })
                .collect(),
        }
    }

    /// Number of ring systems.
    pub fn len(&self) -> usize {
        self.systems.len()
    }

    /// Returns `true` if the molecule has no ring systems.
    pub fn is_empty(&self) -> bool {
        self.systems.is_empty()
    }

    /// Iterates the ring systems, ordered by their sites.
    pub fn iter(&self) -> impl Iterator<Item = &RingSystem> + '_ {
        self.systems.iter()
    }
}

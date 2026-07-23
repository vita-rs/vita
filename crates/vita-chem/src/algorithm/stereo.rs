//! Stereochemistry over declared configurations, which never name a
//! handedness.
//!
//! A configuration is a coset: an ordered neighbour list under the symmetry
//! group of its geometry. On that footing [`perceive`] reads
//! [`StereoConfigurations`] off coordinates, viewable as
//! [`WithStereoConfigurations`]; [`stereocenters`] finds the
//! [`Stereocenters`] where more than one configuration is realisable;
//! [`form`] canonicalizes constitution and stereochemistry into one
//! [`StereoForm`] and names the [`StereoRelationship`] between two
//! molecules; [`stereoisomers`] enumerates the [`Stereoisomers`] no
//! symmetry equates; [`consistency`] reconciles declarations against
//! stereogenic loci — unspecified and overspecified alike — as a
//! [`StereoConsistency`].

mod center;
mod enumerate;
mod identity;
mod perceive;
mod validate;

pub use center::{Stereocenters, stereocenters};
pub use enumerate::{Stereoisomers, stereoisomers};
pub use identity::{StereoForm, StereoRelationship, form};
pub use perceive::{StereoConfigurations, WithStereoConfigurations, perceive};
pub use validate::{StereoConsistency, consistency};

use vita_core::SiteId;

use crate::algorithm::canonical::{Canonical, canonicalize};
use crate::algorithm::utils::{FxHashMap, FxHashSet};
use crate::{
    BondId, BondOrder, HasBondOrders, HasBonds, HasStereoConfigurations, StereoConfiguration,
    StereoDescriptor, StereoKind, StereoLocus,
};

/// A stereogenic frame located in the graph: the atoms that pin it, and the
/// substituents its geometry arranges.
///
/// A site's frame is the atom and its neighbours. An edge's or axis's is the two
/// termini of its rigid double-bond chain and the two substituents each bears,
/// walked out from the anchor — so a plain double bond and a long cumulene resolve
/// alike. The substituents of a two-ended frame are grouped by end.
struct Frame {
    anchors: Vec<SiteId>,
    substituents: Vec<SiteId>,
}

/// Locates the stereogenic frame a `locus` names, or `None` if the graph does not
/// realise one — a branched cumulene, or a terminus without its two substituents.
///
/// Stereochemistry across a bond or axis is a rigidity phenomenon: the frame is the
/// maximal chain of cumulated double bonds, so its termini — where the substituents
/// hang — are found by following those bonds, not the graph's plain connectivity.
fn frame<M: HasBondOrders>(mol: &M, locus: StereoLocus) -> Option<Frame> {
    match locus {
        StereoLocus::Site(site) => Some(Frame {
            anchors: vec![site],
            substituents: mol.neighbors(site).collect(),
        }),
        StereoLocus::Bond(bond) => {
            let (first, second) = mol.bond_endpoints(bond);
            poles(mol, walk(mol, first, second)?, walk(mol, second, first)?)
        }
        StereoLocus::Axis(center) => {
            let mut chain = double_neighbours(mol, center, center);
            let (Some(first), Some(second), None) = (chain.next(), chain.next(), chain.next())
            else {
                return None;
            };
            poles(mol, walk(mol, first, center)?, walk(mol, second, center)?)
        }
    }
}

/// The frame of an edge or axis from its two termini: each terminus is an anchor and
/// bears two substituents, grouped by end — `None` unless each bears exactly two.
fn poles<M: HasBondOrders>(mol: &M, first: SiteId, second: SiteId) -> Option<Frame> {
    let first_subs = terminal_substituents(mol, first);
    let second_subs = terminal_substituents(mol, second);
    if first_subs.len() != 2 || second_subs.len() != 2 {
        return None;
    }
    Some(Frame {
        anchors: vec![first, second],
        substituents: vec![first_subs[0], first_subs[1], second_subs[0], second_subs[1]],
    })
}

/// Follows the double-bond chain from `start`, entered from `came_from`, to its
/// terminus — the atom with no onward double bond — or `None` if the chain branches.
fn walk<M: HasBondOrders>(mol: &M, start: SiteId, came_from: SiteId) -> Option<SiteId> {
    let (mut previous, mut current) = (came_from, start);
    loop {
        let mut chain = double_neighbours(mol, current, previous);
        match (chain.next(), chain.next()) {
            (None, _) => return Some(current),
            (Some(next), None) => (previous, current) = (current, next),
            (Some(_), Some(_)) => return None,
        }
    }
}

/// The double-bonded neighbours of `site`, excluding `exclude`.
fn double_neighbours<M: HasBondOrders>(
    mol: &M,
    site: SiteId,
    exclude: SiteId,
) -> impl Iterator<Item = SiteId> + '_ {
    mol.bonds_of(site).filter_map(move |(bond, other)| {
        (other != exclude && mol.bond_order(bond) == BondOrder::Double).then_some(other)
    })
}

/// The substituents of a chain terminus: its neighbours off the double-bond chain.
fn terminal_substituents<M: HasBondOrders>(mol: &M, terminus: SiteId) -> Vec<SiteId> {
    mol.bonds_of(terminus)
        .filter_map(|(bond, other)| (mol.bond_order(bond) != BondOrder::Double).then_some(other))
        .collect()
}

/// Every locus a molecule could bear stereochemistry at — each site as a coordination
/// centre and as an allene axis, each bond as a double bond. The caller's `candidate`
/// then says which the geometry admits.
fn candidate_loci<M: HasBonds>(mol: &M) -> impl Iterator<Item = StereoLocus> + '_ {
    mol.sites()
        .map(StereoLocus::Site)
        .chain(mol.sites().map(StereoLocus::Axis))
        .chain(mol.bonds().map(StereoLocus::Bond))
}

/// The stereo an atom carries under a labeling: the sorted descriptors of the
/// configurations anchored on it, a coloring that refines the labeling and appears
/// in the canonical form.
type Signal = Vec<StereoDescriptor>;

/// A canonical labeling refined by stereochemistry: the graph and the caller's
/// coloring, further split by the configurations laid over it.
type Configured<VK, EK> = Canonical<(VK, Signal), EK>;

/// A site's symmetry class under a labeling: the least canonical rank in its orbit,
/// so interchangeable atoms share a class.
fn class<VK, EK>(canon: &Canonical<VK, EK>, site: SiteId) -> usize {
    canon
        .orbit(site)
        .expect("the site is in the molecule")
        .iter()
        .map(|member| {
            canon
                .rank(member)
                .expect("an orbit member is in the molecule")
        })
        .min()
        .expect("an orbit is non-empty")
}

/// The stereo signal every anchor carries under a labeling: each configuration's
/// descriptor relative to the current symmetry classes, mirrored when `reflect`,
/// filed against the atoms its locus pins.
fn signals<M, VK, EK>(
    mol: &M,
    canon: &Configured<VK, EK>,
    reflect: bool,
) -> FxHashMap<SiteId, Signal>
where
    M: HasStereoConfigurations,
{
    let mut signal: FxHashMap<SiteId, Signal> = FxHashMap::default();
    let mut file = |anchor: SiteId, descriptor: StereoDescriptor| {
        signal.entry(anchor).or_default().push(descriptor);
    };
    for config in mol.stereo_configurations() {
        let descriptor = config.descriptor(|site| class(canon, site));
        let descriptor = if reflect {
            descriptor.mirror()
        } else {
            descriptor
        };
        match config.locus() {
            StereoLocus::Site(site) | StereoLocus::Axis(site) => file(site, descriptor),
            StereoLocus::Bond(bond) => {
                let (a, b) = mol.bond_endpoints(bond);
                file(a, descriptor);
                file(b, descriptor);
            }
        }
    }
    for descriptors in signal.values_mut() {
        descriptors.sort_unstable();
    }
    signal
}

/// The caller's colouring `site_key` refined by a stereo `signal` — each site keyed by
/// its own colour and the descriptors filed on it, so stereochemistry splits the
/// classes a bare constitution leaves merged.
fn refined_key<'a, VK>(
    signal: &'a FxHashMap<SiteId, Signal>,
    site_key: &'a impl Fn(SiteId) -> VK,
) -> impl Fn(SiteId) -> (VK, Signal) + 'a {
    move |site| {
        (
            site_key(site),
            signal.get(&site).cloned().unwrap_or_default(),
        )
    }
}

/// A canonical labeling with stereochemistry folded in, refined to a fixpoint.
///
/// Each configuration reduces to a descriptor relative to the current symmetry
/// classes, the descriptors colour their anchors, and the labeling is recomputed
/// until the classes settle — a fixpoint that resolves even a pseudo-asymmetric
/// centre, whose stereogenicity turns on the configurations of its neighbours.
/// `reflect` folds in each configuration's mirror image instead, giving the
/// enantiomer's labeling.
fn refined<M, VK, EK>(
    mol: &M,
    site_key: &impl Fn(SiteId) -> VK,
    bond_key: &impl Fn(BondId) -> EK,
    reflect: bool,
) -> Configured<VK, EK>
where
    M: HasStereoConfigurations,
    VK: Ord,
    EK: Ord,
{
    let mut signal: FxHashMap<SiteId, Signal> = FxHashMap::default();
    let mut classes = 0;
    loop {
        let canon = canonicalize(mol, refined_key(&signal, site_key), bond_key);
        let count = canon.orbits().count();
        if count == classes {
            return canon;
        }
        classes = count;
        signal = signals(mol, &canon, reflect);
    }
}

/// The settled stereo signal: each atom's descriptors under the symmetry classes the
/// stereochemistry itself refines to. This is the coloring a stereogenic-unit
/// detection must use, so a pseudo-asymmetric centre's inequivalent neighbours are
/// told apart.
fn settle<M, VK, EK>(
    mol: &M,
    site_key: &impl Fn(SiteId) -> VK,
    bond_key: &impl Fn(BondId) -> EK,
) -> FxHashMap<SiteId, Signal>
where
    M: HasStereoConfigurations,
    VK: Ord,
    EK: Ord,
{
    signals(mol, &refined(mol, site_key, bond_key, false), false)
}

/// The distinct configurations `kind` realises at `locus` over `subs`, whose symmetry
/// classes are given by `class`.
///
/// Enumerates the neighbour orderings the geometry admits and keeps one per distinct
/// [`StereoDescriptor`] under `class` — collapsing both the rotations a configuration
/// treats as equivalent and the swaps of interchangeable substituents.
fn configurations(
    locus: StereoLocus,
    kind: StereoKind,
    subs: &[SiteId],
    class: impl Fn(SiteId) -> usize,
) -> Vec<StereoConfiguration> {
    let mut seen: FxHashSet<StereoDescriptor> = FxHashSet::default();
    let mut result = Vec::new();
    for ordering in presentations(kind, subs) {
        let Some(config) = StereoConfiguration::new(locus, kind, ordering) else {
            continue;
        };
        if seen.insert(config.descriptor(&class)) {
            result.push(config);
        }
    }
    result
}

/// The distinct configurations `locus` of `kind` realises in a molecule, under a
/// coloring `site_key` that already carries any stereo refinement.
///
/// Locates and individualises the frame, colours the graph so interchangeable
/// substituents share a class, and returns one configuration per stereoisomer those
/// classes admit — empty if the graph realises no frame of the right size. The locus
/// is stereogenic exactly when more than one results; an enumeration branches over
/// each.
fn realisable<M, VK, EK>(
    mol: &M,
    locus: StereoLocus,
    kind: StereoKind,
    site_key: &impl Fn(SiteId) -> VK,
    bond_key: &impl Fn(BondId) -> EK,
) -> Vec<StereoConfiguration>
where
    M: HasBondOrders,
    VK: Ord,
    EK: Ord,
{
    if !locus.anchors(kind) {
        return Vec::new();
    }
    let Some(frame) = frame(mol, locus) else {
        return Vec::new();
    };
    if frame.substituents.len() != kind.slot_count() {
        return Vec::new();
    }
    let canon = canonicalize(
        mol,
        |site| (frame_mark(site, &frame.anchors), site_key(site)),
        bond_key,
    );
    configurations(locus, kind, &frame.substituents, |site| class(&canon, site))
}

/// A distinct individualisation mark for each anchor of a frame — its position in
/// `anchors` plus one, or zero for a site outside it — so a canonical labeling holds
/// the frame fixed and tells its members apart. A site marks itself, an edge or axis
/// its two termini.
fn frame_mark(site: SiteId, anchors: &[SiteId]) -> u8 {
    anchors
        .iter()
        .position(|&anchor| anchor == site)
        .map_or(0, |index| index as u8 + 1)
}

/// The neighbour orderings `kind` admits over `subs`: every ordering for a centre,
/// only the within-end orderings for a two-ended edge — a substituent cannot cross
/// from one end to the other.
fn presentations(kind: StereoKind, subs: &[SiteId]) -> Vec<Vec<SiteId>> {
    let per_end = subs.len() / kind.ends();
    subs.chunks(per_end)
        .fold(vec![Vec::new()], |orderings, end| {
            let within = permutations(end);
            orderings
                .iter()
                .flat_map(|prefix| {
                    within.iter().map(move |ordering| {
                        let mut next = prefix.clone();
                        next.extend_from_slice(ordering);
                        next
                    })
                })
                .collect()
        })
}

/// Every permutation of `items`.
fn permutations<T: Clone>(items: &[T]) -> Vec<Vec<T>> {
    let mut buffer = items.to_vec();
    let mut result = Vec::new();
    permute(&mut buffer, 0, &mut result);
    result
}

/// Appends to `result` every permutation of `items` that holds its first `start`
/// elements fixed.
fn permute<T: Clone>(items: &mut [T], start: usize, result: &mut Vec<Vec<T>>) {
    if start == items.len() {
        result.push(items.to_vec());
        return;
    }
    for i in start..items.len() {
        items.swap(start, i);
        permute(items, start + 1, result);
        items.swap(start, i);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondOrder::{Double, Single};
    use crate::{BondId, HasBonds};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
        configs: Vec<StereoConfiguration>,
    }

    impl Mol {
        fn with(mut self, configs: impl IntoIterator<Item = StereoConfiguration>) -> Self {
            self.configs = configs.into_iter().collect();
            self
        }
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasBonds for Mol {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.endpoints[i]
        }
    }

    impl HasBondOrders for Mol {
        fn bond_order(&self, bond: BondId) -> BondOrder {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.orders[i]
        }
    }

    impl HasStereoConfigurations for Mol {
        fn stereo_configurations(&self) -> impl Iterator<Item = StereoConfiguration> + '_ {
            self.configs.iter().cloned()
        }
    }

    fn mol(sites: &[u32], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: sites.iter().map(|&id| s(id)).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c, _)| (s(a), s(c))).collect(),
            orders: bonds.iter().map(|&(_, _, _, order)| order).collect(),
            configs: Vec::new(),
        }
    }

    fn center() -> Mol {
        mol(
            &[1, 2, 3, 4, 5],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
            ],
        )
    }

    fn trigonal() -> Mol {
        mol(
            &[1, 2, 3, 4],
            &[(1, 1, 2, Single), (2, 1, 3, Single), (3, 1, 4, Single)],
        )
    }

    fn alkene() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6],
            &[
                (1, 1, 2, Double),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 2, 5, Single),
                (5, 2, 6, Single),
            ],
        )
    }

    fn short_alkene() -> Mol {
        mol(
            &[1, 2, 3, 4, 5],
            &[
                (1, 1, 2, Double),
                (2, 1, 3, Single),
                (3, 2, 4, Single),
                (4, 2, 5, Single),
            ],
        )
    }

    fn butatriene() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6, 7, 8],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Double),
                (3, 3, 4, Double),
                (4, 1, 5, Single),
                (5, 1, 6, Single),
                (6, 4, 7, Single),
                (7, 4, 8, Single),
            ],
        )
    }

    fn branched() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Double),
                (3, 2, 4, Double),
                (4, 1, 5, Single),
                (5, 1, 6, Single),
            ],
        )
    }

    fn allene() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6, 7],
            &[
                (1, 1, 2, Double),
                (2, 2, 3, Double),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
                (5, 3, 6, Single),
                (6, 3, 7, Single),
            ],
        )
    }

    fn pseudo_asymmetric() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6, 7, 8, 9],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 2, 4, Single),
                (4, 2, 5, Single),
                (5, 2, 6, Single),
                (6, 3, 7, Single),
                (7, 3, 8, Single),
                (8, 3, 9, Single),
            ],
        )
        .with([
            StereoConfiguration::new(
                StereoLocus::Site(s(2)),
                StereoKind::Tetrahedral,
                [s(4), s(5), s(6), s(1)],
            )
            .unwrap(),
            StereoConfiguration::new(
                StereoLocus::Site(s(3)),
                StereoKind::Tetrahedral,
                [s(8), s(7), s(9), s(1)],
            )
            .unwrap(),
        ])
    }

    #[test]
    fn a_site_frame_is_the_atom_and_its_neighbors() {
        let located = frame(&center(), StereoLocus::Site(s(1))).unwrap();
        assert_eq!(located.anchors, vec![s(1)]);
        assert_eq!(located.substituents, vec![s(2), s(3), s(4), s(5)]);
    }

    #[test]
    fn a_double_bond_frame_pairs_each_terminus_with_its_substituents() {
        let located = frame(&alkene(), StereoLocus::Bond(b(1))).unwrap();
        assert_eq!(located.anchors, vec![s(1), s(2)]);
        assert_eq!(located.substituents, vec![s(3), s(4), s(5), s(6)]);
    }

    #[test]
    fn a_cumulene_bond_frame_walks_out_to_the_chain_termini() {
        let located = frame(&butatriene(), StereoLocus::Bond(b(2))).unwrap();
        assert_eq!(located.anchors, vec![s(1), s(4)]);
        assert_eq!(located.substituents, vec![s(5), s(6), s(7), s(8)]);
    }

    #[test]
    fn a_double_bond_frame_is_none_when_a_terminus_lacks_two_substituents() {
        assert!(frame(&short_alkene(), StereoLocus::Bond(b(1))).is_none());
    }

    #[test]
    fn a_branched_double_bond_chain_has_no_frame() {
        assert!(frame(&branched(), StereoLocus::Bond(b(1))).is_none());
    }

    #[test]
    fn an_allene_axis_frame_is_the_two_termini_and_their_substituents() {
        let located = frame(&allene(), StereoLocus::Axis(s(2))).unwrap();
        assert_eq!(located.anchors, vec![s(1), s(3)]);
        assert_eq!(located.substituents, vec![s(4), s(5), s(6), s(7)]);
    }

    #[test]
    fn an_axis_frame_is_none_off_a_cumulene_center() {
        assert!(frame(&alkene(), StereoLocus::Axis(s(1))).is_none());
    }

    #[test]
    fn candidate_loci_offer_each_site_as_a_centre_and_an_axis_and_each_bond() {
        let loci: Vec<StereoLocus> = candidate_loci(&center()).collect();
        assert_eq!(loci.len(), 5 + 5 + 4);
        assert!(loci.contains(&StereoLocus::Site(s(1))));
        assert!(loci.contains(&StereoLocus::Axis(s(1))));
        assert!(loci.contains(&StereoLocus::Bond(b(1))));
    }

    #[test]
    fn a_reflected_signal_is_the_mirror_of_the_plain_signal() {
        let mol = center().with([StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(2), s(3), s(4), s(5)],
        )
        .unwrap()]);
        let canon = refined(&mol, &|site: SiteId| site.get(), &|_: BondId| 0u8, false);
        let plain = signals(&mol, &canon, false);
        let mirror = signals(&mol, &canon, true);
        assert_eq!(mirror[&s(1)][0], plain[&s(1)][0].mirror());
    }

    #[test]
    fn refined_without_configurations_matches_the_bare_canonical_form() {
        let mol = center();
        let site_key = |_: SiteId| 0u8;
        let bond_key = |_: BondId| 0u8;
        let bare = canonicalize(&mol, site_key, bond_key);
        let refined_form = refined(&mol, &site_key, &bond_key, false);
        assert_eq!(refined_form.orbits().count(), bare.orbits().count());
    }

    #[test]
    fn refined_splits_a_symmetry_that_stereochemistry_breaks() {
        let mol = pseudo_asymmetric();
        let site_key = |site: SiteId| match site.get() {
            1 => 0u8,
            2 | 3 => 1,
            4 | 7 => 2,
            5 | 8 => 3,
            6 | 9 => 4,
            _ => unreachable!(),
        };
        let bond_key = |_: BondId| 0u8;
        let bare = canonicalize(&mol, site_key, bond_key);
        let refined_form = refined(&mol, &site_key, &bond_key, false);
        assert_eq!(class(&bare, s(2)), class(&bare, s(3)));
        assert_ne!(class(&refined_form, s(2)), class(&refined_form, s(3)));
    }

    #[test]
    fn settle_files_a_site_configuration_against_its_site() {
        let mol = center().with([StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(2), s(3), s(4), s(5)],
        )
        .unwrap()]);
        let signal = settle(&mol, &|site: SiteId| site.get(), &|_: BondId| 0u8);
        assert_eq!(signal.len(), 1);
        assert_eq!(signal[&s(1)].len(), 1);
    }

    #[test]
    fn settle_files_a_bond_configuration_against_both_ends() {
        let mol = alkene().with([StereoConfiguration::new(
            StereoLocus::Bond(b(1)),
            StereoKind::CisTrans,
            [s(3), s(4), s(5), s(6)],
        )
        .unwrap()]);
        let signal = settle(&mol, &|site: SiteId| site.get(), &|_: BondId| 0u8);
        assert_eq!(signal.len(), 2);
        assert!(signal.contains_key(&s(1)) && signal.contains_key(&s(2)));
    }

    #[test]
    fn configurations_of_a_center_are_one_per_stereoisomer() {
        let subs = [s(1), s(2), s(3), s(4)];
        let configs = configurations(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            &subs,
            |site| site.get() as usize,
        );
        assert_eq!(configs.len(), StereoKind::Tetrahedral.configuration_count());
    }

    #[test]
    fn configurations_collapse_a_repeated_substituent() {
        let subs = [s(1), s(2), s(3), s(4)];
        let configs = configurations(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            &subs,
            |site| {
                if site == s(1) || site == s(2) {
                    0
                } else {
                    site.get() as usize
                }
            },
        );
        assert_eq!(configs.len(), 1);
    }

    #[test]
    fn an_octahedral_m_a4_b2_has_a_cis_and_a_trans_configuration() {
        let subs = [s(1), s(2), s(3), s(4), s(5), s(6)];
        let configs = configurations(
            StereoLocus::Site(s(1)),
            StereoKind::Octahedral,
            &subs,
            |site| usize::from(site == s(5) || site == s(6)),
        );
        assert_eq!(configs.len(), 2);
    }

    #[test]
    fn realisable_lists_a_stereogenic_centers_configurations() {
        let configs = realisable(
            &center(),
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            &|site: SiteId| site.get(),
            &|_: BondId| 0u8,
        );
        assert_eq!(configs.len(), StereoKind::Tetrahedral.configuration_count());
    }

    #[test]
    fn realisable_collapses_a_symmetric_center_to_one() {
        let configs = realisable(
            &center(),
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            &|_: SiteId| 0u8,
            &|_: BondId| 0u8,
        );
        assert_eq!(configs.len(), 1);
    }

    #[test]
    fn realisable_is_empty_off_an_anchor() {
        let configs = realisable(
            &center(),
            StereoLocus::Bond(b(1)),
            StereoKind::Tetrahedral,
            &|site: SiteId| site.get(),
            &|_: BondId| 0u8,
        );
        assert!(configs.is_empty());
    }

    #[test]
    fn realisable_is_empty_without_a_frame() {
        let configs = realisable(
            &short_alkene(),
            StereoLocus::Bond(b(1)),
            StereoKind::CisTrans,
            &|site: SiteId| site.get(),
            &|_: BondId| 0u8,
        );
        assert!(configs.is_empty());
    }

    #[test]
    fn realisable_is_empty_when_the_substituents_miscount() {
        let configs = realisable(
            &trigonal(),
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            &|site: SiteId| site.get(),
            &|_: BondId| 0u8,
        );
        assert!(configs.is_empty());
    }

    #[test]
    fn frame_mark_numbers_anchors_from_one_and_marks_outsiders_zero() {
        let anchors = [s(4), s(7)];
        assert_eq!(frame_mark(s(4), &anchors), 1);
        assert_eq!(frame_mark(s(7), &anchors), 2);
        assert_eq!(frame_mark(s(1), &anchors), 0);
    }

    #[test]
    fn presentations_of_a_center_are_every_ordering() {
        let subs = [s(1), s(2), s(3), s(4)];
        let orderings = presentations(StereoKind::Tetrahedral, &subs);
        assert_eq!(orderings.len(), 24);
        let mut distinct = orderings.clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(distinct.len(), 24);
    }

    #[test]
    fn presentations_of_an_edge_keep_each_end_within_itself() {
        let subs = [s(1), s(2), s(3), s(4)];
        let orderings = presentations(StereoKind::CisTrans, &subs);
        assert_eq!(orderings.len(), 4);
        for ordering in orderings {
            let mut first = [ordering[0], ordering[1]];
            first.sort_unstable();
            assert_eq!(first, [s(1), s(2)]);
        }
    }
}

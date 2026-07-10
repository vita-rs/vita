mod center;
mod identity;
mod perceive;
mod validate;

pub use center::{Stereocenters, stereocenters};
pub use identity::{StereoForm, StereoRelationship, form};
pub use perceive::{StereoConfigurations, WithStereoConfigurations, perceive};
pub use validate::{StereoConsistency, consistency};

use vita_core::SiteId;

use crate::{BondOrder, HasBondOrders, StereoKind, StereoLocus};

/// The largest neighbour count of any geometry — six, the octahedron or trigonal
/// prism. Every arrangement buffers this many slots; unused ones stay zero.
const MAX_SLOTS: usize = 6;

/// A stereogenic frame located in the graph: the atoms that pin it, and the
/// substituents its geometry arranges.
///
/// A site's frame is the atom and its neighbours. An edge's or axis's is the two
/// termini of its rigid double-bond chain and the two substituents each bears,
/// walked out from the anchor — so a plain double bond and a long cumulene resolve
/// alike. The substituents of a bipartite frame are grouped by end.
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

/// The canonical class of a neighbour arrangement under a geometry's rotation
/// group.
///
/// A configuration's neighbours, ranked and relabelled to their relative order,
/// reduce under the group to the lexicographically least image — a coset
/// representative. Two configurations a ranking cannot tell apart share a token; it
/// is `Copy`, ordered, and hashed, and never leaves the module.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Token([u8; MAX_SLOTS]);

/// The idealised local geometry a [`StereoKind`] denotes, looked up once by
/// [`geometry`]: the proper-rotation group over its neighbour slots, the reflection
/// that mirrors it, the reference directions coordinates align onto, and whether its
/// slots split into two ends. The reduction, mirror, and count that consume it never
/// branch on the kind — a further geometry is new data here, not new logic.
struct Geometry {
    /// The proper-rotation group over the slots — the orderings a configuration
    /// treats as equivalent, the closure of the geometry's rotation generators.
    group: &'static [&'static [u8]],
    /// The slot permutation a reflection induces, carrying a configuration onto its
    /// mirror image: outside [`group`](Self::group) for a chiral geometry, inside it
    /// for an achiral one, so the mirror is well defined on cosets.
    reflection: &'static [u8],
    /// The reference unit directions, slot by slot, a centre's coordinates align
    /// onto; empty for a geometry read from a dihedral instead — a double bond or an
    /// allene.
    directions: &'static [[f64; 3]],
    /// Whether the slots split into two independent ends — a double bond, an allene —
    /// whose substituents cannot cross, as against one freely permuted centre.
    bipartite: bool,
}

impl Geometry {
    /// The token of a configuration whose `neighbors` are ranked by `rank`.
    ///
    /// Reduces the neighbours' *relative* order — their ranks relabelled to `0..k` by
    /// magnitude, not the ranks themselves — so symmetry-equivalent centres, whose
    /// neighbours carry different ranks in the same pattern, reduce alike.
    fn token(&self, neighbors: &[SiteId], rank: impl Fn(SiteId) -> usize) -> Token {
        Token(reduce(relative_order(neighbors, rank), self.group))
    }

    /// The token of the mirror image of the configuration `token` names.
    ///
    /// Applies the reflection and re-reduces: a chiral geometry's token flips to
    /// another class, an achiral one's stays put.
    fn mirror(&self, token: Token) -> Token {
        Token(reduce(apply(self.reflection, token.0), self.group))
    }

    /// The number of distinct configurations the geometry realises given the symmetry
    /// `classes` of its substituents, one per slot — a locus is stereogenic exactly
    /// when this exceeds one, the all-distinct case being [`StereoKind::configuration_count`].
    fn configuration_count(&self, classes: &[usize]) -> usize {
        let mut cosets: Vec<[u8; MAX_SLOTS]> = self
            .arrangements(classes)
            .into_iter()
            .map(|arrangement| reduce(arrangement, self.group))
            .collect();
        cosets.sort_unstable();
        cosets.dedup();
        cosets.len()
    }

    /// Every arrangement of the substituent `classes` the geometry admits: within
    /// each end for a bipartite geometry, over all slots for a centre. Classes are
    /// relabelled to `0..m` so the buffers compare directly.
    fn arrangements(&self, classes: &[usize]) -> Vec<[u8; MAX_SLOTS]> {
        let labels = relabel(classes);
        if self.bipartite {
            let mut result = Vec::with_capacity(4);
            for first in [[labels[0], labels[1]], [labels[1], labels[0]]] {
                for second in [[labels[2], labels[3]], [labels[3], labels[2]]] {
                    let mut slots = [0u8; MAX_SLOTS];
                    slots[..4].copy_from_slice(&[first[0], first[1], second[0], second[1]]);
                    result.push(slots);
                }
            }
            result
        } else {
            let n = labels.len();
            let mut labels = labels;
            labels.sort_unstable();
            let mut result = Vec::new();
            loop {
                let mut slots = [0u8; MAX_SLOTS];
                slots[..n].copy_from_slice(&labels);
                result.push(slots);
                if !next_permutation(&mut labels) {
                    break;
                }
            }
            result
        }
    }
}

/// The geometry each [`StereoKind`] denotes.
///
/// - `Tetrahedral` — A₄ (12), two configurations.
/// - `CisTrans`, `Allene` — `{(), (0 1)(2 3)}` (2) over the two-plus-two ends, told
///   apart by their reflection: the identity for a double bond, the end-swap for an
///   allene.
/// - `SquarePlanar` — D₄ (8), three configurations.
/// - `TrigonalBipyramidal` — D₃ (6), twenty configurations.
/// - `SquarePyramidal` — C₄ (4), thirty configurations.
/// - `Octahedral` — O (24), thirty configurations.
/// - `TrigonalPrismatic` — D₃ (6), one hundred and twenty configurations.
fn geometry(kind: StereoKind) -> &'static Geometry {
    match kind {
        StereoKind::Tetrahedral => &TETRAHEDRAL,
        StereoKind::CisTrans => &CIS_TRANS,
        StereoKind::Allene => &ALLENE,
        StereoKind::SquarePlanar => &SQUARE_PLANAR,
        StereoKind::TrigonalBipyramidal => &TRIGONAL_BIPYRAMIDAL,
        StereoKind::SquarePyramidal => &SQUARE_PYRAMIDAL,
        StereoKind::Octahedral => &OCTAHEDRAL,
        StereoKind::TrigonalPrismatic => &TRIGONAL_PRISMATIC,
    }
}

/// The two-plus-two rotation group the bipartite edge geometries share: the identity
/// and the simultaneous swap of both ends.
const PAIRED: &[&[u8]] = &[&[0, 1, 2, 3], &[1, 0, 3, 2]];

/// √3 / 2 — the in-plane offset of a trigonal vertex.
const R: f64 = 0.866_025_403_784_438_6;

static TETRAHEDRAL: Geometry = Geometry {
    group: &[
        &[0, 1, 2, 3],
        &[0, 2, 3, 1],
        &[0, 3, 1, 2],
        &[1, 0, 3, 2],
        &[1, 2, 0, 3],
        &[1, 3, 2, 0],
        &[2, 0, 1, 3],
        &[2, 1, 3, 0],
        &[2, 3, 0, 1],
        &[3, 0, 2, 1],
        &[3, 1, 0, 2],
        &[3, 2, 1, 0],
    ],
    reflection: &[1, 0, 2, 3],
    directions: &[
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
    ],
    bipartite: false,
};

static CIS_TRANS: Geometry = Geometry {
    group: PAIRED,
    reflection: &[0, 1, 2, 3],
    directions: &[],
    bipartite: true,
};

static ALLENE: Geometry = Geometry {
    group: PAIRED,
    reflection: &[1, 0, 2, 3],
    directions: &[],
    bipartite: true,
};

static SQUARE_PLANAR: Geometry = Geometry {
    group: &[
        &[0, 1, 2, 3],
        &[0, 3, 2, 1],
        &[1, 0, 3, 2],
        &[1, 2, 3, 0],
        &[2, 1, 0, 3],
        &[2, 3, 0, 1],
        &[3, 0, 1, 2],
        &[3, 2, 1, 0],
    ],
    reflection: &[0, 1, 2, 3],
    directions: &[
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ],
    bipartite: false,
};

static TRIGONAL_BIPYRAMIDAL: Geometry = Geometry {
    group: &[
        &[0, 1, 2, 3, 4],
        &[0, 1, 3, 4, 2],
        &[0, 1, 4, 2, 3],
        &[1, 0, 2, 4, 3],
        &[1, 0, 3, 2, 4],
        &[1, 0, 4, 3, 2],
    ],
    reflection: &[1, 0, 2, 3, 4],
    directions: &[
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [-0.5, R, 0.0],
        [-0.5, -R, 0.0],
    ],
    bipartite: false,
};

static SQUARE_PYRAMIDAL: Geometry = Geometry {
    group: &[
        &[0, 1, 2, 3, 4],
        &[0, 2, 3, 4, 1],
        &[0, 3, 4, 1, 2],
        &[0, 4, 1, 2, 3],
    ],
    reflection: &[0, 1, 4, 3, 2],
    directions: &[
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ],
    bipartite: false,
};

static OCTAHEDRAL: Geometry = Geometry {
    group: &[
        &[0, 1, 2, 3, 4, 5],
        &[0, 1, 3, 2, 5, 4],
        &[0, 1, 4, 5, 3, 2],
        &[0, 1, 5, 4, 2, 3],
        &[1, 0, 2, 3, 5, 4],
        &[1, 0, 3, 2, 4, 5],
        &[1, 0, 4, 5, 2, 3],
        &[1, 0, 5, 4, 3, 2],
        &[2, 3, 0, 1, 5, 4],
        &[2, 3, 1, 0, 4, 5],
        &[2, 3, 4, 5, 0, 1],
        &[2, 3, 5, 4, 1, 0],
        &[3, 2, 0, 1, 4, 5],
        &[3, 2, 1, 0, 5, 4],
        &[3, 2, 4, 5, 1, 0],
        &[3, 2, 5, 4, 0, 1],
        &[4, 5, 0, 1, 2, 3],
        &[4, 5, 1, 0, 3, 2],
        &[4, 5, 2, 3, 1, 0],
        &[4, 5, 3, 2, 0, 1],
        &[5, 4, 0, 1, 3, 2],
        &[5, 4, 1, 0, 2, 3],
        &[5, 4, 2, 3, 0, 1],
        &[5, 4, 3, 2, 1, 0],
    ],
    reflection: &[1, 0, 3, 2, 5, 4],
    directions: &[
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ],
    bipartite: false,
};

static TRIGONAL_PRISMATIC: Geometry = Geometry {
    group: &[
        &[0, 1, 2, 3, 4, 5],
        &[1, 2, 0, 4, 5, 3],
        &[2, 0, 1, 5, 3, 4],
        &[3, 5, 4, 0, 2, 1],
        &[4, 3, 5, 1, 0, 2],
        &[5, 4, 3, 2, 1, 0],
    ],
    reflection: &[0, 2, 1, 3, 5, 4],
    directions: &[
        [1.0, 0.0, 1.0],
        [-0.5, R, 1.0],
        [-0.5, -R, 1.0],
        [1.0, 0.0, -1.0],
        [-0.5, R, -1.0],
        [-0.5, -R, -1.0],
    ],
    bipartite: false,
};

/// The least image of `order` over the geometry's rotation group.
fn reduce(order: [u8; MAX_SLOTS], group: &[&[u8]]) -> [u8; MAX_SLOTS] {
    group
        .iter()
        .map(|permutation| apply(permutation, order))
        .min()
        .expect("a group contains the identity")
}

/// `order` read through `permutation`: slot `i` takes `order[permutation[i]]`.
fn apply(permutation: &[u8], order: [u8; MAX_SLOTS]) -> [u8; MAX_SLOTS] {
    let mut image = [0u8; MAX_SLOTS];
    for (slot, &source) in permutation.iter().enumerate() {
        image[slot] = order[source as usize];
    }
    image
}

/// The neighbours' ranks relabelled to their relative order `0..k`, stable under
/// ties so the mapping is a total function of the ranks alone.
fn relative_order(neighbors: &[SiteId], rank: impl Fn(SiteId) -> usize) -> [u8; MAX_SLOTS] {
    let mut ranks = [0usize; MAX_SLOTS];
    for (slot, &neighbor) in neighbors.iter().enumerate() {
        ranks[slot] = rank(neighbor);
    }
    let mut order = [0u8; MAX_SLOTS];
    for slot in 0..neighbors.len() {
        let smaller = (0..neighbors.len())
            .filter(|&other| (ranks[other], other) < (ranks[slot], slot))
            .count();
        order[slot] = smaller as u8;
    }
    order
}

/// Relabels the substituent classes to `0..m` by first appearance in sorted order.
fn relabel(classes: &[usize]) -> Vec<u8> {
    let mut distinct: Vec<usize> = classes.to_vec();
    distinct.sort_unstable();
    distinct.dedup();
    classes
        .iter()
        .map(|class| distinct.binary_search(class).expect("class is present") as u8)
        .collect()
}

/// Advances `slice` to the next lexicographic permutation, `false` at the last.
fn next_permutation(slice: &mut [u8]) -> bool {
    let n = slice.len();
    let Some(pivot) = (1..n).rev().find(|&i| slice[i - 1] < slice[i]) else {
        return false;
    };
    let successor = (pivot..n)
        .rev()
        .find(|&i| slice[i] > slice[pivot - 1])
        .expect("a successor exists past the pivot");
    slice.swap(pivot - 1, successor);
    slice[pivot..].reverse();
    true
}

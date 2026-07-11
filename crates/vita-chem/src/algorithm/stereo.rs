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

    const KINDS: [StereoKind; 8] = [
        StereoKind::Tetrahedral,
        StereoKind::CisTrans,
        StereoKind::Allene,
        StereoKind::SquarePlanar,
        StereoKind::TrigonalBipyramidal,
        StereoKind::SquarePyramidal,
        StereoKind::Octahedral,
        StereoKind::TrigonalPrismatic,
    ];

    fn identity(n: usize) -> Vec<u8> {
        (0..n as u8).collect()
    }

    fn is_permutation(perm: &[u8], n: usize) -> bool {
        let mut sorted = perm.to_vec();
        sorted.sort_unstable();
        sorted == identity(n)
    }

    fn compose(after: &[u8], before: &[u8]) -> Vec<u8> {
        before.iter().map(|&i| after[i as usize]).collect()
    }

    fn permutations(n: usize) -> Vec<Vec<u8>> {
        let mut result = Vec::new();
        permute(&mut identity(n), 0, &mut result);
        result
    }

    fn permute(slice: &mut [u8], start: usize, out: &mut Vec<Vec<u8>>) {
        if start == slice.len() {
            out.push(slice.to_vec());
            return;
        }
        for i in start..slice.len() {
            slice.swap(start, i);
            permute(slice, start + 1, out);
            slice.swap(start, i);
        }
    }

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    fn signed_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        dot(
            a,
            [
                b[1] * c[2] - b[2] * c[1],
                b[2] * c[0] - b[0] * c[2],
                b[0] * c[1] - b[1] * c[0],
            ],
        )
    }

    fn preserves_angles(directions: &[[f64; 3]], perm: &[u8]) -> bool {
        let n = directions.len();
        (0..n).all(|i| {
            (0..n).all(|j| {
                let moved = dot(directions[perm[i] as usize], directions[perm[j] as usize]);
                (moved - dot(directions[i], directions[j])).abs() < 1e-9
            })
        })
    }

    fn spanning_triple(directions: &[[f64; 3]]) -> Option<(usize, usize, usize)> {
        let n = directions.len();
        (0..n)
            .flat_map(|a| (a + 1..n).flat_map(move |b| (b + 1..n).map(move |c| (a, b, c))))
            .find(|&(a, b, c)| {
                signed_volume(directions[a], directions[b], directions[c]).abs() > 1e-9
            })
    }

    fn is_rotation(directions: &[[f64; 3]], perm: &[u8]) -> bool {
        if !preserves_angles(directions, perm) {
            return false;
        }
        match spanning_triple(directions) {
            None => true,
            Some((a, b, c)) => {
                let before = signed_volume(directions[a], directions[b], directions[c]);
                let after = signed_volume(
                    directions[perm[a] as usize],
                    directions[perm[b] as usize],
                    directions[perm[c] as usize],
                );
                (before > 0.0) == (after > 0.0)
            }
        }
    }

    fn ranked(neighbors: &[SiteId]) -> impl Fn(SiteId) -> usize + '_ {
        move |site| neighbors.iter().position(|&s| s == site).unwrap()
    }

    fn distinct_token(kind: StereoKind) -> Token {
        let neighbors: Vec<SiteId> = (1..=kind.slot_count() as u32).map(s).collect();
        geometry(kind).token(&neighbors, ranked(&neighbors))
    }

    struct Mol {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
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

    fn mol(sites: &[u32], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: sites.iter().map(|&id| s(id)).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c, _)| (s(a), s(c))).collect(),
            orders: bonds.iter().map(|&(_, _, _, order)| order).collect(),
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

    #[test]
    fn a_single_element_has_no_next_permutation() {
        assert!(!next_permutation(&mut [0u8]));
    }

    #[test]
    fn next_permutation_walks_every_permutation_in_lexicographic_order() {
        let mut slice = [0u8, 1, 2];
        let mut seen = vec![slice.to_vec()];
        while next_permutation(&mut slice) {
            seen.push(slice.to_vec());
        }
        assert_eq!(
            seen,
            vec![
                vec![0, 1, 2],
                vec![0, 2, 1],
                vec![1, 0, 2],
                vec![1, 2, 0],
                vec![2, 0, 1],
                vec![2, 1, 0],
            ],
        );
    }

    #[test]
    fn next_permutation_skips_the_repeats_of_a_multiset() {
        let mut slice = [0u8, 0, 1];
        let mut seen = vec![slice.to_vec()];
        while next_permutation(&mut slice) {
            seen.push(slice.to_vec());
        }
        assert_eq!(seen, vec![vec![0, 0, 1], vec![0, 1, 0], vec![1, 0, 0]]);
    }

    #[test]
    fn applying_the_identity_leaves_the_order_unchanged() {
        let order = [3, 1, 4, 1, 5, 9];
        assert_eq!(apply(&[0, 1, 2, 3, 4, 5], order), order);
    }

    #[test]
    fn apply_reads_each_slot_through_the_permutation() {
        assert_eq!(
            apply(&[1, 0, 2, 3, 4, 5], [3, 1, 4, 1, 5, 9]),
            [1, 3, 4, 1, 5, 9]
        );
    }

    #[test]
    fn relative_order_ranks_the_neighbors_by_magnitude() {
        let neighbors = [s(1), s(2), s(3), s(4)];
        let rank = |site| {
            if site == s(2) || site == s(4) {
                10
            } else if site == s(3) {
                20
            } else {
                30
            }
        };
        assert_eq!(relative_order(&neighbors, rank), [3, 0, 2, 1, 0, 0]);
    }

    #[test]
    fn relative_order_depends_only_on_the_relative_ranks() {
        let neighbors = [s(1), s(2), s(3), s(4)];
        let compact = |site| neighbors.iter().position(|&x| x == site).unwrap();
        let spread = |site| 100 * neighbors.iter().position(|&x| x == site).unwrap();
        assert_eq!(
            relative_order(&neighbors, compact),
            relative_order(&neighbors, spread),
        );
    }

    #[test]
    fn relabel_maps_classes_to_dense_labels_in_sorted_order() {
        assert_eq!(relabel(&[7, 3, 7, 9]), vec![1, 0, 1, 2]);
    }

    #[test]
    fn relabel_gives_identical_classes_one_label() {
        assert_eq!(relabel(&[5, 5, 5]), vec![0, 0, 0]);
    }

    #[test]
    fn every_geometry_group_permutes_its_slots() {
        for kind in KINDS {
            for &element in geometry(kind).group {
                assert!(is_permutation(element, kind.slot_count()));
            }
        }
    }

    #[test]
    fn every_geometry_group_contains_the_identity() {
        for kind in KINDS {
            let identity = identity(kind.slot_count());
            assert!(geometry(kind).group.contains(&identity.as_slice()));
        }
    }

    #[test]
    fn every_geometry_group_is_closed_under_composition() {
        for kind in KINDS {
            let group = geometry(kind).group;
            for &g in group {
                for &h in group {
                    let product = compose(g, h);
                    assert!(group.contains(&product.as_slice()));
                }
            }
        }
    }

    #[test]
    fn every_geometry_group_is_the_rotation_group_of_its_directions() {
        for kind in KINDS {
            let directions = geometry(kind).directions;
            if directions.is_empty() {
                continue;
            }
            let mut rotations: Vec<Vec<u8>> = permutations(directions.len())
                .into_iter()
                .filter(|perm| is_rotation(directions, perm))
                .collect();
            rotations.sort_unstable();
            let mut group: Vec<Vec<u8>> = geometry(kind).group.iter().map(|g| g.to_vec()).collect();
            group.sort_unstable();
            assert_eq!(group, rotations);
        }
    }

    #[test]
    fn every_geometry_admits_the_kinds_configuration_count() {
        for kind in KINDS {
            let all_distinct: Vec<usize> = (0..kind.slot_count()).collect();
            assert_eq!(
                geometry(kind).configuration_count(&all_distinct),
                kind.configuration_count(),
            );
        }
    }

    #[test]
    fn a_repeated_substituent_leaves_a_tetrahedron_with_one_configuration() {
        assert_eq!(
            geometry(StereoKind::Tetrahedral).configuration_count(&[0, 0, 1, 2]),
            1,
        );
    }

    #[test]
    fn an_octahedral_m_a4_b2_has_a_cis_and_a_trans_configuration() {
        assert_eq!(
            geometry(StereoKind::Octahedral).configuration_count(&[0, 0, 0, 0, 1, 1]),
            2,
        );
    }

    #[test]
    fn every_geometry_reflection_permutes_its_slots() {
        for kind in KINDS {
            assert!(is_permutation(geometry(kind).reflection, kind.slot_count()));
        }
    }

    #[test]
    fn a_geometry_reflection_lies_in_its_group_exactly_when_the_kind_is_achiral() {
        for kind in KINDS {
            let reflection = geometry(kind).reflection;
            let in_group = geometry(kind).group.contains(&reflection);
            assert_eq!(in_group, !kind.is_chiral());
        }
    }

    #[test]
    fn every_spatial_reflection_preserves_the_pairwise_angles() {
        for kind in KINDS {
            let directions = geometry(kind).directions;
            if !directions.is_empty() {
                assert!(preserves_angles(directions, geometry(kind).reflection));
            }
        }
    }

    #[test]
    fn only_the_edge_geometries_are_bipartite() {
        for kind in KINDS {
            let expected = matches!(kind, StereoKind::CisTrans | StereoKind::Allene);
            assert_eq!(geometry(kind).bipartite, expected);
        }
    }

    #[test]
    fn a_geometry_has_reference_directions_exactly_when_it_is_not_bipartite() {
        for kind in KINDS {
            assert_eq!(
                geometry(kind).directions.is_empty(),
                geometry(kind).bipartite
            );
        }
    }

    #[test]
    fn present_directions_number_the_slots() {
        for kind in KINDS {
            let geometry = geometry(kind);
            if !geometry.directions.is_empty() {
                assert_eq!(geometry.directions.len(), kind.slot_count());
            }
        }
    }

    #[test]
    fn a_rotation_of_the_neighbors_preserves_the_token() {
        for kind in KINDS {
            let neighbors: Vec<SiteId> = (1..=kind.slot_count() as u32).map(s).collect();
            let token = geometry(kind).token(&neighbors, ranked(&neighbors));
            for &element in geometry(kind).group {
                let rotated: Vec<SiteId> = element.iter().map(|&i| neighbors[i as usize]).collect();
                assert_eq!(geometry(kind).token(&rotated, ranked(&neighbors)), token);
            }
        }
    }

    #[test]
    fn a_reflection_returns_a_token_to_itself() {
        for kind in KINDS {
            let token = distinct_token(kind);
            assert_eq!(geometry(kind).mirror(geometry(kind).mirror(token)), token);
        }
    }

    #[test]
    fn an_achiral_geometry_mirror_fixes_every_token() {
        for kind in [StereoKind::CisTrans, StereoKind::SquarePlanar] {
            let token = distinct_token(kind);
            assert_eq!(geometry(kind).mirror(token), token);
        }
    }

    #[test]
    fn a_chiral_geometry_mirror_alters_a_distinct_token() {
        for kind in KINDS {
            if kind.is_chiral() {
                let token = distinct_token(kind);
                assert_ne!(geometry(kind).mirror(token), token);
            }
        }
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
    fn a_double_bond_frame_is_none_when_a_terminus_lacks_two_substituents() {
        assert!(frame(&short_alkene(), StereoLocus::Bond(b(1))).is_none());
    }

    #[test]
    fn a_cumulene_bond_frame_walks_out_to_the_chain_termini() {
        let located = frame(&butatriene(), StereoLocus::Bond(b(2))).unwrap();
        assert_eq!(located.anchors, vec![s(1), s(4)]);
        assert_eq!(located.substituents, vec![s(5), s(6), s(7), s(8)]);
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
}

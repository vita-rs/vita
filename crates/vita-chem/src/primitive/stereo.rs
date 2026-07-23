use vita_core::SiteId;

use crate::BondId;

/// The largest neighbor count of any geometry — six, the octahedron or trigonal
/// prism. The fixed-size slot buffers hold this many; unused entries stay zero.
const MAX_SLOTS: usize = 6;

/// Where a stereogenic unit is anchored in a molecule.
///
/// Stereochemistry is carried by the arrangement of a site's substituents (a
/// coordination center), of the substituents across a bond (a double bond's two
/// ends), or across an axis (an allene's two termini); `StereoLocus` names that
/// anchor. It labels a unit for reporting and reconciliation and carries no meaning
/// of its own — the arrangement that fixes which orderings are equivalent is the
/// unit's [`StereoKind`].
///
/// A site locus orders before a bond, a bond before an axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoLocus {
    /// A site whose substituents' arrangement is stereogenic.
    Site(SiteId),
    /// A bond whose ends' substituents' arrangement is stereogenic.
    Bond(BondId),
    /// An axis, named by its central site, whose termini's substituents'
    /// arrangement is stereogenic.
    Axis(SiteId),
}

impl StereoLocus {
    /// Whether this locus is the anchor `kind` lives on: a site for a coordination
    /// center, a bond for a double bond, an axis for an allene.
    #[inline]
    pub const fn anchors(self, kind: StereoKind) -> bool {
        matches!(
            (self, kind),
            (
                Self::Site(_),
                StereoKind::Tetrahedral
                    | StereoKind::SquarePlanar
                    | StereoKind::TrigonalBipyramidal
                    | StereoKind::SquarePyramidal
                    | StereoKind::Octahedral
                    | StereoKind::TrigonalPrismatic,
            ) | (Self::Bond(_), StereoKind::CisTrans)
                | (Self::Axis(_), StereoKind::Allene)
        )
    }
}

/// The kind of a stereogenic unit: the idealised local geometry whose rotation
/// group fixes which of its neighbor orderings are equivalent.
///
/// A kind is a pure data key. It selects the permutation group under which a
/// configuration's neighbor ordering reduces, and the geometric reference against
/// which coordinates are perceived.
///
/// The kinds order by neighbor count, then configuration count, then locus —
/// a center before a bond before an axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoKind {
    /// A tetrahedral center (4 neighbors, 2 configurations).
    Tetrahedral,
    /// A double bond (4 neighbors, 2 configurations).
    CisTrans,
    /// An allene axis (4 neighbors, 2 configurations).
    Allene,
    /// A square-planar center (4 neighbors, 3 configurations).
    SquarePlanar,
    /// A trigonal-bipyramidal center (5 neighbors, 20 configurations).
    TrigonalBipyramidal,
    /// A square-pyramidal center (5 neighbors, 30 configurations).
    SquarePyramidal,
    /// An octahedral center (6 neighbors, 30 configurations).
    Octahedral,
    /// A trigonal-prismatic center (6 neighbors, 120 configurations).
    TrigonalPrismatic,
}

impl StereoKind {
    /// Returns the number of neighbor slots the geometry arranges.
    #[inline]
    pub const fn slot_count(self) -> usize {
        geometry(self).slot_count()
    }

    /// Returns the number of distinct configurations the geometry admits — the
    /// stereoisomers its slots realize when every substituent differs.
    #[inline]
    pub const fn configuration_count(self) -> usize {
        geometry(self).configuration_count()
    }

    /// Returns whether the geometry is chiral — whether a configuration and its
    /// mirror image are distinct, no rotation carrying one onto the other.
    #[inline]
    pub const fn is_chiral(self) -> bool {
        geometry(self).is_chiral()
    }

    /// The reference unit directions its coordinates align onto, slot by slot — empty
    /// for a geometry read from a dihedral instead, a double bond or an allene.
    #[inline]
    pub(crate) const fn directions(self) -> &'static [[f64; 3]] {
        geometry(self).directions
    }

    /// The slot permutation a reflection induces, carrying a configuration onto its
    /// mirror image.
    #[inline]
    pub(crate) const fn reflection(self) -> &'static [u8] {
        geometry(self).reflection
    }

    /// The number of independent ends the slots split into — one freely permuted
    /// center, or the two ends of a rigid edge whose substituents cannot cross.
    #[inline]
    pub(crate) const fn ends(self) -> usize {
        geometry(self).ends
    }
}

/// A stereodescriptor: which of a geometry's configurations a unit realizes, read
/// against a ranking of its substituents.
///
/// A configuration is a coset of neighbor orderings under the geometry's rotation
/// group — the orderings a rotation cannot tell apart. Its descriptor is that
/// coset's canonical representative, computed from the neighbors' *ranks* rather
/// than the neighbors themselves, so symmetry-equivalent units reduce alike and
/// two units of the same kind are comparable. It is opaque and relative: it
/// distinguishes a configuration from the geometry's others and from its
/// [`mirror`](Self::mirror), but names none of them — no descriptor is "R" or "Λ".
/// Naming a configuration takes a convention the library does not invent; telling
/// them apart takes only this.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoDescriptor {
    kind: StereoKind,
    coset: [u8; MAX_SLOTS],
}

impl StereoDescriptor {
    /// The kind whose configurations this descriptor selects among.
    #[inline]
    pub const fn kind(&self) -> StereoKind {
        self.kind
    }

    /// The descriptor of the mirror-image configuration.
    ///
    /// A chiral kind's descriptor flips to a distinct one; an achiral kind's is its
    /// own mirror, so the two coincide.
    #[inline]
    pub fn mirror(&self) -> StereoDescriptor {
        let geometry = geometry(self.kind);
        StereoDescriptor {
            kind: self.kind,
            coset: reduce(apply(geometry.reflection, self.coset), geometry.group),
        }
    }
}

/// A stereo configuration at one [`StereoLocus`].
///
/// The *order* of the neighbors fixes the configuration: it is a reference
/// arrangement, and every reordering the unit's [`StereoKind`] treats as equivalent
/// denotes the same one. The library reads this order and never invents it — a
/// source (a wedge, coordinates, a SMILES `@`) already committed to it, and must
/// commit to the convention here. Equal configurations compare equal however their
/// neighbors were listed: a rotation of the reference order is the same
/// configuration, and [`new`](Self::new) stores every configuration in one canonical
/// order so `==`, `Ord`, and `Hash` see through the choice.
///
/// Each kind fills its slots from the neighbor list; a chiral kind's reference order
/// is additionally the one whose positions span a positive signed volume, its mirror
/// the order the kind's reflection gives. An achiral kind — a double bond, a square
/// plane — equates an arrangement with its mirror, so its order fixes no handedness.
///
/// - [`Tetrahedral`](StereoKind::Tetrahedral) — the four substituents
///   `[n0, n1, n2, n3]`; reference `(n1 − n0) · ((n2 − n0) × (n3 − n0))`.
/// - [`CisTrans`](StereoKind::CisTrans) — one end's two, then the other's,
///   `[e1a, e1b, e2a, e2b]`, `e1a` on `e2a`'s side.
/// - [`Allene`](StereoKind::Allene) — one terminus's two, then the other's,
///   `[t1a, t1b, t2a, t2b]`; reference the twist across the axis, its mirror
///   swapping `t1a` and `t1b`.
/// - [`SquarePlanar`](StereoKind::SquarePlanar) — the four in cyclic order
///   `[n0, n1, n2, n3]`, `n0`/`n2` and `n1`/`n3` trans.
/// - [`TrigonalBipyramidal`](StereoKind::TrigonalBipyramidal) — the two axial, then
///   the three equatorial, `[a0, a1, e0, e1, e2]`; reference
///   `(a0 − a1) · ((e1 − e0) × (e2 − e0))`.
/// - [`SquarePyramidal`](StereoKind::SquarePyramidal) — the apical, then the four
///   basal in cyclic order, `[p, b0, b1, b2, b3]`; reference
///   `(p − b0) · ((b1 − b0) × (b2 − b0))`.
/// - [`Octahedral`](StereoKind::Octahedral) — three trans pairs
///   `[n0, n1, n2, n3, n4, n5]`, `n2i` opposite `n2i+1`; reference
///   `(n0 − n1) · ((n2 − n3) × (n4 − n5))`.
/// - [`TrigonalPrismatic`](StereoKind::TrigonalPrismatic) — one triangular face,
///   then the other, `[t0, t1, t2, b0, b1, b2]`, `bi` eclipsing `ti`; reference
///   `(t0 − b0) · ((t1 − t0) × (t2 − t0))`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoConfiguration {
    locus: StereoLocus,
    kind: StereoKind,
    neighbors: Vec<SiteId>,
}

impl StereoConfiguration {
    /// Builds a configuration at `locus` of `kind` from its `neighbors`, given in the
    /// order the kind's contract prescribes.
    ///
    /// Returns `None` unless `locus` [anchors](StereoLocus::anchors) `kind` and its
    /// neighbor count is the kind's [`StereoKind::slot_count`].
    #[inline]
    pub fn new(
        locus: StereoLocus,
        kind: StereoKind,
        neighbors: impl IntoIterator<Item = SiteId>,
    ) -> Option<Self> {
        let neighbors: Vec<SiteId> = neighbors.into_iter().collect();
        (locus.anchors(kind) && neighbors.len() == kind.slot_count()).then(|| Self {
            locus,
            kind,
            neighbors: normalized(kind, neighbors),
        })
    }

    /// The anchor of the stereogenic unit.
    #[inline]
    pub const fn locus(&self) -> StereoLocus {
        self.locus
    }

    /// The idealised geometry fixing which orderings are equivalent.
    #[inline]
    pub const fn kind(&self) -> StereoKind {
        self.kind
    }

    /// The neighbors in canonical reference order (see the type's contract).
    #[inline]
    pub fn neighbors(&self) -> &[SiteId] {
        &self.neighbors
    }

    /// This configuration's [`StereoDescriptor`] under the substituent ranking `rank`.
    ///
    /// Two configurations of the same kind share a descriptor exactly when `rank`
    /// cannot tell their arrangements apart; a chiral configuration's differs from
    /// its mirror's. Pass the symmetry classes to identify a stereoisomer, a total
    /// priority order to name a configuration.
    #[inline]
    pub fn descriptor(&self, rank: impl Fn(SiteId) -> usize) -> StereoDescriptor {
        let geometry = geometry(self.kind);
        StereoDescriptor {
            kind: self.kind,
            coset: reduce(relative_order(&self.neighbors, rank), geometry.group),
        }
    }
}

/// The idealised local geometry a [`StereoKind`] denotes, looked up once by
/// [`geometry`]: the proper-rotation group over its neighbor slots, the reflection
/// that mirrors it, the reference directions coordinates align onto, and how many
/// independent ends its slots split into. Every fact about a kind is a function of
/// this one datum — a further geometry is new data here, not new logic.
struct Geometry {
    /// The proper-rotation group over the slots — the orderings a configuration
    /// treats as equivalent, the closure of the geometry's rotation generators.
    group: &'static [&'static [u8]],
    /// The slot permutation a reflection induces, carrying a configuration onto its
    /// mirror image: outside [`group`](Self::group) for a chiral geometry, inside it
    /// for an achiral one, so the mirror is well defined on cosets.
    reflection: &'static [u8],
    /// The reference unit directions, slot by slot, a center's coordinates align
    /// onto; empty for a geometry read from a dihedral instead — a double bond or an
    /// allene.
    directions: &'static [[f64; 3]],
    /// The number of independent ends the slots split into: one freely permuted
    /// center, or two ends of a rigid frame — a double bond, an allene — whose
    /// substituents cannot cross between them.
    ends: usize,
}

impl Geometry {
    /// The number of neighbor slots — the length of any group element.
    const fn slot_count(&self) -> usize {
        self.reflection.len()
    }

    /// Whether the geometry is chiral: its reflection is not itself a rotation.
    const fn is_chiral(&self) -> bool {
        let mut index = 0;
        while index < self.group.len() {
            if slices_equal(self.group[index], self.reflection) {
                return false;
            }
            index += 1;
        }
        true
    }

    /// The distinct configurations the all-different case admits: the orderings that
    /// respect the ends, quotiented by the rotation group. The presentations number
    /// `ends! · (slots / ends)!^ends` — the ends ordered, then each filled — and the
    /// group acts on them freely when every substituent differs, so the quotient is
    /// their count over the group's order.
    const fn configuration_count(&self) -> usize {
        let per_end = self.slot_count() / self.ends;
        factorial(self.ends) * factorial(per_end).pow(self.ends as u32) / self.group.len()
    }
}

/// The geometry each [`StereoKind`] denotes.
///
/// - `Tetrahedral` — A₄ (12), two configurations.
/// - `CisTrans`, `Allene` — the rigid-frame rotations [`EDGE`] (4) over the
///   two-plus-two ends, told apart by their reflection: an in-plane mirror for the
///   planar double bond, a diagonal one for the perpendicular allene.
/// - `SquarePlanar` — D₄ (8), three configurations.
/// - `TrigonalBipyramidal` — D₃ (6), twenty configurations.
/// - `SquarePyramidal` — C₄ (4), thirty configurations.
/// - `Octahedral` — O (24), thirty configurations.
/// - `TrigonalPrismatic` — D₃ (6), one hundred and twenty configurations.
const fn geometry(kind: StereoKind) -> &'static Geometry {
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

/// The proper-rotation group both two-ended edge geometries share — a double bond,
/// an allene: the identity and the three half-turns of a rigid four-substituent
/// frame (about its axis, and about each of the two axes across it, one of which
/// swaps the ends). Chirality then turns on the reflection alone.
const EDGE: &[&[u8]] = &[&[0, 1, 2, 3], &[1, 0, 3, 2], &[2, 3, 0, 1], &[3, 2, 1, 0]];

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
    ends: 1,
};

static CIS_TRANS: Geometry = Geometry {
    group: EDGE,
    reflection: &[0, 1, 2, 3],
    directions: &[],
    ends: 2,
};

static ALLENE: Geometry = Geometry {
    group: EDGE,
    reflection: &[1, 0, 2, 3],
    directions: &[],
    ends: 2,
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
    ends: 1,
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
    ends: 1,
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
    ends: 1,
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
    ends: 1,
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
    ends: 1,
};

/// The lexicographically least neighbor ordering in a configuration's rotation
/// orbit — the canonical presentation, so equal configurations store equal
/// neighbors and `==`, `Ord`, and `Hash` follow from the stored order. Only proper
/// rotations act, leaving intact the handedness a chiral kind fixes.
fn normalized(kind: StereoKind, neighbors: Vec<SiteId>) -> Vec<SiteId> {
    geometry(kind)
        .group
        .iter()
        .map(|permutation| {
            permutation
                .iter()
                .map(|&slot| neighbors[slot as usize])
                .collect()
        })
        .min()
        .expect("a group contains the identity")
}

/// The least image of `order` over a rotation `group`.
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

/// The neighbors' ranks relabeled to their relative order — each slot the count
/// of neighbors of strictly lower rank.
///
/// Ties collapse: equal ranks share a slot value, so the result is a function of
/// the ranks alone, blind to the order the neighbors were given in. Ranks a
/// symmetry cannot yet tell apart therefore reduce alike, which is what lets a
/// descriptor refine a coloring — where the ranks are the current symmetry classes,
/// not a total labeling — without depending on that labeling.
fn relative_order(neighbors: &[SiteId], rank: impl Fn(SiteId) -> usize) -> [u8; MAX_SLOTS] {
    let mut ranks = [0usize; MAX_SLOTS];
    for (slot, &neighbor) in neighbors.iter().enumerate() {
        ranks[slot] = rank(neighbor);
    }
    let mut order = [0u8; MAX_SLOTS];
    for slot in 0..neighbors.len() {
        order[slot] = (0..neighbors.len())
            .filter(|&other| ranks[other] < ranks[slot])
            .count() as u8;
    }
    order
}

/// `n!`.
const fn factorial(n: usize) -> usize {
    let mut product = 1;
    let mut factor = 2;
    while factor <= n {
        product *= factor;
        factor += 1;
    }
    product
}

/// Whether two slices hold the same bytes — a `const` slice equality, for want of
/// one in the standard library.
const fn slices_equal(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut index = 0;
    while index < a.len() {
        if a[index] != b[index] {
            return false;
        }
        index += 1;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn neighbors(count: usize) -> Vec<SiteId> {
        (1..=count as u32).map(s).collect()
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

    fn rotation_group(directions: &[[f64; 3]]) -> Vec<Vec<u8>> {
        let mut rotations: Vec<Vec<u8>> = permutations(directions.len())
            .into_iter()
            .filter(|perm| is_rotation(directions, perm))
            .collect();
        rotations.sort_unstable();
        rotations
    }

    fn sorted_group(group: &[&[u8]]) -> Vec<Vec<u8>> {
        let mut group: Vec<Vec<u8>> = group.iter().map(|element| element.to_vec()).collect();
        group.sort_unstable();
        group
    }

    #[test]
    fn a_site_anchors_every_coordination_center() {
        for kind in [
            StereoKind::Tetrahedral,
            StereoKind::SquarePlanar,
            StereoKind::TrigonalBipyramidal,
            StereoKind::SquarePyramidal,
            StereoKind::Octahedral,
            StereoKind::TrigonalPrismatic,
        ] {
            assert!(StereoLocus::Site(s(1)).anchors(kind), "{kind:?}");
        }
    }

    #[test]
    fn a_bond_anchors_a_double_bond() {
        assert!(StereoLocus::Bond(b(1)).anchors(StereoKind::CisTrans));
    }

    #[test]
    fn an_axis_anchors_an_allene() {
        assert!(StereoLocus::Axis(s(1)).anchors(StereoKind::Allene));
    }

    #[test]
    fn only_a_site_anchors_a_coordination_center() {
        assert!(!StereoLocus::Bond(b(1)).anchors(StereoKind::Tetrahedral));
        assert!(!StereoLocus::Axis(s(1)).anchors(StereoKind::Tetrahedral));
    }

    #[test]
    fn only_a_bond_anchors_a_double_bond() {
        assert!(!StereoLocus::Site(s(1)).anchors(StereoKind::CisTrans));
        assert!(!StereoLocus::Axis(s(1)).anchors(StereoKind::CisTrans));
    }

    #[test]
    fn only_an_axis_anchors_an_allene() {
        assert!(!StereoLocus::Site(s(1)).anchors(StereoKind::Allene));
        assert!(!StereoLocus::Bond(b(1)).anchors(StereoKind::Allene));
    }

    #[test]
    fn loci_order_by_anchor_then_by_identifier() {
        assert!(StereoLocus::Site(s(9)) < StereoLocus::Bond(b(1)));
        assert!(StereoLocus::Bond(b(9)) < StereoLocus::Axis(s(1)));
        assert!(StereoLocus::Site(s(1)) < StereoLocus::Site(s(2)));
    }

    #[test]
    fn slot_count_is_the_neighbor_count_of_the_geometry() {
        assert_eq!(StereoKind::Tetrahedral.slot_count(), 4);
        assert_eq!(StereoKind::CisTrans.slot_count(), 4);
        assert_eq!(StereoKind::Allene.slot_count(), 4);
        assert_eq!(StereoKind::SquarePlanar.slot_count(), 4);
        assert_eq!(StereoKind::TrigonalBipyramidal.slot_count(), 5);
        assert_eq!(StereoKind::SquarePyramidal.slot_count(), 5);
        assert_eq!(StereoKind::Octahedral.slot_count(), 6);
        assert_eq!(StereoKind::TrigonalPrismatic.slot_count(), 6);
    }

    #[test]
    fn configuration_count_is_the_number_of_distinct_stereoisomers() {
        assert_eq!(StereoKind::Tetrahedral.configuration_count(), 2);
        assert_eq!(StereoKind::CisTrans.configuration_count(), 2);
        assert_eq!(StereoKind::Allene.configuration_count(), 2);
        assert_eq!(StereoKind::SquarePlanar.configuration_count(), 3);
        assert_eq!(StereoKind::TrigonalBipyramidal.configuration_count(), 20);
        assert_eq!(StereoKind::SquarePyramidal.configuration_count(), 30);
        assert_eq!(StereoKind::Octahedral.configuration_count(), 30);
        assert_eq!(StereoKind::TrigonalPrismatic.configuration_count(), 120);
    }

    #[test]
    fn chiral_geometries_are_chiral() {
        assert!(StereoKind::Tetrahedral.is_chiral());
        assert!(StereoKind::Allene.is_chiral());
        assert!(StereoKind::TrigonalBipyramidal.is_chiral());
        assert!(StereoKind::SquarePyramidal.is_chiral());
        assert!(StereoKind::Octahedral.is_chiral());
        assert!(StereoKind::TrigonalPrismatic.is_chiral());
    }

    #[test]
    fn cis_trans_and_square_planar_are_not_chiral() {
        assert!(!StereoKind::CisTrans.is_chiral());
        assert!(!StereoKind::SquarePlanar.is_chiral());
    }

    #[test]
    fn kinds_order_by_neighbor_count_then_configuration_count_then_locus() {
        assert!(StereoKind::SquarePlanar < StereoKind::TrigonalBipyramidal);
        assert!(StereoKind::SquarePyramidal < StereoKind::Octahedral);
        assert!(StereoKind::TrigonalBipyramidal < StereoKind::SquarePyramidal);
        assert!(StereoKind::Octahedral < StereoKind::TrigonalPrismatic);
        assert!(StereoKind::Tetrahedral < StereoKind::CisTrans);
        assert!(StereoKind::CisTrans < StereoKind::Allene);
    }

    #[test]
    fn new_rejects_a_kind_off_its_anchor() {
        assert!(
            StereoConfiguration::new(
                StereoLocus::Bond(b(1)),
                StereoKind::Tetrahedral,
                neighbors(4)
            )
            .is_none()
        );
    }

    #[test]
    fn new_rejects_an_empty_neighbor_list() {
        assert!(
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::Tetrahedral,
                neighbors(0)
            )
            .is_none()
        );
    }

    #[test]
    fn new_rejects_more_neighbors_than_the_slot_count() {
        assert!(
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::Tetrahedral,
                neighbors(5)
            )
            .is_none()
        );
    }

    #[test]
    fn locus_returns_the_anchor() {
        let config = StereoConfiguration::new(
            StereoLocus::Site(s(7)),
            StereoKind::Tetrahedral,
            neighbors(4),
        )
        .unwrap();
        assert_eq!(config.locus(), StereoLocus::Site(s(7)));
    }

    #[test]
    fn kind_returns_the_geometry() {
        let config = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            neighbors(4),
        )
        .unwrap();
        assert_eq!(config.kind(), StereoKind::Tetrahedral);
    }

    #[test]
    fn neighbors_are_the_configuration_substituents() {
        let config = StereoConfiguration::new(
            StereoLocus::Site(s(9)),
            StereoKind::Tetrahedral,
            [s(4), s(2), s(3), s(1)],
        )
        .unwrap();
        let mut present = config.neighbors().to_vec();
        present.sort_unstable();
        assert_eq!(present, vec![s(1), s(2), s(3), s(4)]);
    }

    #[test]
    fn a_rotation_of_the_neighbors_yields_an_equal_configuration() {
        let reference = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(1), s(2), s(3), s(4)],
        )
        .unwrap();
        let rotated = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(2), s(1), s(4), s(3)],
        )
        .unwrap();
        assert_eq!(reference, rotated);
    }

    #[test]
    fn a_reflection_of_the_neighbors_yields_a_different_configuration() {
        let reference = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(1), s(2), s(3), s(4)],
        )
        .unwrap();
        let mirror = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(2), s(1), s(3), s(4)],
        )
        .unwrap();
        assert_ne!(reference, mirror);
    }

    #[test]
    fn configurations_differ_when_their_locus_or_kind_or_neighbors_differ() {
        let base = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            neighbors(4),
        )
        .unwrap();
        assert_ne!(
            base,
            StereoConfiguration::new(
                StereoLocus::Site(s(2)),
                StereoKind::Tetrahedral,
                neighbors(4),
            )
            .unwrap(),
        );
        assert_ne!(
            base,
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::SquarePlanar,
                neighbors(4),
            )
            .unwrap(),
        );
        assert_ne!(
            base,
            StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                StereoKind::Tetrahedral,
                [s(1), s(2), s(3), s(5)],
            )
            .unwrap(),
        );
    }

    #[test]
    fn a_descriptor_carries_its_configurations_kind() {
        let config = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::SquarePlanar,
            neighbors(4),
        )
        .unwrap();
        assert_eq!(
            config.descriptor(|site| site.get() as usize).kind(),
            StereoKind::SquarePlanar,
        );
    }

    #[test]
    fn a_chiral_descriptor_differs_from_its_mirror() {
        let config = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            neighbors(4),
        )
        .unwrap();
        let descriptor = config.descriptor(|site| site.get() as usize);
        assert_ne!(descriptor, descriptor.mirror());
    }

    #[test]
    fn an_achiral_descriptor_is_its_own_mirror() {
        let config =
            StereoConfiguration::new(StereoLocus::Bond(b(1)), StereoKind::CisTrans, neighbors(4))
                .unwrap();
        let descriptor = config.descriptor(|site| site.get() as usize);
        assert_eq!(descriptor, descriptor.mirror());
    }

    #[test]
    fn a_descriptor_ignores_a_rotation_of_the_neighbors() {
        let rank = |site: SiteId| site.get() as usize;
        let reference = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(1), s(2), s(3), s(4)],
        )
        .unwrap();
        let rotated = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(3), s(4), s(1), s(2)],
        )
        .unwrap();
        assert_eq!(reference.descriptor(rank), rotated.descriptor(rank));
    }

    #[test]
    fn a_descriptor_ignores_the_order_of_equally_ranked_neighbors() {
        let rank = |site: SiteId| match site.get() {
            2 | 3 => 0,
            other => other as usize,
        };
        let reference = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(1), s(2), s(3), s(4)],
        )
        .unwrap();
        let swapped = StereoConfiguration::new(
            StereoLocus::Site(s(1)),
            StereoKind::Tetrahedral,
            [s(1), s(3), s(2), s(4)],
        )
        .unwrap();
        assert_eq!(reference.descriptor(rank), swapped.descriptor(rank));
    }

    #[test]
    fn every_geometry_group_permutes_its_slots() {
        for kind in KINDS {
            for &element in geometry(kind).group {
                assert!(is_permutation(element, kind.slot_count()), "{kind:?}");
            }
        }
    }

    #[test]
    fn every_geometry_group_contains_the_identity() {
        for kind in KINDS {
            let identity = identity(kind.slot_count());
            assert!(
                geometry(kind).group.contains(&identity.as_slice()),
                "{kind:?}"
            );
        }
    }

    #[test]
    fn every_geometry_group_is_closed_under_composition() {
        for kind in KINDS {
            let group = geometry(kind).group;
            for &g in group {
                for &h in group {
                    let product = compose(g, h);
                    assert!(group.contains(&product.as_slice()), "{kind:?}");
                }
            }
        }
    }

    #[test]
    fn a_centers_group_is_the_rotation_group_of_its_reference_directions() {
        for kind in KINDS {
            let directions = geometry(kind).directions;
            if directions.is_empty() {
                continue;
            }
            assert_eq!(
                sorted_group(geometry(kind).group),
                rotation_group(directions),
                "{kind:?}",
            );
        }
    }

    #[test]
    fn the_edge_group_is_the_rotation_group_of_an_allene_frame() {
        let allene = [
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
            [1.0, 0.0, 1.0],
            [1.0, 0.0, -1.0],
        ];
        assert_eq!(sorted_group(EDGE), rotation_group(&allene));
    }

    #[test]
    fn every_geometry_reflection_permutes_its_slots() {
        for kind in KINDS {
            assert!(
                is_permutation(geometry(kind).reflection, kind.slot_count()),
                "{kind:?}"
            );
        }
    }

    #[test]
    fn a_centers_reflection_preserves_the_pairwise_angles() {
        for kind in KINDS {
            let directions = geometry(kind).directions;
            if !directions.is_empty() {
                assert!(
                    preserves_angles(directions, geometry(kind).reflection),
                    "{kind:?}"
                );
            }
        }
    }

    #[test]
    fn a_two_ended_geometry_has_no_reference_directions() {
        for kind in KINDS {
            assert_eq!(
                geometry(kind).directions.is_empty(),
                geometry(kind).ends == 2,
                "{kind:?}"
            );
        }
    }

    #[test]
    fn present_directions_number_the_slots() {
        for kind in KINDS {
            let directions = geometry(kind).directions;
            if !directions.is_empty() {
                assert_eq!(directions.len(), kind.slot_count(), "{kind:?}");
            }
        }
    }

    #[test]
    fn only_the_edge_geometries_have_two_ends() {
        for kind in KINDS {
            let two_ended = matches!(kind, StereoKind::CisTrans | StereoKind::Allene);
            assert_eq!(geometry(kind).ends == 2, two_ended, "{kind:?}");
        }
    }
}

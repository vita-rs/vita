/// The idealized arrangement of a site's substituents about it.
///
/// A geometry is fixed by where the substituents sit, not by what put them
/// there: lone pairs choose the arrangement but are no part of it, so one
/// geometry serves however many the site carries. Distortion is the rule; a
/// real structure is described by the arrangement it lies nearest.
///
/// The members are the polyhedral symbols of *Nomenclature of Inorganic
/// Chemistry* (IUPAC Recommendations 2005), Table IR-9.2, up to coordination
/// number six: a bound of this crate's, not the standard's, whose table runs to
/// nine — enumerating a geometry's slot orderings costs `n!`, which six still
/// affords. The set grows as the standard is extended.
///
/// The geometries order by coordination number, then by configuration count,
/// and where those agree as Table IR-9.2 lists them.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CoordinationGeometry {
    /// Two substituents opposed — *L-2* (1 configuration).
    Linear,
    /// Two substituents subtending less than a straight angle — *A-2* (1
    /// configuration).
    Angular,
    /// Three substituents about a coplanar site — *TP-3* (1 configuration).
    TrigonalPlanar,
    /// Three substituents about an apical site — *TPY-3* (2 configurations).
    TrigonalPyramidal,
    /// Two substituents opposed and one across them — *TS-3* (3 configurations).
    TShaped,
    /// Four substituents at the vertices of a tetrahedron — *T-4* (2
    /// configurations).
    Tetrahedral,
    /// Four substituents squared about a coplanar site — *SP-4* (3 configurations).
    SquarePlanar,
    /// Four substituents squared about a site off their plane — *SPY-4* (6
    /// configurations).
    PyramidalizedSquare,
    /// Two substituents opposed and two across them — *SS-4* (12 configurations).
    Seesaw,
    /// Five substituents at the vertices of a trigonal bipyramid — *TBPY-5* (20
    /// configurations).
    TrigonalBipyramidal,
    /// Five substituents at the vertices of a square pyramid — *SPY-5* (30
    /// configurations).
    SquarePyramidal,
    /// Six substituents at the vertices of an octahedron — *OC-6* (30
    /// configurations).
    Octahedral,
    /// Six substituents at the vertices of a trigonal prism — *TPR-6* (120
    /// configurations).
    TrigonalPrismatic,
}

impl CoordinationGeometry {
    /// Returns the number of slots the geometry arranges — its coordination
    /// number.
    #[inline]
    pub const fn slot_count(self) -> usize {
        idealization(self).slot_count()
    }

    /// Returns the number of distinct configurations the geometry admits — the
    /// stereoisomers its slots realize when every substituent differs.
    ///
    /// One means the geometry equates every ordering of its substituents, and
    /// so carries no stereochemistry however it is substituted.
    #[inline]
    pub const fn configuration_count(self) -> usize {
        idealization(self).configuration_count()
    }

    /// Returns whether the geometry is chiral — whether a configuration and its
    /// mirror image are distinct, no rotation carrying one onto the other.
    #[inline]
    pub const fn is_chiral(self) -> bool {
        idealization(self).is_chiral()
    }
}

/// The idealization a [`CoordinationGeometry`] denotes, looked up once by
/// [`idealization`]: the proper-rotation group over its slots, and the
/// reflection that mirrors it. Every fact about a geometry is a function of this
/// one datum — a further geometry is new data here, not new logic.
struct Idealization {
    /// The proper-rotation group over the slots — the orderings a
    /// configuration treats as equivalent, the closure of the geometry's
    /// rotation generators.
    group: &'static [&'static [u8]],
    /// The slot permutation a reflection induces, carrying a configuration onto
    /// its mirror image: outside [`group`](Self::group) for a chiral geometry,
    /// inside it for an achiral one, so the mirror is well defined on cosets.
    reflection: &'static [u8],
}

impl Idealization {
    /// The number of slots — the length of any group element.
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

    /// The distinct configurations the all-different case admits: the orderings
    /// of the slots quotiented by the rotation group, which acts on them
    /// freely when every substituent differs.
    const fn configuration_count(&self) -> usize {
        factorial(self.slot_count()) / self.group.len()
    }
}

/// The idealization each [`CoordinationGeometry`] denotes.
///
/// - `Linear`, `Angular` — S₂ (2), one configuration.
/// - `TrigonalPlanar` — D₃ (6), one configuration.
/// - `TrigonalPyramidal` — C₃ (3), two configurations.
/// - `TShaped` — C₂ (2), three configurations.
/// - `Tetrahedral` — T (12), two configurations.
/// - `SquarePlanar` — D₄ (8), three configurations.
/// - `PyramidalizedSquare` — C₄ (4), six configurations.
/// - `Seesaw` — C₂ (2), twelve configurations.
/// - `TrigonalBipyramidal` — D₃ (6), twenty configurations.
/// - `SquarePyramidal` — C₄ (4), thirty configurations.
/// - `Octahedral` — O (24), thirty configurations.
/// - `TrigonalPrismatic` — D₃ (6), one hundred and twenty configurations.
const fn idealization(geometry: CoordinationGeometry) -> &'static Idealization {
    match geometry {
        CoordinationGeometry::Linear => &LINEAR,
        CoordinationGeometry::Angular => &ANGULAR,
        CoordinationGeometry::TrigonalPlanar => &TRIGONAL_PLANAR,
        CoordinationGeometry::TrigonalPyramidal => &TRIGONAL_PYRAMIDAL,
        CoordinationGeometry::TShaped => &T_SHAPED,
        CoordinationGeometry::Tetrahedral => &TETRAHEDRAL,
        CoordinationGeometry::SquarePlanar => &SQUARE_PLANAR,
        CoordinationGeometry::PyramidalizedSquare => &PYRAMIDALIZED_SQUARE,
        CoordinationGeometry::Seesaw => &SEESAW,
        CoordinationGeometry::TrigonalBipyramidal => &TRIGONAL_BIPYRAMIDAL,
        CoordinationGeometry::SquarePyramidal => &SQUARE_PYRAMIDAL,
        CoordinationGeometry::Octahedral => &OCTAHEDRAL,
        CoordinationGeometry::TrigonalPrismatic => &TRIGONAL_PRISMATIC,
    }
}

/// The proper-rotation group both two-substituent geometries share: the identity, and
/// the half-turn about the bisector of the two directions, which exchanges them
/// whatever angle they subtend. Neither geometry is therefore stereogenic.
const PAIR: &[&[u8]] = &[&[0, 1], &[1, 0]];

static LINEAR: Idealization = Idealization {
    group: PAIR,
    reflection: &[0, 1],
};

static ANGULAR: Idealization = Idealization {
    group: PAIR,
    reflection: &[0, 1],
};

static TRIGONAL_PLANAR: Idealization = Idealization {
    group: &[
        &[0, 1, 2],
        &[0, 2, 1],
        &[1, 0, 2],
        &[1, 2, 0],
        &[2, 0, 1],
        &[2, 1, 0],
    ],
    reflection: &[0, 1, 2],
};

static TRIGONAL_PYRAMIDAL: Idealization = Idealization {
    group: &[&[0, 1, 2], &[1, 2, 0], &[2, 0, 1]],
    reflection: &[0, 2, 1],
};

static T_SHAPED: Idealization = Idealization {
    group: &[&[0, 1, 2], &[1, 0, 2]],
    reflection: &[0, 1, 2],
};

static TETRAHEDRAL: Idealization = Idealization {
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
};

static SQUARE_PLANAR: Idealization = Idealization {
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
};

static PYRAMIDALIZED_SQUARE: Idealization = Idealization {
    group: &[&[0, 1, 2, 3], &[1, 2, 3, 0], &[2, 3, 0, 1], &[3, 0, 1, 2]],
    reflection: &[0, 3, 2, 1],
};

static SEESAW: Idealization = Idealization {
    group: &[&[0, 1, 2, 3], &[1, 0, 3, 2]],
    reflection: &[0, 1, 3, 2],
};

static TRIGONAL_BIPYRAMIDAL: Idealization = Idealization {
    group: &[
        &[0, 1, 2, 3, 4],
        &[0, 1, 3, 4, 2],
        &[0, 1, 4, 2, 3],
        &[1, 0, 2, 4, 3],
        &[1, 0, 3, 2, 4],
        &[1, 0, 4, 3, 2],
    ],
    reflection: &[1, 0, 2, 3, 4],
};

static SQUARE_PYRAMIDAL: Idealization = Idealization {
    group: &[
        &[0, 1, 2, 3, 4],
        &[0, 2, 3, 4, 1],
        &[0, 3, 4, 1, 2],
        &[0, 4, 1, 2, 3],
    ],
    reflection: &[0, 1, 4, 3, 2],
};

static OCTAHEDRAL: Idealization = Idealization {
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
};

static TRIGONAL_PRISMATIC: Idealization = Idealization {
    group: &[
        &[0, 1, 2, 3, 4, 5],
        &[1, 2, 0, 4, 5, 3],
        &[2, 0, 1, 5, 3, 4],
        &[3, 5, 4, 0, 2, 1],
        &[4, 3, 5, 1, 0, 2],
        &[5, 4, 3, 2, 1, 0],
    ],
    reflection: &[0, 2, 1, 3, 5, 4],
};

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

/// Whether two slices hold the same bytes — a `const` slice equality, for want
/// of one in the standard library.
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

    use super::CoordinationGeometry::*;

    const GEOMETRIES: [CoordinationGeometry; 13] = [
        Linear,
        Angular,
        TrigonalPlanar,
        TrigonalPyramidal,
        TShaped,
        Tetrahedral,
        SquarePlanar,
        PyramidalizedSquare,
        Seesaw,
        TrigonalBipyramidal,
        SquarePyramidal,
        Octahedral,
        TrigonalPrismatic,
    ];

    const R: f64 = 0.866_025_403_784_438_6;

    fn placement(geometry: CoordinationGeometry) -> Vec<[f64; 3]> {
        match geometry {
            Linear => vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
            Angular => angular(1.2),
            TrigonalPlanar => vec![[1.0, 0.0, 0.0], [-0.5, R, 0.0], [-0.5, -R, 0.0]],
            TrigonalPyramidal => vec![[1.0, -1.0, -1.0], [-1.0, -1.0, 1.0], [-1.0, 1.0, -1.0]],
            TShaped => vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]],
            Tetrahedral => vec![
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [-1.0, 1.0, -1.0],
            ],
            SquarePlanar => vec![
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
            ],
            PyramidalizedSquare => pyramidalized_square(1.2),
            Seesaw => vec![
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [1.0, 0.0, 0.0],
                [-0.5, R, 0.0],
            ],
            TrigonalBipyramidal => vec![
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [1.0, 0.0, 0.0],
                [-0.5, R, 0.0],
                [-0.5, -R, 0.0],
            ],
            SquarePyramidal => vec![
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
            ],
            Octahedral => vec![
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ],
            TrigonalPrismatic => vec![
                [1.0, 0.0, 1.0],
                [-0.5, R, 1.0],
                [-0.5, -R, 1.0],
                [1.0, 0.0, -1.0],
                [-0.5, R, -1.0],
                [-0.5, -R, -1.0],
            ],
        }
    }

    fn angular(angle: f64) -> Vec<[f64; 3]> {
        let (half_sine, half_cosine) = (angle / 2.0).sin_cos();
        vec![
            [half_sine, 0.0, half_cosine],
            [-half_sine, 0.0, half_cosine],
        ]
    }

    fn pyramidalized_square(angle: f64) -> Vec<[f64; 3]> {
        let (sine, cosine) = angle.sin_cos();
        vec![
            [sine, 0.0, cosine],
            [0.0, sine, cosine],
            [-sine, 0.0, cosine],
            [0.0, -sine, cosine],
        ]
    }

    fn identity(n: usize) -> Vec<u8> {
        (0..n as u8).collect()
    }

    fn is_permutation(permutation: &[u8], n: usize) -> bool {
        let mut sorted = permutation.to_vec();
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

    fn preserves_angles(placement: &[[f64; 3]], permutation: &[u8]) -> bool {
        let n = placement.len();
        (0..n).all(|i| {
            (0..n).all(|j| {
                let moved = dot(
                    placement[permutation[i] as usize],
                    placement[permutation[j] as usize],
                );
                (moved - dot(placement[i], placement[j])).abs() < 1e-9
            })
        })
    }

    fn spanning_triple(placement: &[[f64; 3]]) -> Option<(usize, usize, usize)> {
        let n = placement.len();
        (0..n)
            .flat_map(|a| (a + 1..n).flat_map(move |b| (b + 1..n).map(move |c| (a, b, c))))
            .find(|&(a, b, c)| signed_volume(placement[a], placement[b], placement[c]).abs() > 1e-9)
    }

    fn is_rotation(placement: &[[f64; 3]], permutation: &[u8]) -> bool {
        if !preserves_angles(placement, permutation) {
            return false;
        }
        match spanning_triple(placement) {
            None => true,
            Some((a, b, c)) => {
                let before = signed_volume(placement[a], placement[b], placement[c]);
                let after = signed_volume(
                    placement[permutation[a] as usize],
                    placement[permutation[b] as usize],
                    placement[permutation[c] as usize],
                );
                (before > 0.0) == (after > 0.0)
            }
        }
    }

    fn rotation_group(placement: &[[f64; 3]]) -> Vec<Vec<u8>> {
        let mut rotations: Vec<Vec<u8>> = permutations(placement.len())
            .into_iter()
            .filter(|permutation| is_rotation(placement, permutation))
            .collect();
        rotations.sort_unstable();
        rotations
    }

    fn sorted_group(geometry: CoordinationGeometry) -> Vec<Vec<u8>> {
        let mut group: Vec<Vec<u8>> = idealization(geometry)
            .group
            .iter()
            .map(|element| element.to_vec())
            .collect();
        group.sort_unstable();
        group
    }

    #[test]
    fn slot_count_is_the_substituent_count_of_the_geometry() {
        assert_eq!(Linear.slot_count(), 2);
        assert_eq!(Angular.slot_count(), 2);
        assert_eq!(TrigonalPlanar.slot_count(), 3);
        assert_eq!(TrigonalPyramidal.slot_count(), 3);
        assert_eq!(TShaped.slot_count(), 3);
        assert_eq!(Tetrahedral.slot_count(), 4);
        assert_eq!(SquarePlanar.slot_count(), 4);
        assert_eq!(PyramidalizedSquare.slot_count(), 4);
        assert_eq!(Seesaw.slot_count(), 4);
        assert_eq!(TrigonalBipyramidal.slot_count(), 5);
        assert_eq!(SquarePyramidal.slot_count(), 5);
        assert_eq!(Octahedral.slot_count(), 6);
        assert_eq!(TrigonalPrismatic.slot_count(), 6);
    }

    #[test]
    fn configuration_count_is_the_number_of_distinct_stereoisomers() {
        assert_eq!(Linear.configuration_count(), 1);
        assert_eq!(Angular.configuration_count(), 1);
        assert_eq!(TrigonalPlanar.configuration_count(), 1);
        assert_eq!(TrigonalPyramidal.configuration_count(), 2);
        assert_eq!(TShaped.configuration_count(), 3);
        assert_eq!(Tetrahedral.configuration_count(), 2);
        assert_eq!(SquarePlanar.configuration_count(), 3);
        assert_eq!(PyramidalizedSquare.configuration_count(), 6);
        assert_eq!(Seesaw.configuration_count(), 12);
        assert_eq!(TrigonalBipyramidal.configuration_count(), 20);
        assert_eq!(SquarePyramidal.configuration_count(), 30);
        assert_eq!(Octahedral.configuration_count(), 30);
        assert_eq!(TrigonalPrismatic.configuration_count(), 120);
    }

    #[test]
    fn chiral_geometries_are_chiral() {
        assert!(TrigonalPyramidal.is_chiral());
        assert!(Tetrahedral.is_chiral());
        assert!(PyramidalizedSquare.is_chiral());
        assert!(Seesaw.is_chiral());
        assert!(TrigonalBipyramidal.is_chiral());
        assert!(SquarePyramidal.is_chiral());
        assert!(Octahedral.is_chiral());
        assert!(TrigonalPrismatic.is_chiral());
    }

    #[test]
    fn achiral_geometries_are_not_chiral() {
        assert!(!Linear.is_chiral());
        assert!(!Angular.is_chiral());
        assert!(!TrigonalPlanar.is_chiral());
        assert!(!TShaped.is_chiral());
        assert!(!SquarePlanar.is_chiral());
    }

    #[test]
    fn geometries_order_by_slot_count_then_configuration_count() {
        assert!(GEOMETRIES.is_sorted());
        assert!(
            GEOMETRIES.is_sorted_by_key(|geometry| (
                geometry.slot_count(),
                geometry.configuration_count()
            ))
        );
    }

    #[test]
    fn every_group_permutes_its_slots() {
        for geometry in GEOMETRIES {
            for &element in idealization(geometry).group {
                assert!(
                    is_permutation(element, geometry.slot_count()),
                    "{geometry:?}"
                );
            }
        }
    }

    #[test]
    fn every_group_contains_the_identity() {
        for geometry in GEOMETRIES {
            let identity = identity(geometry.slot_count());
            assert!(
                idealization(geometry).group.contains(&identity.as_slice()),
                "{geometry:?}"
            );
        }
    }

    #[test]
    fn every_group_is_closed_under_composition() {
        for geometry in GEOMETRIES {
            let group = idealization(geometry).group;
            for &g in group {
                for &h in group {
                    let product = compose(g, h);
                    assert!(group.contains(&product.as_slice()), "{geometry:?}");
                }
            }
        }
    }

    #[test]
    fn every_group_is_the_rotation_group_of_the_geometry() {
        for geometry in GEOMETRIES {
            assert_eq!(
                sorted_group(geometry),
                rotation_group(&placement(geometry)),
                "{geometry:?}"
            );
        }
    }

    #[test]
    fn a_family_has_one_group_throughout() {
        for angle in [0.6, 1.2, 1.8, 3.0] {
            assert_eq!(rotation_group(&angular(angle)), sorted_group(Angular));
        }
        for angle in [0.3, 0.9, 1.4] {
            assert_eq!(
                rotation_group(&pyramidalized_square(angle)),
                sorted_group(PyramidalizedSquare)
            );
        }
    }

    #[test]
    fn every_reflection_permutes_its_slots() {
        for geometry in GEOMETRIES {
            assert!(
                is_permutation(idealization(geometry).reflection, geometry.slot_count()),
                "{geometry:?}"
            );
        }
    }

    #[test]
    fn every_reflection_preserves_the_pairwise_angles() {
        for geometry in GEOMETRIES {
            assert!(
                preserves_angles(&placement(geometry), idealization(geometry).reflection),
                "{geometry:?}"
            );
        }
    }
}

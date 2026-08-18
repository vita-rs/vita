/// The symmetry an idealized arrangement of substituents fixes, and the reference
/// frame realizing it — the one datum every fact about an arrangement is a function
/// of, so a further arrangement is new data here, not new logic.
mod idealization {
    /// The largest slot count any arrangement holds — six, the octahedron or
    /// trigonal prism. The fixed-size slot buffers hold this many; unused entries
    /// stay zero.
    pub const MAX_SLOTS: usize = 6;

    /// The symmetry of an arrangement over its slots.
    pub struct Symmetry {
        /// The proper-rotation group over the slots — the orderings a configuration
        /// treats as equivalent, the closure of the arrangement's rotation
        /// generators.
        pub group: &'static [&'static [u8]],
        /// The slot permutation a reflection induces, carrying a configuration onto
        /// its mirror image: outside `group` for a chiral arrangement, inside it for
        /// an achiral one, so the mirror is well defined on cosets.
        pub reflection: &'static [u8],
    }

    impl Symmetry {
        /// The number of slots — the length of any group element.
        pub const fn slot_count(&self) -> usize {
            self.reflection.len()
        }

        /// Whether the arrangement is chiral: its reflection is not itself a
        /// rotation.
        pub const fn is_chiral(&self) -> bool {
            let mut index = 0;
            while index < self.group.len() {
                if slices_equal(self.group[index], self.reflection) {
                    return false;
                }
                index += 1;
            }
            true
        }

        /// The distinct configurations the all-different case admits over `ends`
        /// independent ends: the orderings that respect the ends, quotiented by the
        /// rotation group. The presentations number `ends! · (slots / ends)!^ends` —
        /// the ends ordered, then each filled — and the group acts on them freely
        /// when every substituent differs, so the quotient is their count over the
        /// group's order.
        pub const fn configuration_count(&self, ends: usize) -> usize {
            let per_end = self.slot_count() / ends;
            factorial(ends) * factorial(per_end).pow(ends as u32) / self.group.len()
        }
    }

    /// A symmetry together with the directions realizing it, which a coordination
    /// sphere has and a rigid two-ended frame does not.
    pub struct Idealization {
        /// The symmetry over the slots.
        pub symmetry: Symmetry,
        /// The reference directions, slot by slot, a site's coordinates align onto.
        ///
        /// Not normalized: each arrangement is stated in whichever vectors put it
        /// plainest, so a consumer comparing angles normalizes them first.
        pub directions: &'static [[f64; 3]],
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

    /// Reading a symmetry off a placement of the slots in space — the independent
    /// oracle every declared [`Symmetry`] is checked against.
    #[cfg(test)]
    pub mod oracle {
        /// The identity permutation of `n` slots.
        pub fn identity(n: usize) -> Vec<u8> {
            (0..n as u8).collect()
        }

        /// Whether `permutation` permutes `n` slots.
        pub fn is_permutation(permutation: &[u8], n: usize) -> bool {
            let mut sorted = permutation.to_vec();
            sorted.sort_unstable();
            sorted == identity(n)
        }

        /// `before` followed by `after`.
        pub fn compose(after: &[u8], before: &[u8]) -> Vec<u8> {
            before.iter().map(|&index| after[index as usize]).collect()
        }

        /// A group's elements, sorted, so two groups compare as sets.
        pub fn sorted_group(group: &[&[u8]]) -> Vec<Vec<u8>> {
            let mut sorted: Vec<Vec<u8>> = group.iter().map(|element| element.to_vec()).collect();
            sorted.sort_unstable();
            sorted
        }

        /// The permutations of `placement` that are proper rotations of it — the
        /// rotation group the placement realizes, derived from the geometry alone.
        pub fn rotation_group(placement: &[[f64; 3]]) -> Vec<Vec<u8>> {
            let mut rotations: Vec<Vec<u8>> = permutations(placement.len())
                .into_iter()
                .filter(|permutation| is_rotation(placement, permutation))
                .collect();
            rotations.sort_unstable();
            rotations
        }

        /// Whether `permutation` leaves every pairwise angle of `placement` intact.
        pub fn preserves_angles(placement: &[[f64; 3]], permutation: &[u8]) -> bool {
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

        /// Whether `permutation` is a proper rotation of `placement`: an isometry that
        /// also keeps the orientation of a spanning triple.
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

        /// The first slot triple `placement` spans with nonzero volume, or `None` if it
        /// is coplanar with the origin.
        fn spanning_triple(placement: &[[f64; 3]]) -> Option<(usize, usize, usize)> {
            let n = placement.len();
            (0..n)
                .flat_map(|a| (a + 1..n).flat_map(move |b| (b + 1..n).map(move |c| (a, b, c))))
                .find(|&(a, b, c)| {
                    signed_volume(placement[a], placement[b], placement[c]).abs() > 1e-9
                })
        }

        /// The volume the three vectors span, signed by their handedness.
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

        /// The inner product of two vectors.
        fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
            a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        }

        /// Every permutation of `n` slots.
        fn permutations(n: usize) -> Vec<Vec<u8>> {
            let mut result = Vec::new();
            permute(&mut identity(n), 0, &mut result);
            result
        }

        /// Appends to `result` every permutation of `slice` that holds its first
        /// `start` elements fixed.
        fn permute(slice: &mut [u8], start: usize, result: &mut Vec<Vec<u8>>) {
            if start == slice.len() {
                result.push(slice.to_vec());
                return;
            }
            for index in start..slice.len() {
                slice.swap(start, index);
                permute(slice, start + 1, result);
                slice.swap(start, index);
            }
        }
    }
}

/// The arrangements a coordination sphere takes, named as the standard names them.
mod coordination_geometry {
    use super::idealization::{Idealization, Symmetry};

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
        /// Every geometry, in declaration order.
        pub(crate) const ALL: [Self; 13] = [
            Self::Linear,
            Self::Angular,
            Self::TrigonalPlanar,
            Self::TrigonalPyramidal,
            Self::TShaped,
            Self::Tetrahedral,
            Self::SquarePlanar,
            Self::PyramidalizedSquare,
            Self::Seesaw,
            Self::TrigonalBipyramidal,
            Self::SquarePyramidal,
            Self::Octahedral,
            Self::TrigonalPrismatic,
        ];

        /// Returns the number of slots the geometry arranges — its coordination
        /// number.
        #[inline]
        pub const fn slot_count(self) -> usize {
            self.symmetry().slot_count()
        }

        /// Returns the number of distinct configurations the geometry admits — the
        /// stereoisomers its slots realize when every substituent differs.
        ///
        /// One means the geometry equates every ordering of its substituents, and
        /// so carries no stereochemistry however it is substituted.
        #[inline]
        pub const fn configuration_count(self) -> usize {
            self.symmetry().configuration_count(1)
        }

        /// Returns whether the geometry is chiral — whether a configuration and its
        /// mirror image are distinct, no rotation carrying one onto the other.
        #[inline]
        pub const fn is_chiral(self) -> bool {
            self.symmetry().is_chiral()
        }

        /// The reference directions its substituents align onto, slot by slot.
        ///
        /// Not normalized — each arrangement is stated in whichever vectors put it
        /// plainest, so a consumer comparing angles normalizes them first.
        #[inline]
        pub(crate) const fn directions(self) -> &'static [[f64; 3]] {
            self.idealization().directions
        }

        /// The idealization the geometry denotes.
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
        ///
        /// The reference directions of the geometries a parent polyhedron subsumes
        /// are that parent's, less the vertices the substituents leave vacant.
        pub(super) const fn idealization(self) -> &'static Idealization {
            match self {
                Self::Linear => &LINEAR,
                Self::Angular => &ANGULAR,
                Self::TrigonalPlanar => &TRIGONAL_PLANAR,
                Self::TrigonalPyramidal => &TRIGONAL_PYRAMIDAL,
                Self::TShaped => &T_SHAPED,
                Self::Tetrahedral => &TETRAHEDRAL,
                Self::SquarePlanar => &SQUARE_PLANAR,
                Self::PyramidalizedSquare => &PYRAMIDALIZED_SQUARE,
                Self::Seesaw => &SEESAW,
                Self::TrigonalBipyramidal => &TRIGONAL_BIPYRAMIDAL,
                Self::SquarePyramidal => &SQUARE_PYRAMIDAL,
                Self::Octahedral => &OCTAHEDRAL,
                Self::TrigonalPrismatic => &TRIGONAL_PRISMATIC,
            }
        }

        /// The symmetry the geometry fixes.
        const fn symmetry(self) -> &'static Symmetry {
            &self.idealization().symmetry
        }
    }

    /// The proper-rotation group both two-substituent geometries share: the identity,
    /// and the half-turn about the bisector of the two directions, which exchanges
    /// them whatever angle they subtend. Neither geometry is therefore stereogenic.
    const PAIR: &[&[u8]] = &[&[0, 1], &[1, 0]];

    /// √3 / 2 — the in-plane offset of a trigonal vertex.
    const R: f64 = 0.866_025_403_784_438_6;

    static LINEAR: Idealization = Idealization {
        symmetry: Symmetry {
            group: PAIR,
            reflection: &[0, 1],
        },
        directions: &[[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
    };

    static ANGULAR: Idealization = Idealization {
        symmetry: Symmetry {
            group: PAIR,
            reflection: &[0, 1],
        },
        directions: &[[1.0, 1.0, 1.0], [1.0, -1.0, -1.0]],
    };

    static TRIGONAL_PLANAR: Idealization = Idealization {
        symmetry: Symmetry {
            group: &[
                &[0, 1, 2],
                &[0, 2, 1],
                &[1, 0, 2],
                &[1, 2, 0],
                &[2, 0, 1],
                &[2, 1, 0],
            ],
            reflection: &[0, 1, 2],
        },
        directions: &[[1.0, 0.0, 0.0], [-0.5, R, 0.0], [-0.5, -R, 0.0]],
    };

    static TRIGONAL_PYRAMIDAL: Idealization = Idealization {
        symmetry: Symmetry {
            group: &[&[0, 1, 2], &[1, 2, 0], &[2, 0, 1]],
            reflection: &[0, 2, 1],
        },
        directions: &[[1.0, -1.0, -1.0], [-1.0, -1.0, 1.0], [-1.0, 1.0, -1.0]],
    };

    static T_SHAPED: Idealization = Idealization {
        symmetry: Symmetry {
            group: &[&[0, 1, 2], &[1, 0, 2]],
            reflection: &[0, 1, 2],
        },
        directions: &[[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]],
    };

    static TETRAHEDRAL: Idealization = Idealization {
        symmetry: Symmetry {
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
        },
        directions: &[
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
        ],
    };

    static SQUARE_PLANAR: Idealization = Idealization {
        symmetry: Symmetry {
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
        },
        directions: &[
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
    };

    static PYRAMIDALIZED_SQUARE: Idealization = Idealization {
        symmetry: Symmetry {
            group: &[&[0, 1, 2, 3], &[1, 2, 3, 0], &[2, 3, 0, 1], &[3, 0, 1, 2]],
            reflection: &[0, 3, 2, 1],
        },
        directions: &[
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [-1.0, 0.0, 1.0],
            [0.0, -1.0, 1.0],
        ],
    };

    static SEESAW: Idealization = Idealization {
        symmetry: Symmetry {
            group: &[&[0, 1, 2, 3], &[1, 0, 3, 2]],
            reflection: &[0, 1, 3, 2],
        },
        directions: &[
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-0.5, R, 0.0],
        ],
    };

    static TRIGONAL_BIPYRAMIDAL: Idealization = Idealization {
        symmetry: Symmetry {
            group: &[
                &[0, 1, 2, 3, 4],
                &[0, 1, 3, 4, 2],
                &[0, 1, 4, 2, 3],
                &[1, 0, 2, 4, 3],
                &[1, 0, 3, 2, 4],
                &[1, 0, 4, 3, 2],
            ],
            reflection: &[1, 0, 2, 3, 4],
        },
        directions: &[
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-0.5, R, 0.0],
            [-0.5, -R, 0.0],
        ],
    };

    static SQUARE_PYRAMIDAL: Idealization = Idealization {
        symmetry: Symmetry {
            group: &[
                &[0, 1, 2, 3, 4],
                &[0, 2, 3, 4, 1],
                &[0, 3, 4, 1, 2],
                &[0, 4, 1, 2, 3],
            ],
            reflection: &[0, 1, 4, 3, 2],
        },
        directions: &[
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
    };

    static OCTAHEDRAL: Idealization = Idealization {
        symmetry: Symmetry {
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
        },
        directions: &[
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
    };

    static TRIGONAL_PRISMATIC: Idealization = Idealization {
        symmetry: Symmetry {
            group: &[
                &[0, 1, 2, 3, 4, 5],
                &[1, 2, 0, 4, 5, 3],
                &[2, 0, 1, 5, 3, 4],
                &[3, 5, 4, 0, 2, 1],
                &[4, 3, 5, 1, 0, 2],
                &[5, 4, 3, 2, 1, 0],
            ],
            reflection: &[0, 2, 1, 3, 5, 4],
        },
        directions: &[
            [1.0, 0.0, 1.0],
            [-0.5, R, 1.0],
            [-0.5, -R, 1.0],
            [1.0, 0.0, -1.0],
            [-0.5, R, -1.0],
            [-0.5, -R, -1.0],
        ],
    };

    #[cfg(test)]
    mod tests {
        use super::*;

        use super::super::idealization::{MAX_SLOTS, oracle::*};
        use super::CoordinationGeometry::*;

        const GEOMETRIES: [CoordinationGeometry; 13] = CoordinationGeometry::ALL;

        fn group(geometry: CoordinationGeometry) -> Vec<Vec<u8>> {
            sorted_group(geometry.idealization().symmetry.group)
        }

        fn reflection(geometry: CoordinationGeometry) -> &'static [u8] {
            geometry.idealization().symmetry.reflection
        }

        fn directions(geometry: CoordinationGeometry) -> &'static [[f64; 3]] {
            geometry.idealization().directions
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
            assert!(GEOMETRIES.is_sorted_by_key(|geometry| (
                geometry.slot_count(),
                geometry.configuration_count()
            )));
        }

        #[test]
        fn no_geometry_holds_more_slots_than_the_buffer() {
            for geometry in GEOMETRIES {
                assert!(geometry.slot_count() <= MAX_SLOTS, "{geometry:?}");
            }
        }

        #[test]
        fn every_group_permutes_its_slots() {
            for geometry in GEOMETRIES {
                for element in group(geometry) {
                    assert!(
                        is_permutation(&element, geometry.slot_count()),
                        "{geometry:?}"
                    );
                }
            }
        }

        #[test]
        fn every_group_contains_the_identity() {
            for geometry in GEOMETRIES {
                let identity = identity(geometry.slot_count());
                assert!(group(geometry).contains(&identity), "{geometry:?}");
            }
        }

        #[test]
        fn every_group_is_closed_under_composition() {
            for geometry in GEOMETRIES {
                let elements = group(geometry);
                for g in &elements {
                    for h in &elements {
                        assert!(elements.contains(&compose(g, h)), "{geometry:?}");
                    }
                }
            }
        }

        #[test]
        fn every_group_is_the_rotation_group_of_the_reference_directions() {
            for geometry in GEOMETRIES {
                assert_eq!(
                    group(geometry),
                    rotation_group(directions(geometry)),
                    "{geometry:?}"
                );
            }
        }

        #[test]
        fn a_family_has_one_group_throughout() {
            for angle in [0.6, 1.2, 1.8, 3.0] {
                assert_eq!(rotation_group(&angular(angle)), group(Angular));
            }
            for angle in [0.3, 0.9, 1.4] {
                assert_eq!(
                    rotation_group(&pyramidalized_square(angle)),
                    group(PyramidalizedSquare)
                );
            }
        }

        #[test]
        fn every_reflection_permutes_its_slots() {
            for geometry in GEOMETRIES {
                assert!(
                    is_permutation(reflection(geometry), geometry.slot_count()),
                    "{geometry:?}"
                );
            }
        }

        #[test]
        fn every_reflection_preserves_the_pairwise_angles() {
            for geometry in GEOMETRIES {
                assert!(
                    preserves_angles(directions(geometry), reflection(geometry)),
                    "{geometry:?}"
                );
            }
        }

        #[test]
        fn reference_directions_number_the_slots() {
            for geometry in GEOMETRIES {
                assert_eq!(
                    directions(geometry).len(),
                    geometry.slot_count(),
                    "{geometry:?}"
                );
            }
        }
    }
}

/// Stereochemistry's vocabulary: where a stereogenic unit sits, which arrangement it
/// carries, and which of that arrangement's configurations it realizes.
mod stereo {
    use vita_core::SiteId;

    use super::coordination_geometry::CoordinationGeometry;
    use super::idealization::{MAX_SLOTS, Symmetry};
    use crate::BondId;

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
                (Self::Site(_), StereoKind::Center(_))
                    | (Self::Bond(_), StereoKind::Bond)
                    | (Self::Axis(_), StereoKind::Axis)
            )
        }
    }

    /// A [`CoordinationGeometry`] that admits more than one configuration.
    ///
    /// The geometries a rotation equates every ordering of —
    /// [`Linear`](CoordinationGeometry::Linear),
    /// [`Angular`](CoordinationGeometry::Angular),
    /// [`TrigonalPlanar`](CoordinationGeometry::TrigonalPlanar) — are excluded: they
    /// carry no stereochemistry however they are substituted, so no configuration
    /// distinguishes anything at a site bearing one.
    ///
    /// Admitting configurations is a property of the geometry, not of a site: a
    /// carbon is [`Tetrahedral`](CoordinationGeometry::Tetrahedral) in methane too.
    /// Whether more than one is realizable turns on the substituents, which only
    /// [`stereocenters`](crate::stereo::stereocenters) settles.
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct StereogenicGeometry(CoordinationGeometry);

    impl StereogenicGeometry {
        /// Constructs a `StereogenicGeometry` from `geometry`, returning `None` if it
        /// admits a single configuration.
        #[inline]
        pub const fn new(geometry: CoordinationGeometry) -> Option<Self> {
            if geometry.configuration_count() > 1 {
                Some(Self(geometry))
            } else {
                None
            }
        }

        /// Returns the geometry.
        #[inline]
        pub const fn geometry(self) -> CoordinationGeometry {
            self.0
        }
    }

    /// The kind of a stereogenic unit: the arrangement its anchor carries, whose
    /// rotation group fixes which of its neighbor orderings are equivalent.
    ///
    /// A kind is a pure data key. It selects the permutation group under which a
    /// configuration's neighbor ordering reduces, and the geometric reference against
    /// which coordinates are perceived. A center's arrangement is whichever
    /// coordination geometry it takes; a bond's and an axis's are fixed — four
    /// substituents on a rigid two-ended frame, told apart by the mirror each admits.
    ///
    /// The centers order among themselves by their geometry, then the bond, then the
    /// axis.
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub enum StereoKind {
        /// A coordination center of the given geometry.
        Center(StereogenicGeometry),
        /// A double bond (4 neighbors, 2 configurations).
        Bond,
        /// An allene axis (4 neighbors, 2 configurations).
        Axis,
    }

    impl StereoKind {
        /// Returns the number of neighbor slots the arrangement holds.
        #[inline]
        pub const fn slot_count(self) -> usize {
            self.symmetry().slot_count()
        }

        /// Returns the number of distinct configurations the arrangement admits — the
        /// stereoisomers its slots realize when every substituent differs.
        #[inline]
        pub const fn configuration_count(self) -> usize {
            self.symmetry().configuration_count(self.ends())
        }

        /// Returns whether the arrangement is chiral — whether a configuration and its
        /// mirror image are distinct, no rotation carrying one onto the other.
        #[inline]
        pub const fn is_chiral(self) -> bool {
            self.symmetry().is_chiral()
        }

        /// The reference directions its coordinates align onto, slot by slot — empty
        /// for an arrangement read from a dihedral instead, a bond or an axis.
        ///
        /// Not normalized — see [`CoordinationGeometry::directions`].
        #[inline]
        pub(crate) const fn directions(self) -> &'static [[f64; 3]] {
            match self {
                Self::Center(stereogenic) => stereogenic.geometry().directions(),
                Self::Bond | Self::Axis => &[],
            }
        }

        /// The slot permutation a reflection induces, carrying a configuration onto its
        /// mirror image.
        #[inline]
        pub(crate) const fn reflection(self) -> &'static [u8] {
            self.symmetry().reflection
        }

        /// The number of independent ends the slots split into — one freely permuted
        /// center, or the two ends of a rigid frame whose substituents cannot cross.
        #[inline]
        pub(crate) const fn ends(self) -> usize {
            match self {
                Self::Center(_) => 1,
                Self::Bond | Self::Axis => 2,
            }
        }

        /// The symmetry the arrangement fixes.
        const fn symmetry(self) -> &'static Symmetry {
            match self {
                Self::Center(stereogenic) => &stereogenic.geometry().idealization().symmetry,
                Self::Bond => &BOND,
                Self::Axis => &AXIS,
            }
        }
    }

    /// A stereodescriptor: which of an arrangement's configurations a unit realizes,
    /// read against a ranking of its substituents.
    ///
    /// A configuration is a coset of neighbor orderings under the arrangement's
    /// rotation group — the orderings a rotation cannot tell apart. Its descriptor is
    /// that coset's canonical representative, computed from the neighbors' *ranks*
    /// rather than the neighbors themselves, so symmetry-equivalent units reduce alike
    /// and two units of the same kind are comparable. It is opaque and relative: it
    /// distinguishes a configuration from the arrangement's others and from its
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
        pub const fn kind(self) -> StereoKind {
            self.kind
        }

        /// The descriptor of the mirror-image configuration.
        ///
        /// A chiral kind's descriptor flips to a distinct one; an achiral kind's is its
        /// own mirror, so the two coincide.
        #[inline]
        pub fn mirror(self) -> StereoDescriptor {
            let symmetry = self.kind.symmetry();
            StereoDescriptor {
                kind: self.kind,
                coset: reduce(apply(symmetry.reflection, self.coset), symmetry.group),
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
    /// Below, `a` is the anchoring site, and a center's slots are those of the
    /// coordination geometry named.
    ///
    /// - [`TrigonalPyramidal`](CoordinationGeometry::TrigonalPyramidal) — the three
    ///   substituents `[n0, n1, n2]`; reference
    ///   `(n0 − a) · ((n1 − a) × (n2 − a))`.
    /// - [`TShaped`](CoordinationGeometry::TShaped) — the two opposed substituents,
    ///   then the one across them, `[t0, t1, c]`.
    /// - [`Tetrahedral`](CoordinationGeometry::Tetrahedral) — the four substituents
    ///   `[n0, n1, n2, n3]`; reference `(n1 − n0) · ((n2 − n0) × (n3 − n0))`.
    /// - [`SquarePlanar`](CoordinationGeometry::SquarePlanar) — the four in cyclic
    ///   order `[n0, n1, n2, n3]`, `n0`/`n2` and `n1`/`n3` trans.
    /// - [`PyramidalizedSquare`](CoordinationGeometry::PyramidalizedSquare) — the four
    ///   in cyclic order `[n0, n1, n2, n3]`, `n0`/`n2` and `n1`/`n3` across; reference
    ///   `(n0 − a) · ((n1 − a) × (n2 − a))`.
    /// - [`Seesaw`](CoordinationGeometry::Seesaw) — the two axial, then the two
    ///   equatorial, `[a0, a1, e0, e1]`; reference `(a0 − a1) · ((e0 − a) × (e1 − a))`.
    /// - [`TrigonalBipyramidal`](CoordinationGeometry::TrigonalBipyramidal) — the two
    ///   axial, then the three equatorial, `[a0, a1, e0, e1, e2]`; reference
    ///   `(a0 − a1) · ((e1 − e0) × (e2 − e0))`.
    /// - [`SquarePyramidal`](CoordinationGeometry::SquarePyramidal) — the apical, then
    ///   the four basal in cyclic order, `[p, b0, b1, b2, b3]`; reference
    ///   `(p − b0) · ((b1 − b0) × (b2 − b0))`.
    /// - [`Octahedral`](CoordinationGeometry::Octahedral) — three trans pairs
    ///   `[n0, n1, n2, n3, n4, n5]`, `n2i` opposite `n2i+1`; reference
    ///   `(n0 − n1) · ((n2 − n3) × (n4 − n5))`.
    /// - [`TrigonalPrismatic`](CoordinationGeometry::TrigonalPrismatic) — one
    ///   triangular face, then the other, `[t0, t1, t2, b0, b1, b2]`, `bi` eclipsing
    ///   `ti`; reference `(t0 − b0) · ((t1 − t0) × (t2 − t0))`.
    /// - [`Bond`](StereoKind::Bond) — one end's two, then the other's,
    ///   `[e1a, e1b, e2a, e2b]`, `e1a` on `e2a`'s side.
    /// - [`Axis`](StereoKind::Axis) — one terminus's two, then the other's,
    ///   `[t1a, t1b, t2a, t2b]`; reference the twist across the axis, its mirror
    ///   swapping `t1a` and `t1b`.
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
                neighbors: normalized(kind, &neighbors),
            })
        }

        /// The anchor of the stereogenic unit.
        #[inline]
        pub const fn locus(&self) -> StereoLocus {
            self.locus
        }

        /// The arrangement fixing which orderings are equivalent.
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
            StereoDescriptor {
                kind: self.kind,
                coset: reduce(
                    relative_order(&self.neighbors, rank),
                    self.kind.symmetry().group,
                ),
            }
        }
    }

    /// The proper-rotation group both rigid two-ended frames share — a double bond, an
    /// allene: the identity and the three half-turns of a rigid four-substituent frame
    /// (about its axis, and about each of the two axes across it, one of which swaps
    /// the ends). Chirality then turns on the reflection alone.
    const EDGE: &[&[u8]] = &[&[0, 1, 2, 3], &[1, 0, 3, 2], &[2, 3, 0, 1], &[3, 2, 1, 0]];

    /// A double bond: the in-plane mirror leaves the frame alone, so it is achiral.
    static BOND: Symmetry = Symmetry {
        group: EDGE,
        reflection: &[0, 1, 2, 3],
    };

    /// An allene axis: the diagonal mirror swaps one terminus's substituents, so it is
    /// chiral.
    static AXIS: Symmetry = Symmetry {
        group: EDGE,
        reflection: &[1, 0, 2, 3],
    };

    /// The lexicographically least neighbor ordering in a configuration's rotation
    /// orbit — the canonical presentation, so equal configurations store equal
    /// neighbors and `==`, `Ord`, and `Hash` follow from the stored order. Only proper
    /// rotations act, leaving intact the handedness a chiral kind fixes.
    fn normalized(kind: StereoKind, neighbors: &[SiteId]) -> Vec<SiteId> {
        kind.symmetry()
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

    #[cfg(test)]
    mod tests {
        use super::*;

        use super::super::idealization::oracle::*;
        use super::CoordinationGeometry::*;

        const DEGENERATE: [CoordinationGeometry; 3] = [Linear, Angular, TrigonalPlanar];

        const STEREOGENIC: [CoordinationGeometry; 10] = [
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

        fn s(n: u32) -> SiteId {
            SiteId::new(n).unwrap()
        }

        fn b(n: u32) -> BondId {
            BondId::new(n).unwrap()
        }

        fn center(geometry: CoordinationGeometry) -> StereoKind {
            StereoKind::Center(
                StereogenicGeometry::new(geometry).expect("the geometry is stereogenic"),
            )
        }

        fn kinds() -> Vec<StereoKind> {
            STEREOGENIC
                .into_iter()
                .map(center)
                .chain([StereoKind::Bond, StereoKind::Axis])
                .collect()
        }

        fn neighbors(count: usize) -> Vec<SiteId> {
            (1..=count as u32).map(s).collect()
        }

        #[test]
        fn a_site_anchors_a_coordination_center() {
            for geometry in STEREOGENIC {
                assert!(
                    StereoLocus::Site(s(1)).anchors(center(geometry)),
                    "{geometry:?}"
                );
            }
        }

        #[test]
        fn a_bond_anchors_a_double_bond() {
            assert!(StereoLocus::Bond(b(1)).anchors(StereoKind::Bond));
        }

        #[test]
        fn an_axis_anchors_an_allene() {
            assert!(StereoLocus::Axis(s(1)).anchors(StereoKind::Axis));
        }

        #[test]
        fn only_a_site_anchors_a_coordination_center() {
            assert!(!StereoLocus::Bond(b(1)).anchors(center(Tetrahedral)));
            assert!(!StereoLocus::Axis(s(1)).anchors(center(Tetrahedral)));
        }

        #[test]
        fn only_a_bond_anchors_a_double_bond() {
            assert!(!StereoLocus::Site(s(1)).anchors(StereoKind::Bond));
            assert!(!StereoLocus::Axis(s(1)).anchors(StereoKind::Bond));
        }

        #[test]
        fn only_an_axis_anchors_an_allene() {
            assert!(!StereoLocus::Site(s(1)).anchors(StereoKind::Axis));
            assert!(!StereoLocus::Bond(b(1)).anchors(StereoKind::Axis));
        }

        #[test]
        fn loci_order_by_anchor_then_by_identifier() {
            assert!(StereoLocus::Site(s(9)) < StereoLocus::Bond(b(1)));
            assert!(StereoLocus::Bond(b(9)) < StereoLocus::Axis(s(1)));
            assert!(StereoLocus::Site(s(1)) < StereoLocus::Site(s(2)));
        }

        #[test]
        fn new_admits_every_geometry_with_more_than_one_configuration() {
            for geometry in STEREOGENIC {
                assert_eq!(
                    StereogenicGeometry::new(geometry).map(StereogenicGeometry::geometry),
                    Some(geometry),
                );
            }
        }

        #[test]
        fn new_rejects_a_geometry_with_a_single_configuration() {
            for geometry in DEGENERATE {
                assert_eq!(StereogenicGeometry::new(geometry), None, "{geometry:?}");
            }
        }

        #[test]
        fn stereogenic_geometries_order_as_their_geometries_do() {
            let mut ordered: Vec<StereogenicGeometry> = STEREOGENIC
                .into_iter()
                .filter_map(StereogenicGeometry::new)
                .collect();
            ordered.sort_unstable();
            assert!(ordered.is_sorted_by_key(|stereogenic| stereogenic.geometry()));
        }

        #[test]
        fn slot_count_is_the_neighbor_count_of_the_arrangement() {
            assert_eq!(center(Tetrahedral).slot_count(), 4);
            assert_eq!(center(Octahedral).slot_count(), 6);
            assert_eq!(StereoKind::Bond.slot_count(), 4);
            assert_eq!(StereoKind::Axis.slot_count(), 4);
        }

        #[test]
        fn a_centers_slot_count_is_its_geometrys() {
            for geometry in STEREOGENIC {
                assert_eq!(
                    center(geometry).slot_count(),
                    geometry.slot_count(),
                    "{geometry:?}"
                );
            }
        }

        #[test]
        fn configuration_count_is_the_number_of_distinct_stereoisomers() {
            assert_eq!(center(Tetrahedral).configuration_count(), 2);
            assert_eq!(center(TrigonalPrismatic).configuration_count(), 120);
            assert_eq!(StereoKind::Bond.configuration_count(), 2);
            assert_eq!(StereoKind::Axis.configuration_count(), 2);
        }

        #[test]
        fn a_centers_configuration_count_is_its_geometrys() {
            for geometry in STEREOGENIC {
                assert_eq!(
                    center(geometry).configuration_count(),
                    geometry.configuration_count(),
                    "{geometry:?}"
                );
            }
        }

        #[test]
        fn an_allene_is_chiral_and_a_double_bond_is_not() {
            assert!(StereoKind::Axis.is_chiral());
            assert!(!StereoKind::Bond.is_chiral());
        }

        #[test]
        fn a_centers_chirality_is_its_geometrys() {
            for geometry in STEREOGENIC {
                assert_eq!(
                    center(geometry).is_chiral(),
                    geometry.is_chiral(),
                    "{geometry:?}"
                );
            }
        }

        #[test]
        fn a_center_has_one_end_and_a_rigid_frame_has_two() {
            assert_eq!(center(Tetrahedral).ends(), 1);
            assert_eq!(StereoKind::Bond.ends(), 2);
            assert_eq!(StereoKind::Axis.ends(), 2);
        }

        #[test]
        fn a_centers_reference_directions_are_its_geometrys() {
            for geometry in STEREOGENIC {
                assert_eq!(
                    center(geometry).directions().len(),
                    geometry.slot_count(),
                    "{geometry:?}"
                );
            }
        }

        #[test]
        fn a_rigid_frame_has_no_reference_directions() {
            assert!(StereoKind::Bond.directions().is_empty());
            assert!(StereoKind::Axis.directions().is_empty());
        }

        #[test]
        fn kinds_order_by_geometry_then_bond_then_axis() {
            assert!(center(Tetrahedral) < center(Octahedral));
            assert!(center(TrigonalPrismatic) < StereoKind::Bond);
            assert!(StereoKind::Bond < StereoKind::Axis);
        }

        #[test]
        fn the_edge_group_permutes_its_slots() {
            for &element in EDGE {
                assert!(is_permutation(element, 4));
            }
            assert!(is_permutation(BOND.reflection, 4));
            assert!(is_permutation(AXIS.reflection, 4));
        }

        #[test]
        fn the_edge_group_contains_the_identity() {
            assert!(EDGE.contains(&identity(4).as_slice()));
        }

        #[test]
        fn the_edge_group_is_closed_under_composition() {
            for &g in EDGE {
                for &h in EDGE {
                    assert!(EDGE.contains(&compose(g, h).as_slice()));
                }
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
        fn new_rejects_a_kind_off_its_anchor() {
            assert!(
                StereoConfiguration::new(
                    StereoLocus::Bond(b(1)),
                    center(Tetrahedral),
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
                    center(Tetrahedral),
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
                    center(Tetrahedral),
                    neighbors(5)
                )
                .is_none()
            );
        }

        #[test]
        fn locus_returns_the_anchor() {
            let config = StereoConfiguration::new(
                StereoLocus::Site(s(7)),
                center(Tetrahedral),
                neighbors(4),
            )
            .unwrap();
            assert_eq!(config.locus(), StereoLocus::Site(s(7)));
        }

        #[test]
        fn kind_returns_the_arrangement() {
            let config = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
                neighbors(4),
            )
            .unwrap();
            assert_eq!(config.kind(), center(Tetrahedral));
        }

        #[test]
        fn neighbors_are_the_configuration_substituents() {
            let config = StereoConfiguration::new(
                StereoLocus::Site(s(9)),
                center(Tetrahedral),
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
                center(Tetrahedral),
                [s(1), s(2), s(3), s(4)],
            )
            .unwrap();
            let rotated = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
                [s(2), s(1), s(4), s(3)],
            )
            .unwrap();
            assert_eq!(reference, rotated);
        }

        #[test]
        fn a_reflection_of_the_neighbors_yields_a_different_configuration() {
            let reference = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
                [s(1), s(2), s(3), s(4)],
            )
            .unwrap();
            let mirror = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
                [s(2), s(1), s(3), s(4)],
            )
            .unwrap();
            assert_ne!(reference, mirror);
        }

        #[test]
        fn configurations_differ_when_their_locus_or_kind_or_neighbors_differ() {
            let base = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
                neighbors(4),
            )
            .unwrap();
            assert_ne!(
                base,
                StereoConfiguration::new(
                    StereoLocus::Site(s(2)),
                    center(Tetrahedral),
                    neighbors(4),
                )
                .unwrap(),
            );
            assert_ne!(
                base,
                StereoConfiguration::new(
                    StereoLocus::Site(s(1)),
                    center(SquarePlanar),
                    neighbors(4),
                )
                .unwrap(),
            );
            assert_ne!(
                base,
                StereoConfiguration::new(
                    StereoLocus::Site(s(1)),
                    center(Tetrahedral),
                    [s(1), s(2), s(3), s(5)],
                )
                .unwrap(),
            );
        }

        #[test]
        fn a_descriptor_carries_its_configurations_kind() {
            let config = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(SquarePlanar),
                neighbors(4),
            )
            .unwrap();
            assert_eq!(
                config.descriptor(|site| site.get() as usize).kind(),
                center(SquarePlanar),
            );
        }

        #[test]
        fn a_chiral_descriptor_differs_from_its_mirror() {
            let config = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
                neighbors(4),
            )
            .unwrap();
            let descriptor = config.descriptor(|site| site.get() as usize);
            assert_ne!(descriptor, descriptor.mirror());
        }

        #[test]
        fn an_achiral_descriptor_is_its_own_mirror() {
            let config =
                StereoConfiguration::new(StereoLocus::Bond(b(1)), StereoKind::Bond, neighbors(4))
                    .unwrap();
            let descriptor = config.descriptor(|site| site.get() as usize);
            assert_eq!(descriptor, descriptor.mirror());
        }

        #[test]
        fn a_descriptor_ignores_a_rotation_of_the_neighbors() {
            let rank = |site: SiteId| site.get() as usize;
            let reference = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
                [s(1), s(2), s(3), s(4)],
            )
            .unwrap();
            let rotated = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
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
                center(Tetrahedral),
                [s(1), s(2), s(3), s(4)],
            )
            .unwrap();
            let swapped = StereoConfiguration::new(
                StereoLocus::Site(s(1)),
                center(Tetrahedral),
                [s(1), s(3), s(2), s(4)],
            )
            .unwrap();
            assert_eq!(reference.descriptor(rank), swapped.descriptor(rank));
        }

        #[test]
        fn no_kind_holds_more_slots_than_the_buffer() {
            for kind in kinds() {
                assert!(kind.slot_count() <= MAX_SLOTS, "{kind:?}");
            }
        }
    }
}

pub use coordination_geometry::CoordinationGeometry;
pub use stereo::{
    StereoConfiguration, StereoDescriptor, StereoKind, StereoLocus, StereogenicGeometry,
};

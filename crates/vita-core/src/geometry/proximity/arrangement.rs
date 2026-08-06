use crate::tensor::{Point3, Vector3};
use crate::units::length::{Angstrom, Length, LengthUnit};
use crate::units::volume::CubicAngstrom;
use crate::{HasLattice, HasPositions, Lattice, Quantity, Scalar, SiteId};

/// The sites arranged by where they lie, so that what falls within a cutoff of a
/// place can be asked of them.
///
/// One arrangement answers every such question: which sites lie [`near`](Self::near) a
/// point, which [`pairs`](Self::pairs) of them lie within reach of one another, and —
/// as [`Neighbors`](super::Neighbors) — that same relation indexed by site. What it
/// reads itself comes with the displacement the cutoff was weighed against, so a
/// reading that scales by separation measures nothing twice; the relation it gathers
/// keeps the pairs alone.
///
/// Obtain via [`arrangement`] or [`periodic_arrangement`].
#[derive(Clone, Debug)]
pub struct Arrangement<V> {
    sites: Vec<SiteId>,
    places: Vec<Point3<V>>,
    keys: Vec<[u32; 3]>,
    starts: Vec<u32>,
    cutoff: V,
    frame: Frame<V>,
}

/// How an arrangement lays its cells and measures across them: from the least corner
/// the sites occupy, or around the lattice they repeat under.
#[derive(Clone, Debug)]
enum Frame<V> {
    Corner(Point3<V>),
    Torus(Torus<V>),
}

/// What a lattice adds to an arrangement: the cells it is divided into, the
/// translation each carry names, and the images no step can name.
#[derive(Clone, Debug)]
struct Torus<V> {
    lattice: Lattice<V>,
    dims: [u32; 3],
    carried: [Vector3<V>; 27],
    spares: Vec<Vector3<V>>,
}

impl<V: Scalar> Arrangement<V> {
    /// The cutoff the arrangement was laid for, in unit `U`.
    ///
    /// The reach of every question it answers: the cells are no narrower than it, so
    /// a wider one would outrun them.
    pub fn cutoff<U: LengthUnit>(&self) -> Length<V, U> {
        Length::<V, Angstrom>::new(self.cutoff).to::<U>()
    }

    /// Iterates the sites no further from `point` than the cutoff, each with the
    /// displacement from the point to it.
    ///
    /// A site standing at the point is among them. Yielded in the order the sites are
    /// arranged in, which is settled by the system and the cutoff but is not the
    /// ascending one [`Neighbors`](super::Neighbors) answers in.
    pub fn near<U: LengthUnit>(
        &self,
        point: Point3<Length<V, U>>,
    ) -> impl Iterator<Item = (SiteId, Vector3<Length<V, U>>)> + '_ {
        let point = self.settled(point.map(|value| value.to::<Angstrom>().value()));
        let cell = celled(&self.frame, self.cutoff, point);
        around()
            .filter_map(move |offset| self.stepped(cell, offset))
            .flat_map(move |(target, carry)| {
                let (from, to) = self.block(target);
                (from..to).map(move |slot| (slot, carry))
            })
            .filter_map(move |(slot, carry)| {
                let displacement = self.nearest(self.places[slot] - point, carry)?;
                Some((self.sites[slot], from_angstroms(displacement)))
            })
    }

    /// Iterates the pairs of sites no further apart than the cutoff, each once as
    /// `(a, b)` with `a` before `b`, and the displacement from `a` to `b`.
    ///
    /// The displacement is the one the cutoff was weighed against, so a reading that
    /// scales by separation measures nothing twice. Nothing is held between one pair
    /// and the next, so folding them all away costs no memory beyond the arrangement
    /// itself. Yielded in the order the sites are arranged in — settled by the system
    /// and the cutoff, but not the ascending one [`Neighbors`](super::Neighbors)
    /// answers in.
    pub fn pairs<U: LengthUnit>(
        &self,
    ) -> impl Iterator<Item = (SiteId, SiteId, Vector3<Length<V, U>>)> + '_ {
        self.couples().map(move |(near, far, displacement)| {
            let (near, far) = (self.sites[near], self.sites[far]);
            if near < far {
                (near, far, from_angstroms(displacement))
            } else {
                (far, near, from_angstroms(-displacement))
            }
        })
    }

    /// The sites in the order they are arranged in.
    pub(super) fn sites(&self) -> &[SiteId] {
        &self.sites
    }

    /// Iterates the places within reach of one another, each pair once, by the slot
    /// each holds in the arrangement, with the displacement from the near one to the
    /// far one.
    ///
    /// Every cell is put together with itself and with the cells a forward step
    /// reaches, so that no two cells meet twice.
    pub(super) fn couples(&self) -> impl Iterator<Item = (usize, usize, Vector3<V>)> + '_ {
        (0..self.keys.len())
            .flat_map(move |cell| {
                let (from, to) = self.span(cell);
                let mut sought = cell + 1;
                let alone = (from..to)
                    .flat_map(move |near| (near + 1..to).map(move |far| (near, far, [0; 3])));
                let across = forward()
                    .filter_map(move |offset| self.stepped(self.keys[cell], offset))
                    .flat_map(move |(target, carry)| {
                        let ordered = carry == [0; 3];
                        let found = self.sought(if ordered { sought } else { 0 }, target);
                        if ordered {
                            sought = match found {
                                Ok(index) | Err(index) => index,
                            };
                        }
                        let (start, end) = found.map_or((0, 0), |index| self.span(index));
                        (from..to)
                            .flat_map(move |near| (start..end).map(move |far| (near, far, carry)))
                    });
                alone.chain(across)
            })
            .filter_map(move |(near, far, carry)| {
                let displacement = self.nearest(self.places[far] - self.places[near], carry)?;
                Some((near, far, displacement))
            })
    }

    /// The representative of `point` the arrangement measures from.
    fn settled(&self, point: Point3<V>) -> Point3<V> {
        match &self.frame {
            Frame::Corner(_) => point,
            Frame::Torus(torus) => torus
                .lattice
                .wrap(point.map(Length::<V, Angstrom>::new))
                .map(Quantity::value),
        }
    }

    /// The slots the cell at `index` covers.
    fn span(&self, index: usize) -> (usize, usize) {
        (self.starts[index] as usize, self.starts[index + 1] as usize)
    }

    /// Where the cell `key` sits among the cells from `from` on — found, or where it
    /// would go.
    ///
    /// A step that carries nothing lands on the key it raises, so a run of such steps
    /// reaches cells in ascending order and each may resume from where the last
    /// settled rather than search the whole arrangement again. A step that carries
    /// lands wherever the wrap puts it, and no ordering of the cells can mend that: an
    /// order every step of a torus respected would be invariant under translation,
    /// which a finite group of cells cannot carry.
    fn sought(&self, from: usize, key: [u32; 3]) -> Result<usize, usize> {
        match self.keys[from..].binary_search(&key) {
            Ok(cell) => Ok(from + cell),
            Err(cell) => Err(from + cell),
        }
    }

    /// The slots the cell `key` holds, empty where no site falls in it.
    fn block(&self, key: [u32; 3]) -> (usize, usize) {
        self.sought(0, key).map_or((0, 0), |index| self.span(index))
    }

    /// The cell `offset` steps from `cell` and the whole cells that step carried it
    /// by, or `None` where the step leaves the arrangement.
    fn stepped(&self, cell: [u32; 3], offset: [i32; 3]) -> Option<([u32; 3], [i32; 3])> {
        match &self.frame {
            Frame::Corner(_) => shifted(cell, offset),
            Frame::Torus(torus) => wrapped(cell, offset, torus.dims),
        }
    }

    /// The shortest image of `displacement` no longer than the cutoff, or `None` where
    /// every image is longer.
    ///
    /// Weighed as squares, so that no root is ever taken. The step's own carry names
    /// the only image a divided axis can reach; an undivided one reaches further
    /// through the spares, and the shortest of all of them is the reading, so that a
    /// pair two images bring together is still measured by the nearer.
    fn nearest(&self, displacement: Vector3<V>, carry: [i32; 3]) -> Option<Vector3<V>> {
        let reach = self.cutoff * self.cutoff;
        let Frame::Torus(torus) = &self.frame else {
            return (displacement.norm_squared() <= reach).then_some(displacement);
        };
        let image = displacement + torus.carried[place(carry)];
        let mut shortest = (image.norm_squared(), image);
        for &spare in &torus.spares {
            let image = image + spare;
            let square = image.norm_squared();
            if square < shortest.0 {
                shortest = (square, image);
            }
        }
        (shortest.0 <= reach).then_some(shortest.1)
    }
}

/// The sites arranged so that what lies within `cutoff` of a place can be asked of
/// them.
///
/// A negative cutoff puts nothing within reach; a cutoff of zero leaves only what
/// coincides, at the cost of comparing every pair.
///
/// # Panics
///
/// Panics if any site of [`sites`](crate::HasSites::sites) is absent from
/// [`positions`](HasPositions::positions).
pub fn arrangement<S, V, U>(system: &S, cutoff: Length<V, U>) -> Arrangement<V>
where
    S: HasPositions<V>,
    V: Scalar,
    U: LengthUnit,
{
    let (sites, places) = placed(system);
    let corner = places.iter().copied().reduce(Point3::min);
    let cutoff = cutoff.to::<Angstrom>().value();
    gathered(
        sites,
        places,
        cutoff,
        Frame::Corner(corner.unwrap_or_default()),
    )
}

/// The sites arranged so that what lies within `cutoff` of a place can be asked of
/// them on the torus the lattice defines.
///
/// A negative cutoff puts nothing within reach; a cutoff of zero leaves only what
/// coincides, at the cost of comparing every pair. Every site is taken by its
/// representative inside the cell, so the coordinates need not lie in it.
///
/// # Panics
///
/// Panics if any site of [`sites`](crate::HasSites::sites) is absent from
/// [`positions`](HasPositions::positions).
pub fn periodic_arrangement<S, V, U>(system: &S, cutoff: Length<V, U>) -> Arrangement<V>
where
    S: HasPositions<V> + HasLattice<V>,
    V: Scalar,
    U: LengthUnit,
{
    let (sites, mut places) = placed(system);
    let cutoff = cutoff.to::<Angstrom>().value();
    let lattice = system.lattice();
    for place in &mut places {
        *place = lattice
            .wrap(place.map(Length::<V, Angstrom>::new))
            .map(Quantity::value);
    }
    let dims = divisions(lattice, cutoff);
    let frame = Frame::Torus(Torus {
        lattice,
        dims,
        carried: core::array::from_fn(|place| translation(lattice, steps(place))),
        spares: spared(lattice, dims, cutoff),
    });
    gathered(sites, places, cutoff, frame)
}

/// The sites and their positions in bare ångströms, in the order the system presents
/// them.
fn placed<S, V>(system: &S) -> (Vec<SiteId>, Vec<Point3<V>>)
where
    S: HasPositions<V>,
    V: Scalar,
{
    (
        system.sites().collect(),
        (system.positions::<Angstrom>())
            .map(|position| position.map(Quantity::value))
            .collect(),
    )
}

/// The sites and their places gathered into the cells they fall in, cell by cell.
///
/// The sites of a cell come to lie next to one another, so that a comparison against
/// that cell reads from adjacent memory rather than gathering from across the system,
/// and the arrangement holds one key per cell rather than one per site.
fn gathered<V: Scalar>(
    sites: Vec<SiteId>,
    places: Vec<Point3<V>>,
    cutoff: V,
    frame: Frame<V>,
) -> Arrangement<V> {
    let mut cells: Vec<([u32; 3], u32)> = if cutoff >= V::ZERO {
        (places.iter().enumerate())
            .map(|(slot, &place)| (celled(&frame, cutoff, place), slot as u32))
            .collect()
    } else {
        Vec::new()
    };
    cells.sort_unstable();
    let (mut keys, mut starts) = (Vec::new(), Vec::new());
    for (slot, &(key, _)) in cells.iter().enumerate() {
        if keys.last() != Some(&key) {
            keys.push(key);
            starts.push(slot as u32);
        }
    }
    starts.push(cells.len() as u32);
    Arrangement {
        sites: cells
            .iter()
            .map(|&(_, slot)| sites[slot as usize])
            .collect(),
        places: cells
            .iter()
            .map(|&(_, slot)| places[slot as usize])
            .collect(),
        keys,
        starts,
        cutoff,
        frame,
    }
}

/// The cell `place` falls in.
///
/// Cells are cubes of the cutoff laid from the corner, or the lattice divided `dims`
/// ways along each of its vectors, so that places within reach share a cell or lie in
/// adjacent ones. A cutoff too small to divide by leaves one cell, where the
/// comparison itself does the work.
fn celled<V: Scalar>(frame: &Frame<V>, cutoff: V, place: Point3<V>) -> [u32; 3] {
    match frame {
        Frame::Corner(corner) if cutoff > V::ZERO => ((place - *corner) / cutoff)
            .to_array()
            .map(|step| step.floor().to_f64() as u32),
        Frame::Corner(_) => [0; 3],
        Frame::Torus(torus) => {
            let fraction = torus
                .lattice
                .to_fractional(place.map(Length::<V, Angstrom>::new));
            [0, 1, 2].map(|axis| {
                let step = fraction[axis] * V::from_f64(f64::from(torus.dims[axis]));
                (step.to_f64() as u32).min(torus.dims[axis] - 1)
            })
        }
    }
}

/// The place `steps` holds among the twenty-seven a single step can take, with `-1`,
/// `0`, `1` taken to `0`, `1`, `2` along each axis.
fn place(steps: [i32; 3]) -> usize {
    ((steps[0] + 1) * 9 + (steps[1] + 1) * 3 + steps[2] + 1) as usize
}

/// The steps holding `place`, the inverse of [`place`].
fn steps(place: usize) -> [i32; 3] {
    [place / 9, place / 3 % 3, place % 3].map(|axis| axis as i32 - 1)
}

/// The twenty-seven steps around a cell, its own place among them.
fn around() -> impl Iterator<Item = [i32; 3]> {
    (0..27).map(steps)
}

/// The thirteen steps reaching forward, so that each pair of neighboring cells is
/// examined once; a cell's own sites are paired apart from these.
///
/// Each raises the key of the cell it steps from, so they reach cells in ascending
/// order.
fn forward() -> impl Iterator<Item = [i32; 3]> {
    (place([0; 3]) + 1..27).map(steps)
}

/// The cell's perpendicular width across each of its vectors, in ångströms: the
/// distance between the two faces the other two vectors span.
fn widths<V: Scalar>(lattice: Lattice<V>) -> [V; 3] {
    let (a, b, c) = (
        lattice.a::<Angstrom>().map(Quantity::value),
        lattice.b::<Angstrom>().map(Quantity::value),
        lattice.c::<Angstrom>().map(Quantity::value),
    );
    let volume = lattice.volume::<CubicAngstrom>().value();
    [b.cross(c), c.cross(a), a.cross(b)].map(|face| volume / face.norm())
}

/// The cells to lay along each lattice vector: as many as fit the cutoff across the
/// cell's perpendicular width, or one where fewer than three would fit.
///
/// Two cells to an axis would make each its own forward neighbor once wrapped, and so
/// examine a pair twice; one cell leaves that axis undivided instead. Three or more
/// leave every cell at least as wide as the cutoff and every other image of a site
/// further off than it, so the image a step names is the only one within reach.
fn divisions<V: Scalar>(lattice: Lattice<V>, cutoff: V) -> [u32; 3] {
    widths(lattice).map(|width| {
        if cutoff <= V::ZERO {
            return 1;
        }
        let count = (width / cutoff).to_f64() as u32;
        if count < 3 { 1 } else { count }
    })
}

/// The translation `steps` whole cells along each lattice vector, in bare ångströms.
fn translation<V: Scalar>(lattice: Lattice<V>, steps: [i32; 3]) -> Vector3<V> {
    let count = |axis: usize| V::from_f64(f64::from(steps[axis]));
    (lattice.a::<Angstrom>() * count(0)
        + lattice.b::<Angstrom>() * count(1)
        + lattice.c::<Angstrom>() * count(2))
    .map(Quantity::value)
}

/// The translations to the images no step around the arrangement can name.
///
/// A divided axis names its own image in the step it takes, and no other image of it
/// falls within reach; an undivided axis names none, so every translation along it
/// that could bring a pair together is spared for the comparison itself. Two places
/// inside the cell stand less than one width apart along an axis, so a translation of
/// `n` widths brings them together only while `|n|` falls short of one more than the
/// widths the cutoff spans. Empty where every axis is divided.
fn spared<V: Scalar>(lattice: Lattice<V>, dims: [u32; 3], cutoff: V) -> Vec<Vector3<V>> {
    let widths = widths(lattice);
    let bounds = [0, 1, 2].map(|axis| match dims[axis] {
        1 => (cutoff / widths[axis]).ceil().to_f64() as i32,
        _ => 0,
    });
    let mut spares = Vec::new();
    for a in -bounds[0]..=bounds[0] {
        for b in -bounds[1]..=bounds[1] {
            for c in -bounds[2]..=bounds[2] {
                if [a, b, c] != [0; 3] {
                    spares.push(translation(lattice, [a, b, c]));
                }
            }
        }
    }
    spares
}

/// The cell `offset` steps from `cell`, or `None` where that leaves the arrangement; a
/// corner has no boundary to cross, so nothing is carried.
fn shifted(cell: [u32; 3], offset: [i32; 3]) -> Option<([u32; 3], [i32; 3])> {
    let mut shifted = [0; 3];
    for axis in 0..3 {
        shifted[axis] = cell[axis].checked_add_signed(offset[axis])?;
    }
    Some((shifted, [0; 3]))
}

/// The cell `offset` steps from `cell` around the torus and the whole cells that step
/// carried it by, or `None` where the offset would step along an undivided axis and
/// revisit the cell it came from.
fn wrapped(cell: [u32; 3], offset: [i32; 3], dims: [u32; 3]) -> Option<([u32; 3], [i32; 3])> {
    let (mut wrapped, mut carry) = ([0; 3], [0; 3]);
    for axis in 0..3 {
        if dims[axis] == 1 {
            if offset[axis] != 0 {
                return None;
            }
            continue;
        }
        let span = i64::from(dims[axis]);
        let step = i64::from(cell[axis]) + i64::from(offset[axis]);
        wrapped[axis] = step.rem_euclid(span) as u32;
        carry[axis] = step.div_euclid(span) as i32;
    }
    Some((wrapped, carry))
}

/// A bare ångström displacement wrapped as lengths in unit `U`.
fn from_angstroms<V: Scalar, U: LengthUnit>(displacement: Vector3<V>) -> Vector3<Length<V, U>> {
    displacement.map(|value| Length::<V, Angstrom>::new(value).to::<U>())
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::fixture::{close, configuration, s};
    use crate::geometry::proximity::fixture::{chain, cube, reach, sheared, spot, strewn, sweep};
    use crate::units::length::Nanometer;

    fn paired(arrangement: &Arrangement<f64>) -> Vec<(SiteId, SiteId)> {
        let mut pairs: Vec<(SiteId, SiteId)> = (arrangement.pairs::<Angstrom>())
            .map(|(a, b, _)| (a, b))
            .collect();
        pairs.sort_unstable();
        pairs
    }

    fn nearby(arrangement: &Arrangement<f64>, point: [f64; 3]) -> Vec<SiteId> {
        let mut sites: Vec<SiteId> = arrangement
            .near(spot(point))
            .map(|(site, _)| site)
            .collect();
        sites.sort_unstable();
        sites
    }

    fn shortest(
        points: &[[f64; 3]],
        lattice: Option<Lattice<f64>>,
        a: SiteId,
        b: SiteId,
    ) -> Vector3<Length<f64, Angstrom>> {
        let place = |site: SiteId| Point3::from_array(points[site.get() as usize - 1]);
        let separation = (place(b) - place(a)).map(Length::new);
        lattice.map_or(separation, |lattice| lattice.minimum_image(separation))
    }

    #[test]
    fn an_empty_system_arranges_nothing() {
        let found = arrangement(&configuration(&[]), reach(1.0));
        assert!(paired(&found).is_empty() && nearby(&found, [0.0; 3]).is_empty());
    }

    #[test]
    fn a_lone_site_has_no_pair() {
        assert!(paired(&arrangement(&configuration(&[[0.0; 3]]), reach(1.0))).is_empty());
    }

    #[test]
    fn a_negative_cutoff_puts_nothing_within_reach() {
        let found = arrangement(&chain(), reach(-1.0));
        assert!(paired(&found).is_empty() && nearby(&found, [0.0; 3]).is_empty());
    }

    #[test]
    fn a_cutoff_of_zero_leaves_only_what_coincides() {
        let system = configuration(&[[0.0; 3], [0.0; 3], [1.0, 0.0, 0.0]]);
        assert_eq!(
            paired(&arrangement(&system, reach(0.0))),
            vec![(s(1), s(2))]
        );
    }

    #[test]
    fn pairs_are_the_sites_within_the_cutoff() {
        let found = arrangement(&chain(), reach(1.5));
        assert_eq!(
            paired(&found),
            vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4))]
        );
    }

    #[test]
    fn sites_exactly_at_the_cutoff_are_paired() {
        assert_eq!(paired(&arrangement(&chain(), reach(1.0))).len(), 3);
    }

    #[test]
    fn every_pair_gives_the_lesser_site_first() {
        let found = arrangement(&configuration(&strewn()), reach(3.0));
        assert!(found.pairs::<Angstrom>().all(|(a, b, _)| a < b));
    }

    #[test]
    fn a_pair_carries_the_displacement_from_the_lesser_site_to_the_greater() {
        let found = arrangement(&chain(), reach(1.5));
        let (_, _, displacement) = found.pairs::<Angstrom>().next().unwrap();
        assert!(close(displacement.x, Length::new(1.0)) && close(displacement.y, Length::ZERO));
    }

    #[test]
    fn the_cutoff_is_read_back_in_the_requested_unit() {
        let found = arrangement(&chain(), reach(1.5));
        assert!(close(found.cutoff::<Nanometer>(), Length::new(0.15)));
    }

    #[test]
    fn a_displacement_is_read_in_the_requested_unit() {
        let found = arrangement(&chain(), reach(1.5));
        let (_, _, displacement) = found.pairs::<Nanometer>().next().unwrap();
        assert!(close(displacement.x, Length::new(0.1)));
    }

    #[test]
    fn near_reads_the_point_in_the_requested_unit() {
        let found = arrangement(&chain(), reach(0.5));
        let point = Point3::from_array([0.1, 0.0, 0.0]).map(Length::<f64, Nanometer>::new);
        assert_eq!(
            found.near(point).map(|(site, _)| site).collect::<Vec<_>>(),
            vec![s(2)]
        );
    }

    #[test]
    fn near_lists_the_sites_within_the_cutoff_of_a_point() {
        let found = arrangement(&chain(), reach(1.5));
        assert_eq!(nearby(&found, [2.2, 0.0, 0.0]), vec![s(2), s(3), s(4)]);
    }

    #[test]
    fn a_site_standing_at_the_point_is_near_it() {
        assert_eq!(
            nearby(&arrangement(&chain(), reach(0.5)), [1.0, 0.0, 0.0]),
            vec![s(2)]
        );
    }

    #[test]
    fn near_carries_the_displacement_from_the_point() {
        let found = arrangement(&chain(), reach(0.5));
        let (_, displacement) = found.near(spot([0.7, 0.0, 0.0])).next().unwrap();
        assert!(close(displacement.x, Length::new(0.3)));
    }

    #[test]
    fn a_point_out_of_reach_of_every_site_has_nothing_near() {
        assert!(nearby(&arrangement(&chain(), reach(1.5)), [0.0, 9.0, 0.0]).is_empty());
    }

    #[test]
    fn a_point_before_the_corner_finds_what_lies_near_it() {
        assert_eq!(
            nearby(&arrangement(&chain(), reach(1.5)), [-1.0, 0.0, 0.0]),
            vec![s(1)]
        );
    }

    #[test]
    fn sites_in_diagonally_adjacent_cells_are_paired() {
        let system = configuration(&[[0.0; 3], [1.9, 1.9, 1.9], [2.1, 2.1, 2.1]]);
        assert_eq!(
            paired(&arrangement(&system, reach(2.0))),
            vec![(s(2), s(3))]
        );
    }

    #[test]
    fn pairs_are_independent_of_the_order_of_the_positions() {
        let shuffled = configuration(&[
            [3.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ]);
        let found = arrangement(&shuffled, reach(1.5));
        assert_eq!(
            paired(&found),
            vec![(s(1), s(4)), (s(2), s(3)), (s(2), s(4))]
        );
    }

    #[test]
    fn a_point_is_near_the_sites_across_the_boundary() {
        let found = periodic_arrangement(&cube(&[[0.5, 0.0, 0.0]]), reach(2.0));
        assert_eq!(nearby(&found, [9.8, 0.0, 0.0]), vec![s(1)]);
    }

    #[test]
    fn a_point_outside_the_cell_is_read_from_its_image_inside_it() {
        let found = periodic_arrangement(&cube(&[[0.5, 0.0, 0.0]]), reach(2.0));
        assert_eq!(nearby(&found, [20.3, 0.0, 0.0]), vec![s(1)]);
    }

    #[test]
    fn a_position_outside_the_cell_is_read_from_its_image_inside_it() {
        let system = cube(&[[3.0, 0.0, 0.0], [20.5, 0.0, 0.0]]);
        assert_eq!(
            paired(&periodic_arrangement(&system, reach(3.0))),
            vec![(s(1), s(2))]
        );
    }

    #[test]
    fn a_cutoff_of_zero_leaves_only_what_coincides_across_the_lattice() {
        let system = cube(&[[0.0; 3], [10.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let found = periodic_arrangement(&system, reach(0.0));
        assert_eq!(paired(&found), vec![(s(1), s(2))]);
    }

    #[test]
    fn an_axis_too_short_to_divide_pairs_across_itself_once() {
        let system = cube(&[[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]]);
        assert_eq!(
            paired(&periodic_arrangement(&system, reach(5.0))),
            vec![(s(1), s(2))]
        );
    }

    #[test]
    fn a_pair_two_images_bring_together_is_counted_once() {
        let system = cube(&[[0.0; 3], [5.0, 0.0, 0.0]]);
        assert_eq!(paired(&periodic_arrangement(&system, reach(5.0))).len(), 1);
    }

    #[test]
    fn a_pair_two_images_bring_together_is_measured_by_the_nearer() {
        let system = cube(&[[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]]);
        let found = periodic_arrangement(&system, reach(5.0));
        let (_, _, displacement) = found.pairs::<Angstrom>().next().unwrap();
        assert!(close(displacement.x, Length::new(-1.0)));
    }

    #[test]
    fn a_sheared_cell_is_divided_by_its_perpendicular_widths() {
        let system = sheared(&[[5.0, 5.0, 2.4], [5.0, 5.0, 5.3]]);
        assert_eq!(
            paired(&periodic_arrangement(&system, reach(3.0))),
            vec![(s(1), s(2))]
        );
    }

    #[test]
    fn a_negative_cutoff_puts_nothing_within_reach_across_the_lattice() {
        let system = cube(&[[0.5, 0.0, 0.0], [9.5, 0.0, 0.0]]);
        assert!(paired(&periodic_arrangement(&system, reach(-1.0))).is_empty());
    }

    #[test]
    fn the_arrangement_finds_what_a_full_sweep_would() {
        let points = strewn();
        let system = configuration(&points);
        let swept: Vec<(SiteId, SiteId)> = (1..=points.len() as u32)
            .flat_map(|a| (a + 1..=points.len() as u32).map(move |b| (a, b)))
            .filter(|&(a, b)| {
                let (near, far) = (points[a as usize - 1], points[b as usize - 1]);
                (0..3)
                    .map(|axis| (near[axis] - far[axis]).powi(2))
                    .sum::<f64>()
                    <= 9.0
            })
            .map(|(a, b)| (s(a), s(b)))
            .collect();
        let found = arrangement(&system, reach(3.0));
        assert_eq!(paired(&found), swept);
    }

    #[test]
    fn the_wrapped_arrangement_finds_what_a_full_sweep_would() {
        let system = cube(&strewn());
        let found = periodic_arrangement(&system, reach(3.0));
        assert_eq!(paired(&found), sweep(&system, reach(3.0)));
    }

    #[test]
    fn the_oblique_arrangement_finds_what_a_full_sweep_would() {
        let system = sheared(&strewn());
        let found = periodic_arrangement(&system, reach(3.0));
        assert_eq!(paired(&found), sweep(&system, reach(3.0)));
    }

    #[test]
    fn every_displacement_is_the_shortest_image_of_its_pair() {
        let points = strewn();
        let (cubic, oblique) = (cube(&points), sheared(&points));
        let frames = [
            (arrangement(&configuration(&points), reach(3.0)), None),
            (
                periodic_arrangement(&cubic, reach(3.0)),
                Some(cubic.lattice()),
            ),
            (
                periodic_arrangement(&oblique, reach(3.0)),
                Some(oblique.lattice()),
            ),
        ];
        for (found, lattice) in frames {
            assert!(found.pairs::<Angstrom>().all(|(a, b, displacement)| {
                let shortest = shortest(&points, lattice, a, b);
                (0..3).all(|axis| close(displacement[axis], shortest[axis]))
            }));
        }
    }

    #[test]
    fn near_finds_what_a_full_scan_would() {
        let points = strewn();
        let found = arrangement(&configuration(&points), reach(3.0));
        for probe in [
            [0.0; 3],
            [5.2, 4.7, 6.3],
            [9.9, 0.1, 4.4],
            [-2.0, 3.0, 11.0],
        ] {
            let scanned: Vec<SiteId> = (1..=points.len() as u32)
                .filter(|&index| {
                    let place = points[index as usize - 1];
                    (0..3)
                        .map(|axis| (place[axis] - probe[axis]).powi(2))
                        .sum::<f64>()
                        <= 9.0
                })
                .map(s)
                .collect();
            assert_eq!(nearby(&found, probe), scanned);
        }
    }
}

use std::cmp::Ordering;

use vita_core::tensor::{Point3, Vector3};
use vita_core::units::length::Angstrom;
use vita_core::{HasPositions, Quantity, Scalar, SiteId};

use super::{candidate_loci, frame};
use crate::capability::delegation::forward_capabilities;
use crate::{
    HasBondOrders, HasBonds, HasStereoConfigurations, StereoConfiguration, StereoKind, StereoLocus,
};

/// A set of stereo configurations over a molecule — one per stereogenic unit, ordered
/// by locus.
///
/// Obtain via [`perceive`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoConfigurations {
    configurations: Vec<StereoConfiguration>,
}

impl StereoConfigurations {
    /// Number of configurations.
    pub fn len(&self) -> usize {
        self.configurations.len()
    }

    /// Returns `true` if there are no configurations.
    pub fn is_empty(&self) -> bool {
        self.configurations.is_empty()
    }

    /// Iterates the configurations, ordered by locus.
    pub fn iter(&self) -> impl Iterator<Item = &StereoConfiguration> + '_ {
        self.configurations.iter()
    }

    /// The configuration at `locus`, or `None` if it bears none.
    pub fn get(&self, locus: StereoLocus) -> Option<&StereoConfiguration> {
        self.configurations
            .binary_search_by_key(&locus, StereoConfiguration::locus)
            .ok()
            .map(|index| &self.configurations[index])
    }

    /// The `configurations` gathered into a locus-ordered set.
    pub(super) fn from_configurations(mut configurations: Vec<StereoConfiguration>) -> Self {
        configurations.sort_unstable_by_key(StereoConfiguration::locus);
        StereoConfigurations { configurations }
    }

    /// Binds these configurations to `mol`, yielding a view that implements
    /// [`HasStereoConfigurations`].
    ///
    /// The view borrows both, so `mol` stays immutable while it is held — the stereo
    /// cannot silently fall out of step with the molecule it describes. Use it to feed
    /// a molecule's computed stereo to anything that reads the
    /// [`HasStereoConfigurations`] capability.
    pub fn bind<'a, M: HasBonds>(&'a self, mol: &'a M) -> WithStereoConfigurations<'a, M> {
        WithStereoConfigurations {
            mol,
            stereo_configurations: self,
        }
    }
}

/// A molecule viewed together with a set of [`StereoConfigurations`].
///
/// Answers the stereo configurations from that set and forwards every other core and
/// chem capability to the molecule, so a computed result reads as the
/// [`HasStereoConfigurations`] capability its consumers expect — at no cost beyond the
/// two references it holds.
///
/// Obtain via [`StereoConfigurations::bind`].
pub struct WithStereoConfigurations<'a, M> {
    mol: &'a M,
    stereo_configurations: &'a StereoConfigurations,
}

impl<M> Copy for WithStereoConfigurations<'_, M> {}

impl<M> Clone for WithStereoConfigurations<'_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

forward_capabilities!(
    WithStereoConfigurations,
    mol,
    HasAccelerations,
    HasElements,
    HasIsotopes,
    HasLattice,
    HasMasses,
    HasNetCharge,
    HasPositions,
    HasSites,
    HasVelocities,
    HasAromaticity,
    HasBondOrders,
    HasBonds,
    HasFormalCharges,
    HasPartialCharges,
    HasRadicalElectrons,
);

impl<M: HasBonds> HasStereoConfigurations for WithStereoConfigurations<'_, M> {
    fn stereo_configurations(&self) -> impl Iterator<Item = StereoConfiguration> + '_ {
        self.stereo_configurations.iter().cloned()
    }

    fn stereo_configuration_count(&self) -> usize {
        self.stereo_configurations.len()
    }

    fn stereo_configuration(&self, locus: StereoLocus) -> Option<StereoConfiguration> {
        self.stereo_configurations.get(locus).cloned()
    }
}

/// Perceives the stereo configurations a molecule's coordinates realize.
///
/// For each locus the caller's `candidate` admits, reads a [`StereoConfiguration`]
/// from the geometry: a center by matching the directions to its substituents onto
/// the geometry's reference frame — the pairwise angles fix which substituent fills
/// which slot, one signed volume the handedness — a double bond or allene by the
/// sign of a single invariant across its rigid double-bond chain. A locus whose
/// coordinates are degenerate — coplanar center substituents, an eclipsed double
/// bond — fixes no configuration and is skipped.
///
/// The perceived order follows [`StereoConfiguration`]'s convention, and that
/// convention is load-bearing: any other source — a SMILES `@`, a wedge — must be
/// read into the same one, or its configuration will disagree with what is perceived
/// here.
///
/// # Complexity
///
/// O(V + E) time and O(V + E) space, over the molecule's `V` sites and `E` bonds,
/// assuming [`bonds_of`](HasBonds::bonds_of) and
/// [`neighbors`](HasBonds::neighbors) run in O(degree) — bounded work per
/// candidate locus, its substituents capped at six.
pub fn perceive<M, V>(
    mol: &M,
    candidate: impl Fn(StereoLocus) -> Option<StereoKind>,
) -> StereoConfigurations
where
    M: HasBondOrders + HasPositions<V>,
    V: Scalar,
{
    let position = |site: SiteId| mol.position::<Angstrom>(site).map(|length| length.value());
    let position = &position;
    let candidate = &candidate;

    let perceive = move |locus: StereoLocus| -> Option<StereoConfiguration> {
        let kind = candidate(locus)?;
        let located = frame(mol, locus)?;
        let (anchors, subs) = (&located.anchors, &located.substituents);
        match locus {
            StereoLocus::Site(site) => central(site, kind, subs, position),
            StereoLocus::Bond(_) if kind == StereoKind::CisTrans => {
                cis_trans(locus, anchors, subs, position)
            }
            StereoLocus::Axis(_) if kind == StereoKind::Allene => {
                allene(locus, anchors, subs, position)
            }
            _ => None,
        }
    };

    StereoConfigurations::from_configurations(candidate_loci(mol).filter_map(perceive).collect())
}

/// The configuration of a central geometry from the directions to its substituents,
/// or `None` if the coordinates are degenerate or the count does not fit the geometry.
///
/// Assigns each substituent to the reference slot that reproduces the pairwise angles
/// it subtends — the assignment of least Gram residual — then, for a chiral geometry,
/// reflects the order if a signed volume shows the embedding to be improper.
fn central<V>(
    site: SiteId,
    kind: StereoKind,
    substituents: &[SiteId],
    position: &impl Fn(SiteId) -> Point3<V>,
) -> Option<StereoConfiguration>
where
    V: Scalar,
{
    let reference = kind.directions();
    let n = substituents.len();
    if n != reference.len() {
        return None;
    }
    let center = position(site);
    let observed: Vec<Vector3<V>> = substituents
        .iter()
        .map(|&sub| (position(sub) - center).try_normalize())
        .collect::<Option<_>>()?;
    let reference: Vec<Vector3<V>> = reference
        .iter()
        .map(|&[x, y, z]| {
            Vector3::from_array([V::from_f64(x), V::from_f64(y), V::from_f64(z)]).normalize()
        })
        .collect();

    let mut assignment: Vec<u8> = (0..n as u8).collect();
    let mut order: Option<(V, Vec<u8>)> = None;
    loop {
        let residual = gram_residual(&observed, &reference, &assignment);
        if order.as_ref().is_none_or(|(least, _)| residual < *least) {
            order = Some((residual, assignment.clone()));
        }
        if !next_permutation(&mut assignment) {
            break;
        }
    }
    let mut order = order?.1;

    if kind.is_chiral() {
        let triple = noncoplanar_triple(&reference)?;
        let placed: Vec<Vector3<V>> = order.iter().map(|&s| observed[s as usize]).collect();
        let handedness = orientation(&placed, triple);
        if handedness == V::ZERO {
            return None;
        }
        if handedness.signum() != orientation(&reference, triple).signum() {
            order = kind
                .reflection()
                .iter()
                .map(|&s| order[s as usize])
                .collect();
        }
    }

    StereoConfiguration::new(
        StereoLocus::Site(site),
        kind,
        order.iter().map(|&s| substituents[s as usize]),
    )
}

/// The summed absolute difference between the observed and reference Gram matrices
/// under a slot-to-substituent `assignment` — zero for an exact fit.
fn gram_residual<V: Scalar>(
    observed: &[Vector3<V>],
    reference: &[Vector3<V>],
    assignment: &[u8],
) -> V {
    let mut residual = V::ZERO;
    for (a, &row) in assignment.iter().enumerate() {
        for (b, &col) in assignment.iter().enumerate() {
            let seen = observed[row as usize].dot(observed[col as usize]);
            let want = reference[a].dot(reference[b]);
            residual += (seen - want).abs();
        }
    }
    residual
}

/// The signed volume the vectors at `triple` span — positive for a right-handed
/// triple, its sign flipped by a reflection.
fn orientation<V: Scalar>(vectors: &[Vector3<V>], triple: (usize, usize, usize)) -> V {
    let (a, b, c) = triple;
    vectors[a].dot(vectors[b].cross(vectors[c]))
}

/// The first slot triple the reference vectors span with nonzero volume — the frame
/// against which a chiral center's handedness is read.
fn noncoplanar_triple<V: Scalar>(reference: &[Vector3<V>]) -> Option<(usize, usize, usize)> {
    let n = reference.len();
    let threshold = V::from_f64(1e-3);
    (0..n)
        .flat_map(|a| (a + 1..n).flat_map(move |b| (b + 1..n).map(move |c| (a, b, c))))
        .find(|&triple| orientation(reference, triple).abs() > threshold)
}

/// The configuration of a double bond given its two termini (`anchors`) and their
/// substituents `[e1a, e1b, e2a, e2b]`, or `None` if the double bond is eclipsed.
fn cis_trans<V>(
    locus: StereoLocus,
    anchors: &[SiteId],
    subs: &[SiteId],
    position: &impl Fn(SiteId) -> Point3<V>,
) -> Option<StereoConfiguration>
where
    V: Scalar,
{
    let (first, second) = (anchors[0], anchors[1]);
    let n1 = (position(first) - position(subs[0])).cross(position(second) - position(first));
    let n2 = (position(second) - position(first)).cross(position(subs[2]) - position(second));
    let (near, far) = match n1.dot(n2).partial_cmp(&V::ZERO)? {
        Ordering::Greater => (subs[2], subs[3]),
        Ordering::Less => (subs[3], subs[2]),
        Ordering::Equal => return None,
    };
    StereoConfiguration::new(locus, StereoKind::CisTrans, [subs[0], subs[1], near, far])
}

/// The configuration of an allene given its two termini (`anchors`) and their
/// substituents `[e1a, e1b, e2a, e2b]`, or `None` if they are coplanar with the axis.
fn allene<V>(
    locus: StereoLocus,
    anchors: &[SiteId],
    subs: &[SiteId],
    position: &impl Fn(SiteId) -> Point3<V>,
) -> Option<StereoConfiguration>
where
    V: Scalar,
{
    let axis = position(anchors[1]) - position(anchors[0]);
    let twist = (position(subs[0]) - position(subs[1]))
        .cross(position(subs[2]) - position(subs[3]))
        .dot(axis);
    let ordered = match twist.partial_cmp(&V::ZERO)? {
        Ordering::Greater => [subs[0], subs[1], subs[2], subs[3]],
        Ordering::Less => [subs[1], subs[0], subs[2], subs[3]],
        Ordering::Equal => return None,
    };
    StereoConfiguration::new(locus, StereoKind::Allene, ordered)
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
    use vita_core::units::length::{Length, LengthUnit};

    use crate::BondOrder::{Double, Single};
    use crate::{BondId, BondOrder};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn only(target: StereoLocus, kind: StereoKind) -> impl Fn(StereoLocus) -> Option<StereoKind> {
        move |locus| (locus == target).then_some(kind)
    }

    fn config_at(locus: StereoLocus, kind: StereoKind, order: [u32; 4]) -> StereoConfiguration {
        StereoConfiguration::new(locus, kind, order.map(s)).unwrap()
    }

    fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    fn signed_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
    }

    struct Mol {
        sites: Vec<SiteId>,
        coords: Vec<[f64; 3]>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
    }

    impl Mol {
        fn at(&self, site: SiteId) -> [f64; 3] {
            self.coords[self.sites.iter().position(|&x| x == site).unwrap()]
        }
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasPositions<f64> for Mol {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<f64, U>> {
            let [x, y, z] = self.at(site);
            Point3::new(Length::new(x), Length::new(y), Length::new(z))
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

    fn mol(atoms: &[(u32, [f64; 3])], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, _)| s(id)).collect(),
            coords: atoms.iter().map(|&(_, xyz)| xyz).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, a, c, _)| (s(a), s(c))).collect(),
            orders: bonds.iter().map(|&(_, _, _, order)| order).collect(),
        }
    }

    fn mirrored(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.clone(),
            coords: m.coords.iter().map(|&[x, y, z]| [-x, y, z]).collect(),
            bonds: m.bonds.clone(),
            endpoints: m.endpoints.clone(),
            orders: m.orders.clone(),
        }
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            coords: m.coords.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
            orders: m.orders.iter().rev().copied().collect(),
        }
    }

    fn empty() -> Mol {
        mol(&[], &[])
    }

    fn tetrahedral() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 1.0, 1.0]),
                (3, [1.0, -1.0, -1.0]),
                (4, [-1.0, -1.0, 1.0]),
                (5, [-1.0, 1.0, -1.0]),
            ],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
            ],
        )
    }

    fn planar_center() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 0.0, 0.0]),
                (3, [0.0, 1.0, 0.0]),
                (4, [-1.0, 0.0, 0.0]),
                (5, [0.0, -1.0, 0.0]),
            ],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
            ],
        )
    }

    fn octahedral() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 0.0, 0.0]),
                (3, [-1.0, 0.0, 0.0]),
                (4, [0.0, 1.0, 0.0]),
                (5, [0.0, -1.0, 0.0]),
                (6, [0.0, 0.0, 1.0]),
                (7, [0.0, 0.0, -1.0]),
            ],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
                (5, 1, 6, Single),
                (6, 1, 7, Single),
            ],
        )
    }

    fn square_planar() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 0.0, 0.0]),
                (3, [0.0, 1.0, 0.0]),
                (4, [-1.0, 0.0, 0.0]),
                (5, [0.0, -1.0, 0.0]),
            ],
            &[
                (1, 1, 2, Single),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 1, 5, Single),
            ],
        )
    }

    fn cis_alkene() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 0.0, 0.0]),
                (3, [-0.5, 1.0, 0.0]),
                (4, [-0.5, -1.0, 0.0]),
                (5, [1.5, 1.0, 0.0]),
                (6, [1.5, -1.0, 0.0]),
            ],
            &[
                (1, 1, 2, Double),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 2, 5, Single),
                (5, 2, 6, Single),
            ],
        )
    }

    fn trans_alkene() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 0.0, 0.0]),
                (3, [-0.5, 1.0, 0.0]),
                (4, [-0.5, -1.0, 0.0]),
                (5, [1.5, -1.0, 0.0]),
                (6, [1.5, 1.0, 0.0]),
            ],
            &[
                (1, 1, 2, Double),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 2, 5, Single),
                (5, 2, 6, Single),
            ],
        )
    }

    fn degenerate_double_bond() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 0.0, 0.0]),
                (3, [-1.0, 0.0, 0.0]),
                (4, [-0.5, -1.0, 0.0]),
                (5, [1.5, 1.0, 0.0]),
                (6, [1.5, -1.0, 0.0]),
            ],
            &[
                (1, 1, 2, Double),
                (2, 1, 3, Single),
                (3, 1, 4, Single),
                (4, 2, 5, Single),
                (5, 2, 6, Single),
            ],
        )
    }

    fn allene() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 0.0, 0.0]),
                (3, [2.0, 0.0, 0.0]),
                (4, [-0.5, 1.0, 0.0]),
                (5, [-0.5, -1.0, 0.0]),
                (6, [2.5, 0.0, 1.0]),
                (7, [2.5, 0.0, -1.0]),
            ],
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

    fn planar_allene() -> Mol {
        mol(
            &[
                (1, [0.0, 0.0, 0.0]),
                (2, [1.0, 0.0, 0.0]),
                (3, [2.0, 0.0, 0.0]),
                (4, [-0.5, 1.0, 0.0]),
                (5, [-0.5, -1.0, 0.0]),
                (6, [2.5, 1.0, 0.0]),
                (7, [2.5, -1.0, 0.0]),
            ],
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
    fn empty_molecule_has_no_configurations() {
        let perceived = perceive(&empty(), |_| Some(StereoKind::Tetrahedral));
        assert_eq!(perceived.len(), 0);
        assert!(perceived.is_empty());
    }

    #[test]
    fn a_molecule_the_candidate_rejects_has_no_configurations() {
        assert!(perceive(&tetrahedral(), |_| None).is_empty());
    }

    #[test]
    fn a_tetrahedral_configuration_is_perceived_in_positive_orientation() {
        let molecule = tetrahedral();
        let perceived = perceive(
            &molecule,
            only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral),
        );
        let n = perceived.get(StereoLocus::Site(s(1))).unwrap().neighbors();
        let volume = signed_volume(
            difference(molecule.at(n[1]), molecule.at(n[0])),
            difference(molecule.at(n[2]), molecule.at(n[0])),
            difference(molecule.at(n[3]), molecule.at(n[0])),
        );
        assert!(volume > 0.0);
    }

    #[test]
    fn enantiomeric_coordinates_are_perceived_as_distinct_configurations() {
        let candidate = only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral);
        let right = perceive(&tetrahedral(), &candidate);
        let left = perceive(&mirrored(&tetrahedral()), &candidate);
        assert_ne!(
            right.get(StereoLocus::Site(s(1))),
            left.get(StereoLocus::Site(s(1))),
        );
    }

    #[test]
    fn a_coplanar_center_fixes_no_configuration() {
        let perceived = perceive(
            &planar_center(),
            only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral),
        );
        assert!(perceived.is_empty());
    }

    #[test]
    fn an_octahedral_configuration_is_perceived_in_positive_orientation() {
        let molecule = octahedral();
        let perceived = perceive(
            &molecule,
            only(StereoLocus::Site(s(1)), StereoKind::Octahedral),
        );
        let n = perceived.get(StereoLocus::Site(s(1))).unwrap().neighbors();
        let volume = signed_volume(
            difference(molecule.at(n[0]), molecule.at(n[1])),
            difference(molecule.at(n[2]), molecule.at(n[3])),
            difference(molecule.at(n[4]), molecule.at(n[5])),
        );
        assert!(volume > 0.0);
    }

    #[test]
    fn a_square_planar_configuration_is_perceived() {
        let perceived = perceive(
            &square_planar(),
            only(StereoLocus::Site(s(1)), StereoKind::SquarePlanar),
        );
        let config = perceived.get(StereoLocus::Site(s(1))).unwrap();
        assert_eq!(config.kind(), StereoKind::SquarePlanar);
        assert_eq!(config.neighbors().len(), 4);
    }

    #[test]
    fn a_cis_double_bond_is_perceived_in_reference_order() {
        let perceived = perceive(
            &cis_alkene(),
            only(StereoLocus::Bond(b(1)), StereoKind::CisTrans),
        );
        assert_eq!(
            perceived.get(StereoLocus::Bond(b(1))),
            Some(&config_at(
                StereoLocus::Bond(b(1)),
                StereoKind::CisTrans,
                [3, 4, 5, 6]
            )),
        );
    }

    #[test]
    fn a_trans_double_bond_reverses_the_far_end() {
        let perceived = perceive(
            &trans_alkene(),
            only(StereoLocus::Bond(b(1)), StereoKind::CisTrans),
        );
        assert_eq!(
            perceived.get(StereoLocus::Bond(b(1))),
            Some(&config_at(
                StereoLocus::Bond(b(1)),
                StereoKind::CisTrans,
                [3, 4, 6, 5]
            )),
        );
    }

    #[test]
    fn a_degenerate_double_bond_fixes_no_configuration() {
        let perceived = perceive(
            &degenerate_double_bond(),
            only(StereoLocus::Bond(b(1)), StereoKind::CisTrans),
        );
        assert!(perceived.is_empty());
    }

    #[test]
    fn an_allene_is_perceived_in_reference_order() {
        let perceived = perceive(&allene(), only(StereoLocus::Axis(s(2)), StereoKind::Allene));
        assert_eq!(
            perceived.get(StereoLocus::Axis(s(2))),
            Some(&config_at(
                StereoLocus::Axis(s(2)),
                StereoKind::Allene,
                [4, 5, 6, 7]
            )),
        );
    }

    #[test]
    fn the_opposite_allene_twist_swaps_the_first_terminus() {
        let perceived = perceive(
            &mirrored(&allene()),
            only(StereoLocus::Axis(s(2)), StereoKind::Allene),
        );
        assert_eq!(
            perceived.get(StereoLocus::Axis(s(2))),
            Some(&config_at(
                StereoLocus::Axis(s(2)),
                StereoKind::Allene,
                [5, 4, 6, 7]
            )),
        );
    }

    #[test]
    fn a_planar_allene_fixes_no_configuration() {
        let perceived = perceive(
            &planar_allene(),
            only(StereoLocus::Axis(s(2)), StereoKind::Allene),
        );
        assert!(perceived.is_empty());
    }

    #[test]
    fn count_reports_the_number_of_configurations() {
        let perceived = perceive(
            &tetrahedral(),
            only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral),
        );
        assert_eq!(perceived.len(), 1);
        assert!(!perceived.is_empty());
    }

    #[test]
    fn iter_yields_the_perceived_configurations() {
        let perceived = perceive(
            &tetrahedral(),
            only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral),
        );
        let configs: Vec<&StereoConfiguration> = perceived.iter().collect();
        assert_eq!(configs.len(), 1);
        assert_eq!(configs[0].locus(), StereoLocus::Site(s(1)));
    }

    #[test]
    fn get_returns_the_configuration_at_a_locus() {
        let perceived = perceive(
            &tetrahedral(),
            only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral),
        );
        let config = perceived.get(StereoLocus::Site(s(1))).unwrap();
        assert_eq!(config.kind(), StereoKind::Tetrahedral);
    }

    #[test]
    fn get_is_none_for_a_locus_without_a_configuration() {
        let perceived = perceive(
            &tetrahedral(),
            only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral),
        );
        assert!(perceived.get(StereoLocus::Site(s(99))).is_none());
    }

    #[test]
    fn bound_view_answers_the_stereo_configuration_capability() {
        let molecule = tetrahedral();
        let perceived = perceive(
            &molecule,
            only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral),
        );
        let view = perceived.bind(&molecule);
        assert_eq!(view.stereo_configuration_count(), 1);
        assert!(view.stereo_configuration(StereoLocus::Site(s(1))).is_some());
    }

    #[test]
    fn bound_view_forwards_the_skeleton() {
        let molecule = tetrahedral();
        let perceived = perceive(
            &molecule,
            only(StereoLocus::Site(s(1)), StereoKind::Tetrahedral),
        );
        let view = perceived.bind(&molecule);
        assert_eq!(view.bond_endpoints(b(1)), molecule.bond_endpoints(b(1)));
        assert_eq!(view.bond_count(), molecule.bond_count());
    }

    #[test]
    fn perception_is_independent_of_input_order() {
        let candidate = only(StereoLocus::Bond(b(1)), StereoKind::CisTrans);
        assert_eq!(
            perceive(&cis_alkene(), &candidate),
            perceive(&reversed(&cis_alkene()), &candidate),
        );
    }
}

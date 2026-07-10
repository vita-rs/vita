use std::cmp::Ordering;

use vita_core::tensor::{Point3, Vector3};
use vita_core::units::length::Angstrom;
use vita_core::{HasPositions, Scalar, SiteId};

use super::{frame, geometry, next_permutation};
use crate::capability::delegation::forward_capabilities;
use crate::{
    HasBondOrders, HasBonds, HasStereoConfigurations, StereoConfiguration, StereoKind, StereoLocus,
};

/// The stereo configurations perceived in a molecule — one per stereogenic unit its
/// coordinates resolve, ordered by locus.
///
/// Obtain via [`perceive`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoConfigurations {
    configurations: Vec<StereoConfiguration>,
}

impl StereoConfigurations {
    /// Number of perceived configurations.
    pub fn len(&self) -> usize {
        self.configurations.len()
    }

    /// Returns `true` if no configuration was perceived.
    pub fn is_empty(&self) -> bool {
        self.configurations.is_empty()
    }

    /// Iterates the perceived configurations, ordered by locus.
    pub fn iter(&self) -> impl Iterator<Item = &StereoConfiguration> + '_ {
        self.configurations.iter()
    }

    /// The configuration perceived at `locus`.
    ///
    /// Returns `None` if `locus` is absent or its coordinates fixed no configuration.
    pub fn get(&self, locus: StereoLocus) -> Option<&StereoConfiguration> {
        self.configurations
            .binary_search_by_key(&locus, StereoConfiguration::locus)
            .ok()
            .map(|index| &self.configurations[index])
    }

    /// Binds this perception to `mol`, yielding a view that implements
    /// [`HasStereoConfigurations`].
    ///
    /// The view borrows both, so `mol` stays immutable while it is held — the
    /// perception cannot silently fall out of step with the molecule it describes.
    /// Use it to feed a perceived molecule to anything that reads the
    /// [`HasStereoConfigurations`] capability.
    pub fn bind<'a, M: HasBonds>(&'a self, mol: &'a M) -> WithStereoConfigurations<'a, M> {
        WithStereoConfigurations {
            mol,
            stereo_configurations: self,
        }
    }
}

/// A molecule viewed together with its perceived [`StereoConfigurations`].
///
/// Answers the stereo configurations from the perception and forwards every other
/// core and chem capability to the molecule, so a computed result reads as the
/// [`HasStereoConfigurations`] capability its consumers expect — at no cost beyond
/// the two references it holds.
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
    HasHybridizations,
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

/// Perceives the stereo configurations a molecule's coordinates realise.
///
/// For each locus the caller's `candidate` admits, reads a [`StereoConfiguration`]
/// from the geometry: a centre by matching the directions to its substituents onto
/// the geometry's reference frame — the pairwise angles fix which substituent fills
/// which slot, one signed volume the handedness — a double bond or allene by the
/// sign of a single invariant across its rigid double-bond chain. A locus whose
/// coordinates are degenerate — coplanar centre substituents, an eclipsed double
/// bond — fixes no configuration and is skipped.
///
/// The perceived order follows [`StereoConfiguration`]'s convention, and that
/// convention is load-bearing: any other source — a SMILES `@`, a wedge — must be
/// read into the same one, or its handedness will disagree with what is perceived
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
    let perceive = &perceive;

    let sites = mol
        .sites()
        .filter_map(move |s| perceive(StereoLocus::Site(s)));
    let axes = mol
        .sites()
        .filter_map(move |s| perceive(StereoLocus::Axis(s)));
    let bonds = mol
        .bonds()
        .filter_map(move |b| perceive(StereoLocus::Bond(b)));

    let mut configurations: Vec<StereoConfiguration> = sites.chain(axes).chain(bonds).collect();
    configurations.sort_unstable_by_key(StereoConfiguration::locus);
    StereoConfigurations { configurations }
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
    let reference = geometry(kind).directions;
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
            order = geometry(kind)
                .reflection
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
/// against which a chiral centre's handedness is read.
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

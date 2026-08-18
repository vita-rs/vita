use vita_core::HasElements;

use super::CoordinationGeometries;
use crate::algorithm::valence::steric_numbers;
use crate::{CoordinationGeometry, HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// Perceives the coordination geometry of each of a molecule's sites.
///
/// A site's electron domains take the arrangement that separates them furthest — the
/// parent polyhedron of its [`steric_numbers`] count — and the substituents take
/// whichever of its vertices the lone pairs leave, the lone pairs claiming the roomiest
/// first. Naming what remains is VSEPR, and this is its whole content: a site is read
/// through its domain count and its coordination number, and nothing else about it
/// enters.
///
/// The model reaches the eleven geometries it names; the two a coordination sphere can
/// take without a domain count predicting them —
/// [`PyramidalizedSquare`](CoordinationGeometry::PyramidalizedSquare) and
/// [`TrigonalPrismatic`](CoordinationGeometry::TrigonalPrismatic) — never come out of it.
/// A site is left unanswered where the model says nothing: where [`steric_numbers`]
/// settles no count, where one substituent or none fixes no arrangement, or where the
/// domains outrun the vocabulary.
///
/// # Complexity
///
/// O((V + E) · log (V + E)) time and O(V + E) space, over the molecule's `V` sites and
/// `E` bonds, assuming [`degree`](crate::HasBonds::degree) runs in O(degree); counting the
/// domains dominates.
pub fn perceive<M>(mol: &M) -> CoordinationGeometries
where
    M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons,
{
    CoordinationGeometries::from_pairs(
        steric_numbers(mol)
            .iter()
            .filter_map(|(site, domains)| Some((site, vsepr(domains, mol.degree(site) as u32)?))),
    )
}

/// The geometry VSEPR gives a site holding `domains` electron domains of which
/// `coordination` are substituents, or `None` where the model names none.
///
/// Reading down a column is one parent arrangement losing a vertex to each further lone
/// pair: the trigonal bipyramid empties its equatorial belt before its poles, the
/// octahedron its trans pairs.
const fn vsepr(domains: u32, coordination: u32) -> Option<CoordinationGeometry> {
    match (domains, coordination) {
        (2, 2) => Some(CoordinationGeometry::Linear),
        (3, 3) => Some(CoordinationGeometry::TrigonalPlanar),
        (3, 2) => Some(CoordinationGeometry::Angular),
        (4, 4) => Some(CoordinationGeometry::Tetrahedral),
        (4, 3) => Some(CoordinationGeometry::TrigonalPyramidal),
        (4, 2) => Some(CoordinationGeometry::Angular),
        (5, 5) => Some(CoordinationGeometry::TrigonalBipyramidal),
        (5, 4) => Some(CoordinationGeometry::Seesaw),
        (5, 3) => Some(CoordinationGeometry::TShaped),
        (5, 2) => Some(CoordinationGeometry::Linear),
        (6, 6) => Some(CoordinationGeometry::Octahedral),
        (6, 5) => Some(CoordinationGeometry::SquarePyramidal),
        (6, 4) => Some(CoordinationGeometry::SquarePlanar),
        (6, 3) => Some(CoordinationGeometry::TShaped),
        (6, 2) => Some(CoordinationGeometry::Linear),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasSites, SiteId};

    use crate::BondOrder::Single;
    use crate::{BondId, BondOrder, HasBonds};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasElements for Mol {
        fn element(&self, site: SiteId) -> Element {
            self.elements[self.sites.iter().position(|&x| x == site).unwrap()]
        }
    }

    impl HasBonds for Mol {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
            self.endpoints[self.bonds.iter().position(|&x| x == bond).unwrap()]
        }
    }

    impl HasBondOrders for Mol {
        fn bond_order(&self, bond: BondId) -> BondOrder {
            self.orders[self.bonds.iter().position(|&x| x == bond).unwrap()]
        }
    }

    impl HasFormalCharges for Mol {
        fn formal_charge(&self, _site: SiteId) -> i8 {
            0
        }
    }

    impl HasRadicalElectrons for Mol {
        fn radical_electron(&self, _site: SiteId) -> u8 {
            0
        }
    }

    fn molecule(atoms: &[(u32, &str)], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, _)| s(id)).collect(),
            elements: atoms.iter().map(|&(_, symbol)| elem(symbol)).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, u, v, _)| (s(u), s(v))).collect(),
            orders: bonds.iter().map(|&(.., order)| order).collect(),
        }
    }

    fn centered(element: &str, ligands: u32) -> Mol {
        let mut atoms = vec![(1, element)];
        let mut bonds = Vec::new();
        for ligand in 1..=ligands {
            atoms.push((ligand + 1, "H"));
            bonds.push((ligand, 1, ligand + 1, Single));
        }
        molecule(&atoms, &bonds)
    }

    #[test]
    fn an_empty_molecule_has_no_geometries() {
        assert!(perceive(&molecule(&[], &[])).is_empty());
    }

    #[test]
    fn the_model_names_a_geometry_for_every_count_it_reaches() {
        let named = [
            ((2, 2), CoordinationGeometry::Linear),
            ((3, 3), CoordinationGeometry::TrigonalPlanar),
            ((3, 2), CoordinationGeometry::Angular),
            ((4, 4), CoordinationGeometry::Tetrahedral),
            ((4, 3), CoordinationGeometry::TrigonalPyramidal),
            ((4, 2), CoordinationGeometry::Angular),
            ((5, 5), CoordinationGeometry::TrigonalBipyramidal),
            ((5, 4), CoordinationGeometry::Seesaw),
            ((5, 3), CoordinationGeometry::TShaped),
            ((5, 2), CoordinationGeometry::Linear),
            ((6, 6), CoordinationGeometry::Octahedral),
            ((6, 5), CoordinationGeometry::SquarePyramidal),
            ((6, 4), CoordinationGeometry::SquarePlanar),
            ((6, 3), CoordinationGeometry::TShaped),
            ((6, 2), CoordinationGeometry::Linear),
        ];
        for ((domains, coordination), geometry) in named {
            assert_eq!(
                vsepr(domains, coordination),
                Some(geometry),
                "{domains} domains, {coordination} substituents"
            );
        }
    }

    #[test]
    fn methane_is_tetrahedral_about_its_carbon() {
        assert_eq!(
            perceive(&centered("C", 4)).get(s(1)),
            Some(CoordinationGeometry::Tetrahedral)
        );
    }

    #[test]
    fn a_lone_pair_bends_what_the_domains_would_have_left_straight() {
        assert_eq!(
            perceive(&centered("N", 3)).get(s(1)),
            Some(CoordinationGeometry::TrigonalPyramidal)
        );
        assert_eq!(
            perceive(&centered("O", 2)).get(s(1)),
            Some(CoordinationGeometry::Angular)
        );
    }

    #[test]
    fn a_hypervalent_site_takes_the_arrangement_its_domains_fill() {
        assert_eq!(
            perceive(&centered("S", 6)).get(s(1)),
            Some(CoordinationGeometry::Octahedral)
        );
        assert_eq!(
            perceive(&centered("P", 5)).get(s(1)),
            Some(CoordinationGeometry::TrigonalBipyramidal)
        );
    }

    #[test]
    fn the_model_names_none_outside_its_reach() {
        for domains in 0..=8u32 {
            for coordination in 0..=8u32 {
                let reached = (2..=6).contains(&domains)
                    && (2..=6).contains(&coordination)
                    && coordination <= domains;
                assert_eq!(
                    vsepr(domains, coordination).is_some(),
                    reached,
                    "{domains} domains, {coordination} substituents"
                );
            }
        }
    }

    #[test]
    fn a_site_bearing_one_substituent_is_unanswered() {
        let perceived = perceive(&centered("C", 4));
        assert_eq!(perceived.len(), 1);
        assert_eq!(perceived.get(s(2)), None);
    }

    #[test]
    fn a_site_whose_domains_go_uncounted_is_unanswered() {
        assert!(perceive(&centered("Fe", 6)).is_empty());
    }

    #[test]
    fn the_model_never_names_a_geometry_no_domain_count_predicts() {
        let reached: Vec<CoordinationGeometry> = (2..=6)
            .flat_map(|domains| (2..=6).filter_map(move |count| vsepr(domains, count)))
            .collect();
        assert!(!reached.contains(&CoordinationGeometry::PyramidalizedSquare));
        assert!(!reached.contains(&CoordinationGeometry::TrigonalPrismatic));
    }

    #[test]
    fn perception_is_independent_of_the_order_of_the_sites() {
        let forward = centered("N", 3);
        let backward = Mol {
            sites: forward.sites.iter().rev().copied().collect(),
            elements: forward.elements.iter().rev().copied().collect(),
            bonds: forward.bonds.iter().rev().copied().collect(),
            endpoints: forward.endpoints.iter().rev().copied().collect(),
            orders: forward.orders.iter().rev().copied().collect(),
        };
        assert_eq!(perceive(&forward), perceive(&backward));
    }
}

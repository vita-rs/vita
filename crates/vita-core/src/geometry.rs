pub mod measure;

#[cfg(test)]
mod fixture {
    use crate::tensor::Point3;
    use crate::units::length::{Angstrom, Length, LengthUnit};
    use crate::units::mass::{Dalton, Mass, MassUnit};
    use crate::{HasMasses, HasPositions, HasSites, Quantity, SiteId};

    pub struct System(Vec<([f64; 3], f64)>);

    impl HasSites for System {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            (1..=self.0.len() as u32).map(s)
        }
    }

    impl HasPositions<f64> for System {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<f64, U>> {
            Point3::from_array(self.0[site.get() as usize - 1].0)
                .map(|value| Length::<f64, Angstrom>::new(value).to::<U>())
        }
    }

    impl HasMasses<f64> for System {
        fn mass<U: MassUnit>(&self, site: SiteId) -> Mass<f64, U> {
            Mass::<f64, Dalton>::new(self.0[site.get() as usize - 1].1).to::<U>()
        }
    }

    pub fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    pub fn configuration(points: &[[f64; 3]]) -> System {
        System(points.iter().map(|&point| (point, 1.0)).collect())
    }

    pub fn weighted(points: &[([f64; 3], f64)]) -> System {
        System(points.to_vec())
    }

    pub fn close<Q: Quantity<Value = f64>>(a: Q, b: Q) -> bool {
        (a.value() - b.value()).abs() <= 1e-12
    }
}

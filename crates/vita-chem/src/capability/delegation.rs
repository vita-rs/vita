/// Implements the listed capabilities for a wrapper by delegating each method to
/// its inner.
macro_rules! forward_capabilities {
    ($wrapper:ident, $field:ident, $($capability:ident),+ $(,)?) => {
        $( $crate::capability::delegation::forward_capabilities!(@one $wrapper, $field, $capability); )+
    };

    (@one $wrapper:ident, $field:ident, HasAccelerations) => {
        impl<V: ::vita_core::Scalar, M: ::vita_core::HasAccelerations<V>>
            ::vita_core::HasAccelerations<V> for $wrapper<'_, M>
        {
            fn acceleration<U: ::vita_core::units::acceleration::AccelerationUnit>(
                &self,
                site: ::vita_core::SiteId,
            ) -> ::vita_core::tensor::Vector3<::vita_core::units::acceleration::Acceleration<V, U>> {
                self.$field.acceleration::<U>(site)
            }
            fn accelerations<U: ::vita_core::units::acceleration::AccelerationUnit>(
                &self,
            ) -> impl Iterator<Item = ::vita_core::tensor::Vector3<::vita_core::units::acceleration::Acceleration<V, U>>> + '_
            {
                self.$field.accelerations::<U>()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasElements) => {
        impl<M: ::vita_core::HasElements> ::vita_core::HasElements for $wrapper<'_, M> {
            fn element(&self, site: ::vita_core::SiteId) -> ::vita_core::Element {
                self.$field.element(site)
            }
            fn elements(&self) -> impl Iterator<Item = ::vita_core::Element> + '_ {
                self.$field.elements()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasIsotopes) => {
        impl<M: ::vita_core::HasIsotopes> ::vita_core::HasIsotopes for $wrapper<'_, M> {
            fn isotope(&self, site: ::vita_core::SiteId) -> ::vita_core::Isotope {
                self.$field.isotope(site)
            }
            fn isotopes(&self) -> impl Iterator<Item = ::vita_core::Isotope> + '_ {
                self.$field.isotopes()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasLattice) => {
        impl<V: ::vita_core::Scalar, M: ::vita_core::HasLattice<V>> ::vita_core::HasLattice<V>
            for $wrapper<'_, M>
        {
            fn lattice(&self) -> ::vita_core::Lattice<V> {
                self.$field.lattice()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasMasses) => {
        impl<V: ::vita_core::Scalar, M: ::vita_core::HasMasses<V>> ::vita_core::HasMasses<V>
            for $wrapper<'_, M>
        {
            fn mass<U: ::vita_core::units::mass::MassUnit>(
                &self,
                site: ::vita_core::SiteId,
            ) -> ::vita_core::units::mass::Mass<V, U> {
                self.$field.mass::<U>(site)
            }
            fn masses<U: ::vita_core::units::mass::MassUnit>(
                &self,
            ) -> impl Iterator<Item = ::vita_core::units::mass::Mass<V, U>> + '_ {
                self.$field.masses::<U>()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasNetCharge) => {
        impl<V: ::vita_core::Scalar, M: ::vita_core::HasNetCharge<V>> ::vita_core::HasNetCharge<V>
            for $wrapper<'_, M>
        {
            fn net_charge<U: ::vita_core::units::charge::ChargeUnit>(
                &self,
            ) -> ::vita_core::units::charge::Charge<V, U> {
                self.$field.net_charge::<U>()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasPositions) => {
        impl<V: ::vita_core::Scalar, M: ::vita_core::HasPositions<V>> ::vita_core::HasPositions<V>
            for $wrapper<'_, M>
        {
            fn position<U: ::vita_core::units::length::LengthUnit>(
                &self,
                site: ::vita_core::SiteId,
            ) -> ::vita_core::tensor::Point3<::vita_core::units::length::Length<V, U>> {
                self.$field.position::<U>(site)
            }
            fn positions<U: ::vita_core::units::length::LengthUnit>(
                &self,
            ) -> impl Iterator<Item = ::vita_core::tensor::Point3<::vita_core::units::length::Length<V, U>>> + '_
            {
                self.$field.positions::<U>()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasSites) => {
        impl<M: ::vita_core::HasSites> ::vita_core::HasSites for $wrapper<'_, M> {
            fn sites(&self) -> impl Iterator<Item = ::vita_core::SiteId> + '_ {
                self.$field.sites()
            }
            fn site_count(&self) -> usize {
                self.$field.site_count()
            }
            fn contains_site(&self, site: ::vita_core::SiteId) -> bool {
                self.$field.contains_site(site)
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasVelocities) => {
        impl<V: ::vita_core::Scalar, M: ::vita_core::HasVelocities<V>> ::vita_core::HasVelocities<V>
            for $wrapper<'_, M>
        {
            fn velocity<U: ::vita_core::units::velocity::VelocityUnit>(
                &self,
                site: ::vita_core::SiteId,
            ) -> ::vita_core::tensor::Vector3<::vita_core::units::velocity::Velocity<V, U>> {
                self.$field.velocity::<U>(site)
            }
            fn velocities<U: ::vita_core::units::velocity::VelocityUnit>(
                &self,
            ) -> impl Iterator<Item = ::vita_core::tensor::Vector3<::vita_core::units::velocity::Velocity<V, U>>> + '_
            {
                self.$field.velocities::<U>()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasAromaticity) => {
        impl<M: $crate::HasAromaticity> $crate::HasAromaticity for $wrapper<'_, M> {
            fn is_aromatic(&self, bond: $crate::BondId) -> bool {
                self.$field.is_aromatic(bond)
            }
            fn is_aromatic_site(&self, site: ::vita_core::SiteId) -> bool {
                self.$field.is_aromatic_site(site)
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasBondOrders) => {
        impl<M: $crate::HasBondOrders> $crate::HasBondOrders for $wrapper<'_, M> {
            fn bond_order(&self, bond: $crate::BondId) -> $crate::BondOrder {
                self.$field.bond_order(bond)
            }
            fn bond_orders(&self) -> impl Iterator<Item = $crate::BondOrder> + '_ {
                self.$field.bond_orders()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasBonds) => {
        impl<M: $crate::HasBonds> $crate::HasBonds for $wrapper<'_, M> {
            fn bonds(&self) -> impl Iterator<Item = $crate::BondId> + '_ {
                self.$field.bonds()
            }
            fn bond_endpoints(
                &self,
                bond: $crate::BondId,
            ) -> (::vita_core::SiteId, ::vita_core::SiteId) {
                self.$field.bond_endpoints(bond)
            }
            fn bond_count(&self) -> usize {
                self.$field.bond_count()
            }
            fn contains_bond(&self, bond: $crate::BondId) -> bool {
                self.$field.contains_bond(bond)
            }
            fn bond_between(
                &self,
                a: ::vita_core::SiteId,
                b: ::vita_core::SiteId,
            ) -> Option<$crate::BondId> {
                self.$field.bond_between(a, b)
            }
            fn bonds_of(
                &self,
                site: ::vita_core::SiteId,
            ) -> impl Iterator<Item = ($crate::BondId, ::vita_core::SiteId)> + '_ {
                self.$field.bonds_of(site)
            }
            fn neighbors(
                &self,
                site: ::vita_core::SiteId,
            ) -> impl Iterator<Item = ::vita_core::SiteId> + '_ {
                self.$field.neighbors(site)
            }
            fn degree(&self, site: ::vita_core::SiteId) -> usize {
                self.$field.degree(site)
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasFormalCharges) => {
        impl<M: $crate::HasFormalCharges> $crate::HasFormalCharges for $wrapper<'_, M> {
            fn formal_charge(&self, site: ::vita_core::SiteId) -> i8 {
                self.$field.formal_charge(site)
            }
            fn formal_charges(&self) -> impl Iterator<Item = i8> + '_ {
                self.$field.formal_charges()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasPartialCharges) => {
        impl<V: ::vita_core::Scalar, M: $crate::HasPartialCharges<V>> $crate::HasPartialCharges<V>
            for $wrapper<'_, M>
        {
            fn partial_charge<U: ::vita_core::units::charge::ChargeUnit>(
                &self,
                site: ::vita_core::SiteId,
            ) -> ::vita_core::units::charge::Charge<V, U> {
                self.$field.partial_charge::<U>(site)
            }
            fn partial_charges<U: ::vita_core::units::charge::ChargeUnit>(
                &self,
            ) -> impl Iterator<Item = ::vita_core::units::charge::Charge<V, U>> + '_ {
                self.$field.partial_charges::<U>()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasRadicalElectrons) => {
        impl<M: $crate::HasRadicalElectrons> $crate::HasRadicalElectrons for $wrapper<'_, M> {
            fn radical_electron(&self, site: ::vita_core::SiteId) -> u8 {
                self.$field.radical_electron(site)
            }
            fn radical_electrons(&self) -> impl Iterator<Item = u8> + '_ {
                self.$field.radical_electrons()
            }
        }
    };

    (@one $wrapper:ident, $field:ident, HasStereoConfigurations) => {
        impl<M: $crate::HasStereoConfigurations> $crate::HasStereoConfigurations for $wrapper<'_, M> {
            fn stereo_configurations(
                &self,
            ) -> impl Iterator<Item = $crate::StereoConfiguration> + '_ {
                self.$field.stereo_configurations()
            }
            fn stereo_configuration_count(&self) -> usize {
                self.$field.stereo_configuration_count()
            }
            fn stereo_configuration(
                &self,
                locus: $crate::StereoLocus,
            ) -> Option<$crate::StereoConfiguration> {
                self.$field.stereo_configuration(locus)
            }
        }
    };
}

pub(crate) use forward_capabilities;

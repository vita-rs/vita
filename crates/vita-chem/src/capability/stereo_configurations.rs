use crate::{HasBonds, StereoConfiguration, StereoLocus};

/// Per-unit stereochemistry: the [`StereoConfiguration`] of each stereogenic unit
/// whose arrangement the molecule fixes.
///
/// Access is by iterating
/// [`stereo_configurations`](HasStereoConfigurations::stereo_configurations), or by
/// keyed lookup:
/// [`stereo_configuration`](HasStereoConfigurations::stereo_configuration) maps a
/// [`StereoLocus`] to its configuration, if any.
///
/// # Contract
///
/// At most one configuration is declared per [`StereoLocus`]; a unit left undefined
/// is absent, not defaulted. Each configuration's neighbours are the substituents the
/// locus arranges, ordered as its [`StereoKind`](crate::StereoKind) documents. The
/// order the configurations are yielded in carries no meaning.
pub trait HasStereoConfigurations: HasBonds {
    /// Yields the declared configurations, in no particular order.
    fn stereo_configurations(&self) -> impl Iterator<Item = StereoConfiguration> + '_;

    /// Returns the number of declared configurations.
    ///
    /// The default implementation consumes
    /// [`stereo_configurations`](HasStereoConfigurations::stereo_configurations);
    /// override it when the count is known in `O(1)`.
    #[inline]
    fn stereo_configuration_count(&self) -> usize {
        self.stereo_configurations().count()
    }

    /// Returns the configuration declared at `locus`, or `None` if it bears none.
    ///
    /// The default implementation scans
    /// [`stereo_configurations`](HasStereoConfigurations::stereo_configurations) in
    /// `O(n)`; override it for `O(1)` lookup keyed on the locus.
    #[inline]
    fn stereo_configuration(&self, locus: StereoLocus) -> Option<StereoConfiguration> {
        self.stereo_configurations()
            .find(|config| config.locus() == locus)
    }
}

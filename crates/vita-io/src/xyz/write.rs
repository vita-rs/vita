use std::io::{self, Write};

use vita_core::units::length::LengthUnit;
use vita_core::{HasElements, HasPositions, Scalar};

use super::error::Error;

/// Configuration for writing a system as XYZ.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Config<'a> {
    /// The comment line. Must not contain a line terminator, which XYZ cannot represent.
    pub comment: &'a str,
}

/// Writes one system to `writer` in standard XYZ format.
///
/// Atoms are emitted in the order [`sites`](vita_core::HasSites::sites) yields them,
/// coordinates in unit `U`.
///
/// # Errors
///
/// Returns [`Error::Io`](crate::Error::Io) if `writer` fails or if
/// [`Config::comment`] contains a line terminator (`\n` or `\r`).
pub fn write<V, U, S>(writer: &mut impl Write, source: &S, config: &Config<'_>) -> Result<(), Error>
where
    V: Scalar,
    U: LengthUnit,
    S: HasElements + HasPositions<V>,
{
    if config.comment.contains(['\n', '\r']) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "an XYZ comment must not contain a line terminator",
        )
        .into());
    }

    writeln!(writer, "{}", source.site_count())?;
    writeln!(writer, "{}", config.comment)?;
    for site in source.sites() {
        let position = source.position::<U>(site);
        writeln!(
            writer,
            "{} {} {} {}",
            source.element(site).symbol(),
            position.x.value(),
            position.y.value(),
            position.z.value(),
        )?;
    }
    Ok(())
}

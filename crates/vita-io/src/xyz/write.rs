use std::io::{self, Write};

use vita_core::units::length::LengthUnit;
use vita_core::{HasElements, HasPositions, Scalar};

use super::Error;

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

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::tensor::Point3;
    use vita_core::units::length::{Angstrom, Length, Picometer};
    use vita_core::{Element, HasSites, SiteId};

    use crate::xyz;

    struct Atoms<V: Scalar = f64> {
        elements: Vec<Element>,
        positions: Vec<Point3<Length<V, Angstrom>>>,
    }

    impl<V: Scalar> HasSites for Atoms<V> {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            (1..=self.elements.len() as u32).map(|n| SiteId::new(n).unwrap())
        }
    }

    impl<V: Scalar> HasElements for Atoms<V> {
        fn element(&self, site: SiteId) -> Element {
            self.elements[(site.get() - 1) as usize]
        }
    }

    impl<V: Scalar> HasPositions<V> for Atoms<V> {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<V, U>> {
            self.positions[(site.get() - 1) as usize].map(|length| length.to())
        }
    }

    fn atoms(specs: &[(&str, f64, f64, f64)]) -> Atoms {
        Atoms {
            elements: specs
                .iter()
                .map(|&(symbol, ..)| Element::from_symbol(symbol).unwrap())
                .collect(),
            positions: specs
                .iter()
                .map(|&(_, x, y, z)| Point3::new(Length::new(x), Length::new(y), Length::new(z)))
                .collect(),
        }
    }

    fn written<U: LengthUnit>(source: &Atoms, comment: &str) -> Result<String, Error> {
        let mut buffer = Vec::new();
        write::<f64, U, _>(&mut buffer, source, &Config { comment })?;
        Ok(String::from_utf8(buffer).unwrap())
    }

    #[test]
    fn empty_system_writes_only_the_header_lines() {
        assert_eq!(
            written::<Angstrom>(&atoms(&[]), "note").unwrap(),
            "0\nnote\n"
        );
    }

    #[test]
    fn writes_the_count_comment_then_one_line_per_atom() {
        let source = atoms(&[("H", 0.0, 0.0, 0.0), ("O", 1.5, -2.0, 0.25)]);
        assert_eq!(
            written::<Angstrom>(&source, "water").unwrap(),
            "2\nwater\nH 0 0 0\nO 1.5 -2 0.25\n",
        );
    }

    #[test]
    fn writes_coordinates_in_the_requested_unit() {
        let source = atoms(&[("H", 1.0, 1.5, 0.25)]);
        assert_eq!(
            written::<Picometer>(&source, "c").unwrap(),
            "1\nc\nH 100 150 25\n",
        );
    }

    #[test]
    fn writes_an_f32_source() {
        let source = Atoms::<f32> {
            elements: vec![Element::from_symbol("H").unwrap()],
            positions: vec![Point3::new(
                Length::new(1.5),
                Length::new(0.0),
                Length::new(0.0),
            )],
        };
        let mut buffer = Vec::new();
        write::<f32, Angstrom, _>(&mut buffer, &source, &Config { comment: "c" }).unwrap();
        assert_eq!(String::from_utf8(buffer).unwrap(), "1\nc\nH 1.5 0 0\n");
    }

    #[test]
    fn writes_an_empty_comment_as_a_blank_second_line() {
        let source = atoms(&[("H", 0.0, 0.0, 0.0)]);
        assert_eq!(written::<Angstrom>(&source, "").unwrap(), "1\n\nH 0 0 0\n");
    }

    #[test]
    fn rejects_a_comment_containing_a_newline() {
        let error = written::<Angstrom>(&atoms(&[]), "a\nb").unwrap_err();
        assert!(matches!(error, Error::Io(e) if e.kind() == io::ErrorKind::InvalidInput));
    }

    #[test]
    fn rejects_a_comment_containing_a_carriage_return() {
        let error = written::<Angstrom>(&atoms(&[]), "a\rb").unwrap_err();
        assert!(matches!(error, Error::Io(e) if e.kind() == io::ErrorKind::InvalidInput));
    }

    #[test]
    fn propagates_a_writer_error() {
        struct Failing;
        impl Write for Failing {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "broken"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let source = atoms(&[("H", 0.0, 0.0, 0.0)]);
        let error =
            write::<f64, Angstrom, _>(&mut Failing, &source, &Config { comment: "c" }).unwrap_err();
        assert!(matches!(error, Error::Io(e) if e.kind() == io::ErrorKind::BrokenPipe));
    }

    #[test]
    fn a_system_survives_a_write_then_read_round_trip() {
        let original = xyz::read("2\nwater\nO 0 0 0\nH 0.757 0.586 0\n".as_bytes())
            .system::<f64>()
            .unwrap();
        let mut buffer = Vec::new();
        write::<f64, Angstrom, _>(
            &mut buffer,
            &original,
            &Config {
                comment: original.comment(),
            },
        )
        .unwrap();
        let roundtripped = xyz::read(buffer.as_slice()).system::<f64>().unwrap();
        assert_eq!(roundtripped, original);
    }
}

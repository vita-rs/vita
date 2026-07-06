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
    use crate::xyz;
    use vita_core::tensor::Point3;
    use vita_core::units::length::{Angstrom, Length, LengthUnit, Nanometer};
    use vita_core::{Element, HasSites, SiteId};

    fn text<V: Scalar, U: LengthUnit, S>(source: &S, comment: &str) -> String
    where
        S: HasElements + HasPositions<V>,
    {
        let mut out = Vec::new();
        write::<V, U, S>(&mut out, source, &Config { comment }).unwrap();
        String::from_utf8(out).unwrap()
    }

    struct Argon;

    impl HasSites for Argon {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            std::iter::once(SiteId::new(1).unwrap())
        }
        fn site_count(&self) -> usize {
            1
        }
    }

    impl HasElements for Argon {
        fn element(&self, _: SiteId) -> Element {
            Element::from_symbol("Ar").unwrap()
        }
    }

    impl HasPositions<f64> for Argon {
        fn position<U: LengthUnit>(&self, _: SiteId) -> Point3<Length<f64, U>> {
            Point3::new(
                Length::<f64, Angstrom>::new(1.0),
                Length::new(2.0),
                Length::new(3.0),
            )
            .map(|l| l.to())
        }
    }

    #[test]
    fn writes_frame() {
        let mol = xyz::read("2\ntest\nH 0 0 0\nHe 1 2 3\n".as_bytes())
            .system::<f64>()
            .unwrap();
        assert_eq!(
            text::<_, Angstrom, _>(&mol, "test"),
            "2\ntest\nH 0 0 0\nHe 1 2 3\n",
        );
    }

    #[test]
    fn writes_empty_comment() {
        let mol = xyz::read("1\nc\nH 0 0 0\n".as_bytes())
            .system::<f64>()
            .unwrap();
        assert_eq!(text::<_, Angstrom, _>(&mol, ""), "1\n\nH 0 0 0\n");
    }

    #[test]
    fn writes_zero_atoms() {
        let empty = xyz::read("0\nc\n".as_bytes()).system::<f64>().unwrap();
        assert_eq!(text::<_, Angstrom, _>(&empty, "c"), "0\nc\n");
    }

    #[test]
    fn writes_coordinates_in_specified_unit() {
        let mol = xyz::read("1\nc\nH 10 0 0\n".as_bytes())
            .system::<f64>()
            .unwrap();
        assert_eq!(text::<_, Nanometer, _>(&mol, "c"), "1\nc\nH 1 0 0\n");
    }

    #[test]
    fn writes_f32_coordinates() {
        let mol = xyz::read("1\nc\nH 1.5 0 0\n".as_bytes())
            .system::<f32>()
            .unwrap();
        assert_eq!(text::<_, Angstrom, _>(&mol, "c"), "1\nc\nH 1.5 0 0\n");
    }

    #[test]
    fn writes_any_capability_type() {
        assert_eq!(
            text::<_, Angstrom, _>(&Argon, "noble"),
            "1\nnoble\nAr 1 2 3\n"
        );
    }

    #[test]
    fn round_trips_through_read() {
        let input = "3\nwater\nO 0 0 0\nH 0.757 0.586 0\nH -0.757 0.586 0\n";
        let original = xyz::read(input.as_bytes()).system::<f64>().unwrap();
        let mut out = Vec::new();
        write::<_, Angstrom, _>(
            &mut out,
            &original,
            &Config {
                comment: original.comment(),
            },
        )
        .unwrap();
        assert_eq!(xyz::read(out.as_slice()).system::<f64>().unwrap(), original);
    }

    #[test]
    fn rejects_comment_with_lf() {
        let mut out = Vec::new();
        let err =
            write::<_, Angstrom, _>(&mut out, &Argon, &Config { comment: "a\nb" }).unwrap_err();
        match err {
            Error::Io(e) => assert_eq!(e.kind(), io::ErrorKind::InvalidInput),
            other => panic!("expected an I/O error, got {other:?}"),
        }
        assert!(out.is_empty());
    }

    #[test]
    fn rejects_comment_with_cr() {
        let mut out = Vec::new();
        let err =
            write::<_, Angstrom, _>(&mut out, &Argon, &Config { comment: "a\rb" }).unwrap_err();
        match err {
            Error::Io(e) => assert_eq!(e.kind(), io::ErrorKind::InvalidInput),
            other => panic!("expected an I/O error, got {other:?}"),
        }
        assert!(out.is_empty());
    }

    #[test]
    fn propagates_writer_error() {
        struct FailWriter;
        impl Write for FailWriter {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, ""))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let mol = xyz::read("1\nc\nH 0 0 0\n".as_bytes())
            .system::<f64>()
            .unwrap();
        let err =
            write::<_, Angstrom, _>(&mut FailWriter, &mol, &Config { comment: "c" }).unwrap_err();
        match err {
            Error::Io(e) => assert_eq!(e.kind(), io::ErrorKind::BrokenPipe),
            other => panic!("expected an I/O error, got {other:?}"),
        }
    }
}

/// The hybridization state of a bonded atom.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Hybridization {
    /// Spherical; one electron domain (1s orbital only).
    S,
    /// Linear; two electron domains.
    Sp,
    /// Trigonal planar; three electron domains.
    Sp2,
    /// Tetrahedral; four electron domains.
    Sp3,
    /// Square planar; four electron domains with d-orbital participation.
    Sp2d,
    /// Trigonal bipyramidal; five electron domains.
    Sp3d,
    /// Octahedral; six electron domains.
    Sp3d2,
    /// Pentagonal bipyramidal; seven electron domains.
    Sp3d3,
    /// Eight or more electron domains, f-block chemistry, or any case that
    /// the electron-domain model cannot classify.
    Other,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn copy_and_clone() {
        let a = Hybridization::Sp3;
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn eq() {
        assert_eq!(Hybridization::S, Hybridization::S);
        assert_ne!(Hybridization::S, Hybridization::Sp);
        assert_ne!(Hybridization::Sp, Hybridization::Sp2);
        assert_ne!(Hybridization::Sp2, Hybridization::Sp3);
        assert_ne!(Hybridization::Sp3, Hybridization::Sp2d);
    }

    #[test]
    fn debug() {
        assert_eq!(format!("{:?}", Hybridization::S), "S");
        assert_eq!(format!("{:?}", Hybridization::Sp), "Sp");
        assert_eq!(format!("{:?}", Hybridization::Sp2), "Sp2");
        assert_eq!(format!("{:?}", Hybridization::Sp3), "Sp3");
        assert_eq!(format!("{:?}", Hybridization::Sp2d), "Sp2d");
        assert_eq!(format!("{:?}", Hybridization::Sp3d), "Sp3d");
        assert_eq!(format!("{:?}", Hybridization::Sp3d2), "Sp3d2");
        assert_eq!(format!("{:?}", Hybridization::Sp3d3), "Sp3d3");
        assert_eq!(format!("{:?}", Hybridization::Other), "Other");
    }

    #[test]
    fn all_variants_distinct() {
        use std::collections::HashSet;
        let variants = [
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Sp3d3,
            Hybridization::Other,
        ];
        let set: HashSet<_> = variants.into_iter().collect();
        assert_eq!(set.len(), 9);
    }
}

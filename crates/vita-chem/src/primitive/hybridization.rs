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

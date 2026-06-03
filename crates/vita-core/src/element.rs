use core::fmt;
use core::num::NonZeroU8;

/// A chemical element, identified by its atomic number *Z* — the proton count.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Element(NonZeroU8);

/// Highest atomic number with an IUPAC-approved name; above it names are systematic.
const LAST_NAMED: u8 = 118;

impl Element {
    /// Constructs an element from atomic number `z`, returning `None` if `z` is zero.
    #[inline]
    pub const fn new(z: u8) -> Option<Self> {
        match NonZeroU8::new(z) {
            Some(z) => Some(Self(z)),
            None => None,
        }
    }

    /// Constructs an element from a [`NonZeroU8`] atomic number directly.
    #[inline]
    pub const fn from_nonzero(z: NonZeroU8) -> Self {
        Self(z)
    }

    /// Returns the atomic number *Z* (the proton count), always greater than zero.
    #[inline]
    pub const fn atomic_number(self) -> u8 {
        self.0.get()
    }

    /// Returns the element symbol (e.g. `"C"`, `"Og"`).
    ///
    /// Elements with `Z > 118` return their IUPAC systematic symbol (e.g. `"Uue"`).
    #[inline]
    pub fn symbol(self) -> &'static str {
        let z = self.0.get();
        if z <= LAST_NAMED {
            data::IUPAC_SYMBOLS[z as usize]
        } else {
            data::systematic_symbol(z)
        }
    }

    /// Returns the element name (e.g. `"Carbon"`, `"Oganesson"`).
    ///
    /// Elements with `Z > 118` return their IUPAC systematic name (e.g. `"Ununennium"`).
    #[inline]
    pub fn name(self) -> &'static str {
        let z = self.0.get();
        if z <= LAST_NAMED {
            data::IUPAC_NAMES[z as usize]
        } else {
            data::systematic_name(z)
        }
    }

    /// Returns the period (periodic-table row, `1..=7`), or `None` for `Z > 118`.
    #[inline]
    pub const fn period(self) -> Option<u8> {
        Some(match self.0.get() {
            1..=2 => 1,
            3..=10 => 2,
            11..=18 => 3,
            19..=36 => 4,
            37..=54 => 5,
            55..=86 => 6,
            87..=118 => 7,
            _ => return None,
        })
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

/// IUPAC element names and symbols, and the systematic-name generator.
mod data {
    /// IUPAC-approved symbols, indexed by atomic number; index `0` is an unused placeholder.
    #[rustfmt::skip]
    pub(super) const IUPAC_SYMBOLS: [&str; 119] = [
        "",
        "H",                                                                                                                                                                                      "He",
        "Li", "Be",                                                                                                                                                 "B",  "C",  "N",  "O",  "F",  "Ne",
        "Na", "Mg",                                                                                                                                                 "Al", "Si", "P",  "S",  "Cl", "Ar",
        "K",  "Ca",                                                                                     "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr",                                                                                     "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    ];

    /// IUPAC-approved names, indexed by atomic number; index `0` is an unused placeholder.
    #[rustfmt::skip]
    pub(super) const IUPAC_NAMES: [&str; 119] = [
        "",
        /*    1-10 */ "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon",
        /*   11-20 */ "Sodium", "Magnesium", "Aluminium", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium",
        /*   21-30 */ "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc",
        /*   31-40 */ "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium",
        /*   41-50 */ "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin",
        /*   51-60 */ "Antimony", "Tellurium", "Iodine", "Xenon", "Caesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
        /*   61-70 */ "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium",
        /*   71-80 */ "Lutetium", "Hafnium", "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinum", "Gold", "Mercury",
        /*   81-90 */ "Thallium", "Lead", "Bismuth", "Polonium", "Astatine", "Radon", "Francium", "Radium", "Actinium", "Thorium",
        /*  91-100 */ "Protactinium", "Uranium", "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium", "Einsteinium", "Fermium",
        /* 101-110 */ "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium", "Dubnium", "Seaborgium", "Bohrium", "Hassium", "Meitnerium", "Darmstadtium",
        /* 111-118 */ "Roentgenium", "Copernicium", "Nihonium", "Flerovium", "Moscovium", "Livermorium", "Tennessine", "Oganesson",
    ];

    /// Initial letter of each digit root, for systematic symbols.
    const ROOT_INITIAL: [u8; 10] = [b'n', b'u', b'b', b't', b'q', b'p', b'h', b's', b'o', b'e'];

    /// Digit roots for systematic names: `0..=9`.
    const ROOT: [&str; 10] = [
        "nil", "un", "bi", "tri", "quad", "pent", "hex", "sept", "oct", "enn",
    ];

    /// Lowest atomic number named systematically (one above [`super::LAST_NAMED`]).
    const SYSTEMATIC_LO: u16 = (super::LAST_NAMED as u16) + 1;
    /// Highest representable atomic number ([`u8::MAX`]).
    const SYSTEMATIC_HI: u16 = u8::MAX as u16;
    const SYSTEMATIC_LEN: usize = (SYSTEMATIC_HI - SYSTEMATIC_LO + 1) as usize;
    /// Upper bound on a systematic name length (three four-letter roots plus `"ium"`).
    const NAME_CAP: usize = 16;

    /// Compile-time table of IUPAC systematic names for `Z` in `SYSTEMATIC_LO..=SYSTEMATIC_HI`.
    struct Systematic {
        name: [[u8; NAME_CAP]; SYSTEMATIC_LEN],
        name_len: [u8; SYSTEMATIC_LEN],
        symbol: [[u8; 3]; SYSTEMATIC_LEN],
    }

    static SYSTEMATIC: Systematic = build();

    /// Returns the systematic symbol for `z`, sliced from the static table.
    pub(super) fn systematic_symbol(z: u8) -> &'static str {
        let idx = (z as u16 - SYSTEMATIC_LO) as usize;
        // SAFETY: `build` writes only ASCII digit-root initials.
        unsafe { core::str::from_utf8_unchecked(&SYSTEMATIC.symbol[idx]) }
    }

    /// Returns the systematic name for `z`, sliced from the static table.
    pub(super) fn systematic_name(z: u8) -> &'static str {
        let idx = (z as u16 - SYSTEMATIC_LO) as usize;
        let len = SYSTEMATIC.name_len[idx] as usize;
        // SAFETY: `build` writes only ASCII digit-root letters within `len`.
        unsafe { core::str::from_utf8_unchecked(&SYSTEMATIC.name[idx][..len]) }
    }

    /// Materializes the systematic table at compile time (IUPAC 1979 rules).
    const fn build() -> Systematic {
        let mut t = Systematic {
            name: [[0; NAME_CAP]; SYSTEMATIC_LEN],
            name_len: [0; SYSTEMATIC_LEN],
            symbol: [[0; 3]; SYSTEMATIC_LEN],
        };
        // Every `Z` in this range has exactly three digits.
        let mut z = SYSTEMATIC_LO;
        while z <= SYSTEMATIC_HI {
            let idx = (z - SYSTEMATIC_LO) as usize;
            let d0 = (z / 100) as usize;
            let d1 = (z / 10 % 10) as usize;
            let d2 = (z % 10) as usize;

            // Symbol: three digit-root initials, the first capitalized.
            t.symbol[idx][0] = ascii_upper(ROOT_INITIAL[d0]);
            t.symbol[idx][1] = ROOT_INITIAL[d1];
            t.symbol[idx][2] = ROOT_INITIAL[d2];

            // Name: roots concatenated with elisions, then the `-ium` suffix.
            let mut len = append(&mut t.name[idx], 0, ROOT[d0]);
            // `enn` (9) before `nil` (0) drops one `n`: write `"il"` for the `nil`.
            len = append(
                &mut t.name[idx],
                len,
                if d1 == 0 && d0 == 9 { "il" } else { ROOT[d1] },
            );
            len = append(
                &mut t.name[idx],
                len,
                if d2 == 0 && d1 == 9 { "il" } else { ROOT[d2] },
            );
            // A final `bi` (2) or `tri` (3) drops its trailing `i` before `-ium`.
            len = append(
                &mut t.name[idx],
                len,
                if d2 == 2 || d2 == 3 { "um" } else { "ium" },
            );

            t.name[idx][0] = ascii_upper(t.name[idx][0]);
            t.name_len[idx] = len as u8;
            z += 1;
        }
        t
    }

    /// Appends `s` to `buf` starting at `at`, returning the new length.
    const fn append(buf: &mut [u8; NAME_CAP], mut at: usize, s: &str) -> usize {
        let b = s.as_bytes();
        let mut i = 0;
        while i < b.len() {
            buf[at] = b[i];
            at += 1;
            i += 1;
        }
        at
    }

    /// Uppercases an ASCII lowercase byte; leaves other bytes unchanged.
    const fn ascii_upper(b: u8) -> u8 {
        if b.is_ascii_lowercase() { b - 32 } else { b }
    }
}

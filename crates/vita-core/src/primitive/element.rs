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

    /// Constructs an element from its symbol (e.g. `"C"`, `"Og"`, `"Uue"`),
    /// returning `None` if no element has that symbol.
    ///
    /// The exact, case-sensitive inverse of [`symbol`](Self::symbol).
    #[inline]
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        match symbol.len() {
            1 | 2 => data::named_from_symbol(symbol),
            3 => data::systematic_from_symbol(symbol),
            _ => None,
        }
        .and_then(Self::new)
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

    /// Returns the atomic number of the IUPAC-named element with this `symbol`, if any.
    pub(super) fn named_from_symbol(symbol: &str) -> Option<u8> {
        IUPAC_SYMBOLS[1..]
            .iter()
            .position(|&s| s == symbol)
            .map(|i| i as u8 + 1)
    }

    /// Initial letter of each digit root, for systematic symbols.
    const ROOT_INITIAL: [u8; 10] = *b"nubtqphsoe";

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

    /// Returns the atomic number for a three-letter systematic `symbol`, if it denotes
    /// an element in `SYSTEMATIC_LO..=SYSTEMATIC_HI`.
    pub(super) fn systematic_from_symbol(symbol: &str) -> Option<u8> {
        let &[c0, c1, c2] = symbol.as_bytes() else {
            return None;
        };
        if !c0.is_ascii_uppercase() {
            return None;
        }
        let d0 = digit_of_initial(c0.to_ascii_lowercase())? as u16;
        let d1 = digit_of_initial(c1)? as u16;
        let d2 = digit_of_initial(c2)? as u16;
        let z = 100 * d0 + 10 * d1 + d2;
        if (SYSTEMATIC_LO..=SYSTEMATIC_HI).contains(&z) {
            Some(z as u8)
        } else {
            None
        }
    }

    /// Returns the digit whose systematic root has the lowercase initial `b`.
    fn digit_of_initial(b: u8) -> Option<u8> {
        ROOT_INITIAL.iter().position(|&r| r == b).map(|d| d as u8)
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

#[cfg(test)]
mod tests {
    use super::*;

    fn element(z: u8) -> Element {
        Element::new(z).unwrap()
    }

    #[test]
    fn new_rejects_zero() {
        assert_eq!(Element::new(0), None);
    }

    #[test]
    fn atomic_number_returns_the_proton_count() {
        assert_eq!(element(6).atomic_number(), 6);
    }

    #[test]
    fn from_nonzero_preserves_the_atomic_number() {
        let z = NonZeroU8::new(6).unwrap();
        assert_eq!(Element::from_nonzero(z).atomic_number(), 6);
    }

    #[test]
    fn a_named_element_reports_its_iupac_symbol() {
        assert_eq!(element(1).symbol(), "H");
        assert_eq!(element(6).symbol(), "C");
        assert_eq!(element(118).symbol(), "Og");
    }

    #[test]
    fn a_named_element_reports_its_iupac_name() {
        assert_eq!(element(1).name(), "Hydrogen");
        assert_eq!(element(6).name(), "Carbon");
        assert_eq!(element(118).name(), "Oganesson");
    }

    #[test]
    fn from_symbol_finds_a_named_element() {
        assert_eq!(Element::from_symbol("C"), Some(element(6)));
        assert_eq!(Element::from_symbol("Og"), Some(element(118)));
    }

    #[test]
    fn period_maps_each_periodic_table_row() {
        assert_eq!(element(1).period(), Some(1));
        assert_eq!(element(3).period(), Some(2));
        assert_eq!(element(11).period(), Some(3));
        assert_eq!(element(19).period(), Some(4));
        assert_eq!(element(37).period(), Some(5));
        assert_eq!(element(55).period(), Some(6));
        assert_eq!(element(87).period(), Some(7));
    }

    #[test]
    fn from_symbol_rejects_unknown_symbols() {
        assert_eq!(Element::from_symbol("Xx"), None);
        assert_eq!(Element::from_symbol("Xyz"), None);
        assert_eq!(Element::from_symbol("Nnn"), None);
    }

    #[test]
    fn from_symbol_is_case_sensitive() {
        assert_eq!(Element::from_symbol("c"), None);
        assert_eq!(Element::from_symbol("he"), None);
    }

    #[test]
    fn from_symbol_rejects_the_wrong_length() {
        assert_eq!(Element::from_symbol(""), None);
        assert_eq!(Element::from_symbol("Abcd"), None);
    }

    #[test]
    fn elements_beyond_the_named_range_have_no_period() {
        assert_eq!(element(119).period(), None);
        assert_eq!(element(255).period(), None);
    }

    #[test]
    fn one_is_the_smallest_element() {
        assert_eq!(Element::new(1).unwrap().atomic_number(), 1);
    }

    #[test]
    fn the_largest_atomic_number_is_valid() {
        assert_eq!(Element::new(255).unwrap().atomic_number(), 255);
    }

    #[test]
    fn the_last_named_element_has_period_seven() {
        assert_eq!(element(118).period(), Some(7));
    }

    #[test]
    fn the_first_element_past_the_named_range_is_systematic() {
        assert_eq!(element(119).symbol(), "Uue");
        assert_eq!(element(119).name(), "Ununennium");
    }

    #[test]
    fn a_systematic_element_concatenates_its_digit_roots() {
        assert_eq!(element(120).symbol(), "Ubn");
        assert_eq!(element(120).name(), "Unbinilium");
    }

    #[test]
    fn systematic_name_elides_a_doubled_n_between_enn_and_nil() {
        assert_eq!(element(190).name(), "Unennilium");
    }

    #[test]
    fn systematic_name_drops_the_trailing_i_of_a_final_bi() {
        assert_eq!(element(122).name(), "Unbibium");
    }

    #[test]
    fn systematic_name_drops_the_trailing_i_of_a_final_tri() {
        assert_eq!(element(123).name(), "Unbitrium");
    }

    #[test]
    fn the_largest_element_is_named_systematically() {
        assert_eq!(element(255).symbol(), "Bpp");
        assert_eq!(element(255).name(), "Bipentpentium");
    }

    #[test]
    fn from_symbol_finds_a_systematic_element() {
        assert_eq!(Element::from_symbol("Uue"), Some(element(119)));
        assert_eq!(Element::from_symbol("Bpp"), Some(element(255)));
    }

    #[test]
    fn display_is_the_symbol() {
        assert_eq!(format!("{}", element(6)), "C");
        assert_eq!(format!("{}", element(119)), "Uue");
    }

    #[test]
    fn elements_order_by_atomic_number() {
        assert!(element(1) < element(2));
        assert!(element(8) > element(6));
    }

    #[test]
    fn elements_are_equal_when_their_atomic_numbers_match() {
        assert_eq!(element(6), element(6));
        assert_ne!(element(6), element(7));
    }

    #[test]
    fn from_symbol_inverts_symbol_for_every_element() {
        for z in 1..=u8::MAX {
            let e = element(z);
            assert_eq!(Element::from_symbol(e.symbol()), Some(e), "z = {z}");
        }
    }

    #[test]
    fn option_is_the_same_size_as_the_raw_integer() {
        assert_eq!(
            core::mem::size_of::<Option<Element>>(),
            core::mem::size_of::<u8>()
        );
    }
}

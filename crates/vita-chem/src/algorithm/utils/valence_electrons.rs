use vita_core::Element;

/// The valence-electron count of a main-group `element`: its outer-shell s and p
/// electrons, from one (group 1) to eight (group 18).
///
/// Returns `None` for the d- and f-block, whose valence the s- and p-electron
/// count alone does not fix.
pub fn valence_electrons(element: Element) -> Option<u8> {
    Some(match element.atomic_number() {
        1 | 3 | 11 | 19 | 37 | 55 | 87 => 1,
        2 | 4 | 12 | 20 | 38 | 56 | 88 => 2,
        5 | 13 | 31 | 49 | 81 | 113 => 3,
        6 | 14 | 32 | 50 | 82 | 114 => 4,
        7 | 15 | 33 | 51 | 83 | 115 => 5,
        8 | 16 | 34 | 52 | 84 | 116 => 6,
        9 | 17 | 35 | 53 | 85 | 117 => 7,
        10 | 18 | 36 | 54 | 86 | 118 => 8,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    #[test]
    fn counts_valence_electrons_by_group() {
        assert_eq!(valence_electrons(elem("H")), Some(1));
        assert_eq!(valence_electrons(elem("Be")), Some(2));
        assert_eq!(valence_electrons(elem("B")), Some(3));
        assert_eq!(valence_electrons(elem("C")), Some(4));
        assert_eq!(valence_electrons(elem("N")), Some(5));
        assert_eq!(valence_electrons(elem("O")), Some(6));
        assert_eq!(valence_electrons(elem("F")), Some(7));
        assert_eq!(valence_electrons(elem("Ne")), Some(8));
        assert_eq!(valence_electrons(elem("Pb")), Some(4));
        assert_eq!(valence_electrons(elem("Rn")), Some(8));
    }

    #[test]
    fn has_no_value_outside_the_main_group() {
        assert_eq!(valence_electrons(elem("Fe")), None);
        assert_eq!(valence_electrons(elem("Cu")), None);
        assert_eq!(valence_electrons(elem("La")), None);
        assert_eq!(valence_electrons(elem("U")), None);
    }

    #[test]
    fn every_value_is_between_one_and_eight() {
        for z in 1..=118 {
            if let Some(count) = valence_electrons(Element::new(z).unwrap()) {
                assert!((1..=8).contains(&count), "Z {z} has count {count}");
            }
        }
    }
}

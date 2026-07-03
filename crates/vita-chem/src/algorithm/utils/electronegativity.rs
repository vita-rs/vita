use vita_core::Element;

/// Pauling electronegativity of `element`: its pull on shared bonding electrons,
/// on the dimensionless scale from francium (0.70) to fluorine (3.98).
///
/// Returns `None` only where no single value is established.
pub fn electronegativity(element: Element) -> Option<f64> {
    Some(match element.atomic_number() {
        1 => 2.20,
        3 => 0.98,
        4 => 1.57,
        5 => 2.04,
        6 => 2.55,
        7 => 3.04,
        8 => 3.44,
        9 => 3.98,
        11 => 0.93,
        12 => 1.31,
        13 => 1.61,
        14 => 1.90,
        15 => 2.19,
        16 => 2.58,
        17 => 3.16,
        19 => 0.82,
        20 => 1.00,
        21 => 1.36,
        22 => 1.54,
        23 => 1.63,
        24 => 1.66,
        25 => 1.55,
        26 => 1.83,
        27 => 1.88,
        28 => 1.91,
        29 => 1.90,
        30 => 1.65,
        31 => 1.81,
        32 => 2.01,
        33 => 2.18,
        34 => 2.55,
        35 => 2.96,
        36 => 3.00,
        37 => 0.82,
        38 => 0.95,
        39 => 1.22,
        40 => 1.33,
        41 => 1.60,
        42 => 2.16,
        43 => 1.90,
        44 => 2.20,
        45 => 2.28,
        46 => 2.20,
        47 => 1.93,
        48 => 1.69,
        49 => 1.78,
        50 => 1.96,
        51 => 2.05,
        52 => 2.10,
        53 => 2.66,
        54 => 2.60,
        55 => 0.79,
        56 => 0.89,
        57 => 1.10,
        58 => 1.12,
        59 => 1.13,
        60 => 1.14,
        62 => 1.17,
        64 => 1.20,
        65 => 1.10,
        66 => 1.22,
        67 => 1.23,
        68 => 1.24,
        69 => 1.25,
        71 => 1.27,
        72 => 1.30,
        73 => 1.50,
        74 => 2.36,
        75 => 1.90,
        76 => 2.20,
        77 => 2.20,
        78 => 2.28,
        79 => 2.54,
        80 => 2.00,
        81 => 1.62,
        82 => 2.33,
        83 => 2.02,
        84 => 2.00,
        85 => 2.20,
        87 => 0.70,
        88 => 0.90,
        89 => 1.10,
        90 => 1.30,
        91 => 1.50,
        92 => 1.38,
        93 => 1.36,
        94 => 1.28,
        95 => 1.30,
        96 => 1.30,
        97 => 1.30,
        98 => 1.30,
        99 => 1.30,
        100 => 1.30,
        101 => 1.30,
        102 => 1.30,
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
    fn returns_the_pauling_value() {
        assert_eq!(electronegativity(elem("F")), Some(3.98));
        assert_eq!(electronegativity(elem("Fr")), Some(0.70));
        assert_eq!(electronegativity(elem("H")), Some(2.20));
        assert_eq!(electronegativity(elem("C")), Some(2.55));
        assert_eq!(electronegativity(elem("O")), Some(3.44));
        assert_eq!(electronegativity(elem("Cl")), Some(3.16));
        assert_eq!(electronegativity(elem("Na")), Some(0.93));
        assert_eq!(electronegativity(elem("Fe")), Some(1.83));
        assert_eq!(electronegativity(elem("Kr")), Some(3.00));
        assert_eq!(electronegativity(elem("La")), Some(1.10));
        assert_eq!(electronegativity(elem("U")), Some(1.38));
    }

    #[test]
    fn has_no_value_for_elements_without_an_established_number() {
        assert_eq!(electronegativity(elem("He")), None);
        assert_eq!(electronegativity(elem("Ar")), None);
        assert_eq!(electronegativity(elem("Rn")), None);
        assert_eq!(electronegativity(elem("Pm")), None);
        assert_eq!(electronegativity(elem("Eu")), None);
        assert_eq!(electronegativity(elem("Lr")), None);
    }

    #[test]
    fn every_value_lies_within_the_pauling_scale() {
        for z in 1..=118 {
            if let Some(value) = electronegativity(Element::new(z).unwrap()) {
                assert!(
                    (0.70..=3.98).contains(&value),
                    "Z {z} has value {value} outside the scale",
                );
            }
        }
    }
}

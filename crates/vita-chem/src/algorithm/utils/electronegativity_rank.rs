use vita_core::Element;

/// Electronegativity rank of `element`: its tier in the total order by which
/// the ionic approximation of oxidation states awards bond electrons, a
/// greater rank pulling harder.
///
/// The order follows the IUPAC 2016 oxidation-state recommendation.
/// Returns `None` beyond lawrencium, where no scale has measured values.
pub fn electronegativity_rank(element: Element) -> Option<u8> {
    Some(match element.atomic_number() {
        1 => 69,
        2 => 83,
        3 => 8,
        4 => 39,
        5 => 64,
        6 => 73,
        7 => 80,
        8 => 82,
        9 => 84,
        10 => 85,
        11 => 5,
        12 => 25,
        13 => 43,
        14 => 59,
        15 => 68,
        16 => 75,
        17 => 78,
        18 => 81,
        19 => 4,
        20 => 10,
        21 => 18,
        22 => 30,
        23 => 36,
        24 => 44,
        25 => 48,
        26 => 52,
        27 => 54,
        28 => 58,
        29 => 55,
        30 => 41,
        31 => 49,
        32 => 62,
        33 => 67,
        34 => 72,
        35 => 77,
        36 => 79,
        37 => 3,
        38 => 9,
        39 => 13,
        40 => 27,
        41 => 31,
        42 => 32,
        43 => 34,
        44 => 37,
        45 => 38,
        46 => 40,
        47 => 57,
        48 => 35,
        49 => 45,
        50 => 53,
        51 => 61,
        52 => 65,
        53 => 70,
        54 => 74,
        55 => 1,
        56 => 6,
        57 => 12,
        58 => 13,
        59 => 14,
        60 => 15,
        61 => 14,
        62 => 17,
        63 => 19,
        64 => 19,
        65 => 12,
        66 => 20,
        67 => 21,
        68 => 22,
        69 => 23,
        70 => 12,
        71 => 11,
        72 => 16,
        73 => 28,
        74 => 32,
        75 => 42,
        76 => 44,
        77 => 46,
        78 => 47,
        79 => 60,
        80 => 50,
        81 => 51,
        82 => 56,
        83 => 63,
        84 => 66,
        85 => 71,
        86 => 76,
        87 => 2,
        88 => 7,
        89 => 12,
        90 => 26,
        91 => 33,
        92 => 30,
        93 => 29,
        94 => 24,
        95 => 26,
        96 => 26,
        97 => 26,
        98 => 26,
        99 => 26,
        100 => 26,
        101 => 26,
        102 => 26,
        103 => 26,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    fn rank(symbol: &str) -> u8 {
        electronegativity_rank(elem(symbol)).unwrap()
    }

    #[test]
    fn ranks_every_element_through_lawrencium() {
        for z in 1..=103 {
            assert!(
                electronegativity_rank(Element::new(z).unwrap()).is_some(),
                "Z {z} lacks a rank",
            );
        }
    }

    #[test]
    fn has_no_rank_beyond_lawrencium() {
        for z in 104..=118 {
            assert!(
                electronegativity_rank(Element::new(z).unwrap()).is_none(),
                "Z {z} claims a rank",
            );
        }
    }

    #[test]
    fn the_tiers_fill_the_scale_without_gaps() {
        let ranks: Vec<u8> = (1..=103)
            .map(|z| electronegativity_rank(Element::new(z).unwrap()).unwrap())
            .collect();
        assert!(ranks.iter().all(|&tier| (1..=85).contains(&tier)));
        for tier in 1..=85 {
            assert!(ranks.contains(&tier), "no element holds rank {tier}");
        }
    }

    #[test]
    fn runs_from_caesium_to_neon() {
        assert_eq!(electronegativity_rank(elem("Cs")), Some(1));
        assert_eq!(electronegativity_rank(elem("Ne")), Some(85));
    }

    #[test]
    fn orders_the_organic_elements() {
        assert!(rank("O") > rank("N"));
        assert!(rank("N") > rank("S"));
        assert!(rank("S") > rank("C"));
        assert!(rank("C") > rank("H"));
    }

    #[test]
    fn follows_allen_where_the_scales_disagree() {
        assert!(rank("C") > rank("I"));
        assert!(rank("C") > rank("Se"));
        assert!(rank("C") > rank("Au"));
        assert!(rank("H") > rank("P"));
        assert!(rank("N") > rank("Cl"));
        assert!(rank("Fr") > rank("Cs"));
    }

    #[test]
    fn keeps_published_ties_tied() {
        assert_eq!(rank("Cr"), rank("Os"));
        assert_eq!(rank("Mo"), rank("W"));
        assert_eq!(rank("Ti"), rank("U"));
        assert_eq!(rank("Th"), rank("No"));
        assert_eq!(rank("No"), rank("Lr"));
    }
}

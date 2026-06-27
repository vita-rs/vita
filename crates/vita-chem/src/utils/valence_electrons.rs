use vita_core::Element;

/// Valence (outer-shell s and p) electron count of a main-group element.
///
/// Returns `None` for the d- and f-block, where the count is ambiguous.
pub(crate) fn valence_electrons(element: Element) -> Option<u8> {
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

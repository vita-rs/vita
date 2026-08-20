/// Advances `slice` to the next lexicographic permutation, `false` at the last.
///
/// Starting from a sorted slice, repeated calls therefore visit every permutation
/// exactly once and stop, without materializing them. The last permutation is left
/// as it is rather than wrapped back to the first.
///
/// # Complexity
///
/// O(n) time and O(1) space, over the slice's `n` elements; a sweep from the sorted
/// slice averages O(1) a call, so visiting every permutation costs O(n!) rather than
/// O(n · n!).
pub fn next_permutation(slice: &mut [u8]) -> bool {
    let n = slice.len();
    let Some(pivot) = (1..n).rev().find(|&i| slice[i - 1] < slice[i]) else {
        return false;
    };
    let successor = (pivot..n)
        .rev()
        .find(|&i| slice[i] > slice[pivot - 1])
        .expect("a successor exists past the pivot");
    slice.swap(pivot - 1, successor);
    slice[pivot..].reverse();
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn visited(start: &[u8]) -> Vec<Vec<u8>> {
        let mut slice = start.to_vec();
        let mut seen = vec![slice.clone()];
        while next_permutation(&mut slice) {
            seen.push(slice.clone());
        }
        seen
    }

    #[test]
    fn a_slice_too_short_to_reorder_does_not_advance() {
        assert!(!next_permutation(&mut []));
        assert!(!next_permutation(&mut [0]));
    }

    #[test]
    fn the_last_permutation_does_not_advance() {
        let mut descending = vec![2, 1, 0];
        assert!(!next_permutation(&mut descending));
        assert_eq!(descending, vec![2, 1, 0]);
    }

    #[test]
    fn repeated_elements_are_visited_once_each() {
        assert_eq!(
            visited(&[0, 0, 1]),
            vec![vec![0, 0, 1], vec![0, 1, 0], vec![1, 0, 0]]
        );
    }

    #[test]
    fn advancing_visits_every_permutation_once() {
        let seen = visited(&[0, 1, 2, 3]);
        let mut distinct = seen.clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(seen.len(), 24);
        assert_eq!(distinct.len(), 24);
    }

    #[test]
    fn the_permutations_come_out_in_lexicographic_order() {
        assert!(visited(&[0, 1, 2, 3]).is_sorted());
    }
}

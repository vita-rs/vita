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

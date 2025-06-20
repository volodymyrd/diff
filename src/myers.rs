use std::ops::{Index, IndexMut};

/// Represents a single operation in the diff script.
#[derive(Debug, PartialEq, Eq, Clone)]
pub enum DiffOperation {
    /// A character that is common to both sequences.
    Match(char),
    /// A character inserted into sequence A to match sequence B.
    Insertion(char),
    /// A character deleted from sequence A.
    Deletion(char),
}

/// Implements the Myers's O(ND) diff algorithm to find the shortest edit script
/// between two sequences of characters.
///
/// Returns a tuple: (minimum_edit_distance, diff_operations_list).
///
/// # Arguments
/// * `a` - The first sequence of characters.
/// * `b` - The second sequence of characters.
pub fn myers_diff(a: &[char], b: &[char]) -> (usize, Vec<DiffOperation>) {
    let (n, m) = (a.len(), b.len()); // Length of sequence A

    // --- Special Cases for efficiency and correctness ---
    // Both empty
    if n == 0 && m == 0 {
        return (0, Vec::new());
    }

    // Sequence A is empty, all operations are insertions
    if n == 0 {
        let path: Vec<DiffOperation> = b.iter().map(|&c| DiffOperation::Insertion(c)).collect();
        return (m, path);
    }

    // Sequence B is empty, all operations are deletions
    if m == 0 {
        let path: Vec<DiffOperation> = a.iter().map(|&c| DiffOperation::Deletion(c)).collect();
        return (n, path);
    }

    let mut path = vec![];
    let max_d = max_d(a.len(), b.len());
    let mut vf = V::new(max_d);
    let mut vr = V::new(max_d);
    let min_edit_dist = conquer(a, b, &mut vf, &mut vr, &mut path);

    (min_edit_dist, path)
}

fn conquer(a: &[char], b: &[char], vf: &mut V, vr: &mut V, path: &mut Vec<DiffOperation>) -> usize {
    let mut min_edit_dist = 0;

    // Check for common prefix
    let common_prefix_len = common_prefix_len(a, b);
    if common_prefix_len > 0 {
        for &c in &a[..common_prefix_len] {
            path.push(DiffOperation::Match(c));
        }
    }
    let a = &a[common_prefix_len..];
    let b = &b[common_prefix_len..];

    if a.is_empty() && b.is_empty() {
        return min_edit_dist;
    }

    // Check for common suffix, using tmp path to add it to the end if any
    let mut tmp_path = vec![];
    let common_suffix_len = common_suffix_len(a, b);
    if common_suffix_len > 0 {
        for &c in &a[a.len() - common_suffix_len..] {
            tmp_path.push(DiffOperation::Match(c));
        }
    }
    let a = &a[..a.len() - common_suffix_len];
    let b = &b[..b.len() - common_suffix_len];

    // Sequence A is empty, all operations are insertions
    if a.is_empty() {
        for &c in b {
            path.push(DiffOperation::Insertion(c));
        }
        if !tmp_path.is_empty() {
            path.append(&mut tmp_path);
        }
        return min_edit_dist + b.len();
    }

    // Sequence B is empty, all operations are deletions
    if b.is_empty() {
        for &c in a {
            path.push(DiffOperation::Deletion(c));
        }
        if !tmp_path.is_empty() {
            path.append(&mut tmp_path);
        }
        return min_edit_dist + a.len();
    }

    if let Some((x_start, y_start, x_end, y_end)) = find_middle_snake(a, b, vf, vr) {
        min_edit_dist += conquer(&a[..x_start], &b[..y_start], vf, vr, path);
        for &c in &a[x_start..x_end] {
            path.push(DiffOperation::Match(c));
        }
        min_edit_dist += conquer(&a[x_end..], &b[y_end..], vf, vr, path);
    } else {
        for &c in a {
            path.push(DiffOperation::Deletion(c));
        }
        for &c in b {
            path.push(DiffOperation::Insertion(c));
        }
        min_edit_dist += a.len() + b.len();
    }

    if !tmp_path.is_empty() {
        path.append(&mut tmp_path);
    }

    min_edit_dist
}

/// Calculates the length of the common prefix.
fn common_prefix_len(a: &[char], b: &[char]) -> usize {
    if a.is_empty() || b.is_empty() {
        return 0;
    }
    a.iter().zip(b.iter()).take_while(|x| x.0 == x.1).count()
}

/// Calculates the length of common suffix.
pub fn common_suffix_len(a: &[char], b: &[char]) -> usize {
    if a.is_empty() || b.is_empty() {
        return 0;
    }
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .take_while(|x| x.0 == x.1)
        .count()
}

// A D-path is a path which starts at (0,0) that has exactly D non-diagonal
// edges. All D-paths consist of a (D - 1)-path followed by a non-diagonal edge
// and then a possibly empty sequence of diagonal edges called a snake.

/// The `V` Array is the Furthest Reachable Point.
///
/// For each edit distance d, we don't need to store all paths. We only need to know,
/// for each k-diagonal, what is the furthest x-coordinate (i) we can reach using exactly
/// `d` edit operations.
/// `V[k]` store the maximum x-coordinate (row index i) reached on diagonal k after exactly
/// `d` edit operations. If we are at (i, j), then j is i - k.
/// So, `V[k]` tells us the coordinates of the furthest point reached on diagonal k is
/// (V[k], V[k] - k).
///
/// We can't use a traditional Vec to represent `V` since we use `k` as an index
/// that can take on negative values. So instead `V` is represented as a
/// light-weight wrapper around a Vec plus an `offset` which is the maximum value
/// `k` can take on in order to map negative `k`'s back to a value >= 0.
#[derive(Debug)]
struct V {
    offset: isize,
    v: Vec<usize>,
}

impl V {
    /// Create a new instance of `V`.
    /// `max_d` - the maximum possible edit distance, which is max(|A|, |B|).
    fn new(max_d: usize) -> Self {
        Self {
            offset: max_d as isize,
            v: vec![0; 2 * max_d + 1],
        }
    }
}

impl Index<isize> for V {
    type Output = usize;

    fn index(&self, index: isize) -> &Self::Output {
        &self.v[(index + self.offset) as usize]
    }
}

impl IndexMut<isize> for V {
    fn index_mut(&mut self, index: isize) -> &mut Self::Output {
        &mut self.v[(index + self.offset) as usize]
    }
}

fn max_d(n: usize, m: usize) -> usize {
    (n + m + 1) / 2
}

/// Finds the middle snake, returning (x_start, y_start, x_end, y_end) of the snake.
///
/// See the 4b. A Linear Space Refinement http://www.xmailserver.org/diff2.pdf.
fn find_middle_snake(
    a: &[char],
    b: &[char],
    vf: &mut V,
    vr: &mut V,
) -> Option<(usize, usize, usize, usize)> {
    let n = a.len();
    let m = b.len();

    if n == 0 && m == 0 {
        return None;
    }

    // ∆
    let delta = n as isize - m as isize;
    // 5 => 0x101 & 0x001 = 0x001 = 1 => true
    // 10 => 0x1010 & 0x0001 = 0x000 = 0 => false
    let is_delta_odd = (delta & 1) == 1;

    // Resetting V arrays for each find_middle_snake call.
    vf.v.fill(0);
    vr.v.fill(0);

    // We only need to explore the floor of the average [(a.len + b.len)/2].
    let d_max = max_d(n, m);

    for d in 0..d_max as isize {
        // Forward path
        for k in (-d..=d).step_by(2) {
            // Find the end of the furthest reaching forward D-path in diagonal k:
            // If k = −D or k ≠ D and V[k − 1] < V[k + 1] Then x ← V[k + 1] Else x ← V[k − 1]+1.
            let mut x = if k == -d || k != d && vf[k - 1] < vf[k + 1] {
                vf[k + 1]
            } else {
                vf[k - 1] + 1
            };
            // y ← x − k
            let mut y = (x as isize - k) as usize;

            // The coordinate of the start of a snake
            let (x0, y0) = (x, y);

            // While these sequences are identical, keep moving through the graph with no cost:
            // While x < N and y < M and a x + 1 = by + 1 Do (x,y) ← (x+1,y+1)
            if x < n && y < m {
                let common_prefix_len = common_prefix_len(&a[x..], &b[y..]);
                x += common_prefix_len;
                y += common_prefix_len;
            }

            // This is the new best x value: V[k] ← x
            vf[k] = x;

            // by Lemma 1, the optimal edit script length is odd or even as ∆ is odd or even.
            // Thus, when ∆ is odd, check for overlap only while extending forward paths
            // and when ∆ is even, check for overlap only while extending reverse paths.

            // If ∆ is odd and k ∈ [∆ − (D − 1) ,∆ + (D − 1)] Then
            //  If the path overlaps the furthest reaching reverse (D − 1)-path in diagonal k Then
            //    Length of an SES is 2D−1.
            //    The last snake of the forward path is the middle snake
            if is_delta_odd
                && k >= delta - (d - 1)
                && k <= delta + (d - 1)
                // reverse coordinates relate to forward coordinates by xr=n-x, yr=m-y,
                // kr=xr-yr=n-x-m+y=n-m-(x-y)=∆-k
                && vf[k] + vr[delta-k] >= n
            {
                return Some((x0, y0, x, y));
            }
        }

        // Reverse path (Simulating Reverse as Forward)
        for k in (-d..=d).step_by(2) {
            // Find the end of the furthest reaching reverse D-path in diagonal k+∆.
            let mut x = if k == -d || k != d && vr[k - 1] < vr[k + 1] {
                vr[k + 1]
            } else {
                vr[k - 1] + 1
            };

            // y ← x - kr
            let mut y = (x as isize - k) as usize;

            let (x_end, y_end) = (n - x, m - y);

            if x < n && y < m {
                let suffix_len = common_suffix_len(&a[..n - x], &b[..m - y]);
                x += suffix_len;
                y += suffix_len;
            }

            // This is the new best x value: V[kr] ← x
            vr[k] = x;

            // If ∆ is even and kr=k + ∆ ∈ [−D,D] Then
            //  If the path overlaps the furthest reaching forward D-path in diagonal k+∆ Then
            //   Length of an SES is 2D
            //   The last snake of the reverse path is the middle snake.
            if !is_delta_odd && k >= -d && k <= d && vr[k] + vf[delta - k] >= n {
                return Some((n - x, m - y, x_end, y_end));
            }
        }
    }

    None
}

// --- Unit Tests ---
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_myers_diff_empty_a() {
        let a = string_to_char_vec("");
        let b = string_to_char_vec("abc");
        let (distance, path) = myers_diff(&a, &b);
        assert_eq!(distance, 3);
        assert_eq!(path.len(), 3);
        assert_eq!(path[0], DiffOperation::Insertion('a'));
        assert_eq!(path[1], DiffOperation::Insertion('b'));
        assert_eq!(path[2], DiffOperation::Insertion('c'));
    }

    #[test]
    fn test_myers_diff_empty_b() {
        let a = string_to_char_vec("abc");
        let b = string_to_char_vec("");
        let (distance, path) = myers_diff(&a, &b);
        assert_eq!(distance, 3);
        assert_eq!(path.len(), 3);
        assert_eq!(path[0], DiffOperation::Deletion('a'));
        assert_eq!(path[1], DiffOperation::Deletion('b'));
        assert_eq!(path[2], DiffOperation::Deletion('c'));
    }

    #[test]
    fn test_myers_diff_common_prefix_identical() {
        let a = string_to_char_vec("abca");
        let b = string_to_char_vec("abca");
        let (distance, path) = myers_diff(&a, &b);
        assert_eq!(distance, 0);
        assert_eq!(path.len(), 4);
        assert_eq!(path[0], DiffOperation::Match('a'));
        assert_eq!(path[1], DiffOperation::Match('b'));
        assert_eq!(path[2], DiffOperation::Match('c'));
        assert_eq!(path[3], DiffOperation::Match('a'));
    }

    #[test]
    fn test_myers_diff_common_prefix_b_extends() {
        let a = string_to_char_vec("abca");
        let b = string_to_char_vec("abcaf"); // "f" is inserted
        let (distance, path) = myers_diff(&a, &b);
        assert_eq!(distance, 1);
        assert_eq!(path.len(), 5); // 4 matches + 1 insertion
        assert_eq!(path[0], DiffOperation::Match('a'));
        assert_eq!(path[1], DiffOperation::Match('b'));
        assert_eq!(path[2], DiffOperation::Match('c'));
        assert_eq!(path[3], DiffOperation::Match('a'));
        assert_eq!(path[4], DiffOperation::Insertion('f'));
    }

    #[test]
    fn test_myers_diff_common_prefix_a_extends() {
        let a = string_to_char_vec("abca"); // "a" is deleted
        let b = string_to_char_vec("abc");
        let (distance, path) = myers_diff(&a, &b);
        assert_eq!(distance, 1);
        assert_eq!(path.len(), 4); // 3 matches + 1 deletion
        assert_eq!(path[0], DiffOperation::Match('a'));
        assert_eq!(path[1], DiffOperation::Match('b'));
        assert_eq!(path[2], DiffOperation::Match('c'));
        assert_eq!(path[3], DiffOperation::Deletion('a'));
    }

    #[test]
    fn test_myers_diff_common_suffix_b_prepends_middle() {
        let a = string_to_char_vec("abca");
        let b = string_to_char_vec("fabca");
        let (distance, path) = myers_diff(&a, &b);
        assert_eq!(distance, 1);
        assert_eq!(path.len(), 5);
        assert_eq!(path[0], DiffOperation::Insertion('f'));
        assert_eq!(path[1], DiffOperation::Match('a'));
        assert_eq!(path[2], DiffOperation::Match('b'));
        assert_eq!(path[3], DiffOperation::Match('c'));
        assert_eq!(path[4], DiffOperation::Match('a'));
    }

    #[test]
    fn test_myers_diff_common_suffix_a_prepends_middle() {
        let a = string_to_char_vec("labca");
        let b = string_to_char_vec("abca");
        let (distance, path) = myers_diff(&a, &b);
        assert_eq!(distance, 1);
        assert_eq!(path.len(), 5);
        assert_eq!(path[0], DiffOperation::Deletion('l'));
        assert_eq!(path[1], DiffOperation::Match('a'));
        assert_eq!(path[2], DiffOperation::Match('b'));
        assert_eq!(path[3], DiffOperation::Match('c'));
        assert_eq!(path[4], DiffOperation::Match('a'));
    }

    #[test]
    fn test_myers_diff_common_suffix_middle_differs_stubbed() {
        let a = string_to_char_vec("labca");
        let b = string_to_char_vec("fabca");

        let (distance, path) = myers_diff(&a, &b);

        assert_eq!(distance, 2);
        assert_eq!(path.len(), 6);
        assert_eq!(path[0], DiffOperation::Deletion('l'));
        assert_eq!(path[1], DiffOperation::Insertion('f'));
        assert_eq!(path[2], DiffOperation::Match('a'));
        assert_eq!(path[3], DiffOperation::Match('b'));
        assert_eq!(path[4], DiffOperation::Match('c'));
        assert_eq!(path[5], DiffOperation::Match('a'));
    }

    #[test]
    fn test_myers_diff_example_1() {
        let a = string_to_char_vec("abcabba");
        let b = string_to_char_vec("cbabac");

        let (distance, path) = myers_diff(&a, &b);

        assert_eq!(distance, 5);
        // Total ops: 5 edits + 4 matches = 9 operations.
        assert_eq!(path.len(), 9);
        assert_eq!(path[0], DiffOperation::Deletion('a'));
        assert_eq!(path[1], DiffOperation::Insertion('c'));
        assert_eq!(path[2], DiffOperation::Match('b'));
        assert_eq!(path[3], DiffOperation::Deletion('c'));
        assert_eq!(path[4], DiffOperation::Match('a'));
        assert_eq!(path[5], DiffOperation::Match('b'));
        assert_eq!(path[6], DiffOperation::Deletion('b'));
        assert_eq!(path[7], DiffOperation::Match('a'));
        assert_eq!(path[8], DiffOperation::Insertion('c'));
    }

    #[test]
    fn test_common_prefix_len_both_empty() {
        let a = string_to_char_vec("");
        let b = string_to_char_vec("");
        assert_eq!(common_prefix_len(&a, &b), 0);
    }

    #[test]
    fn test_common_prefix_len_first_empty() {
        let a = string_to_char_vec("");
        let b = string_to_char_vec("abc");
        assert_eq!(common_prefix_len(&a, &b), 0);
    }

    #[test]
    fn test_common_prefix_len_second_empty() {
        let a = string_to_char_vec("abc");
        let b = string_to_char_vec("");
        assert_eq!(common_prefix_len(&a, &b), 0);
    }

    #[test]
    fn test_common_prefix_len_no_common_prefix() {
        let a = string_to_char_vec("abc");
        let b = string_to_char_vec("def");
        assert_eq!(common_prefix_len(&a, &b), 0);
    }

    #[test]
    fn test_common_prefix_len_partial_prefix() {
        let a = string_to_char_vec("apple");
        let b = string_to_char_vec("apply");
        assert_eq!(common_prefix_len(&a, &b), 4); // "appl"
    }

    #[test]
    fn test_common_prefix_len_full_prefix_identical() {
        let a = string_to_char_vec("hello");
        let b = string_to_char_vec("hello");
        assert_eq!(common_prefix_len(&a, &b), 5);
    }

    #[test]
    fn test_common_prefix_len_first_is_prefix_of_second() {
        let a = string_to_char_vec("cat");
        let b = string_to_char_vec("caterpillar");
        assert_eq!(common_prefix_len(&a, &b), 3);
    }

    #[test]
    fn test_common_prefix_len_second_is_prefix_of_first() {
        let a = string_to_char_vec("banana");
        let b = string_to_char_vec("ban");
        assert_eq!(common_prefix_len(&a, &b), 3);
    }

    #[test]
    fn test_common_prefix_len_different_cases() {
        let a = string_to_char_vec("Apple");
        let b = string_to_char_vec("apple");
        assert_eq!(common_prefix_len(&a, &b), 0); // Case-sensitive
    }

    #[test]
    fn test_common_prefix_len_with_special_chars() {
        let a = string_to_char_vec("!@#123");
        let b = string_to_char_vec("!@#abc");
        assert_eq!(common_prefix_len(&a, &b), 3);
    }

    #[test]
    fn test_common_suffix_len_both_empty() {
        let a = string_to_char_vec("");
        let b = string_to_char_vec("");
        assert_eq!(common_suffix_len(&a, &b), 0);
    }

    #[test]
    fn test_common_suffix_len_first_empty() {
        let a = string_to_char_vec("");
        let b = string_to_char_vec("abc");
        assert_eq!(common_suffix_len(&a, &b), 0);
    }

    #[test]
    fn test_common_suffix_len_second_empty() {
        let a = string_to_char_vec("abc");
        let b = string_to_char_vec("");
        assert_eq!(common_suffix_len(&a, &b), 0);
    }

    #[test]
    fn test_common_suffix_len_no_common_suffix() {
        let a = string_to_char_vec("abc");
        let b = string_to_char_vec("def");
        assert_eq!(common_suffix_len(&a, &b), 0);
    }

    #[test]
    fn test_common_suffix_len_partial_suffix() {
        let a = string_to_char_vec("apple");
        let b = string_to_char_vec("grapple");
        assert_eq!(common_suffix_len(&a, &b), 5); // "apple"
    }

    #[test]
    fn test_common_suffix_len_partial_suffix_variant() {
        let a = string_to_char_vec("testing");
        let b = string_to_char_vec("flying");
        assert_eq!(common_suffix_len(&a, &b), 3); // "ing"
    }

    #[test]
    fn test_common_suffix_len_full_suffix_identical() {
        let a = string_to_char_vec("hello");
        let b = string_to_char_vec("hello");
        assert_eq!(common_suffix_len(&a, &b), 5);
    }

    #[test]
    fn test_common_suffix_len_first_is_suffix_of_second() {
        let a = string_to_char_vec("ion");
        let b = string_to_char_vec("action");
        assert_eq!(common_suffix_len(&a, &b), 3);
    }

    #[test]
    fn test_common_suffix_len_second_is_suffix_of_first() {
        let a = string_to_char_vec("nation");
        let b = string_to_char_vec("ion");
        assert_eq!(common_suffix_len(&a, &b), 3);
    }

    #[test]
    fn test_common_suffix_len_different_cases() {
        let a = string_to_char_vec("Ending");
        let b = string_to_char_vec("ending");
        assert_eq!(common_suffix_len(&a, &b), 5); // Case-sensitive
    }

    #[test]
    fn test_common_suffix_len_with_special_chars() {
        let a = string_to_char_vec("abc!@#");
        let b = string_to_char_vec("123!@#");
        assert_eq!(common_suffix_len(&a, &b), 3);
    }

    #[test]
    fn test_common_suffix_len_one_char_common() {
        let a = string_to_char_vec("ax");
        let b = string_to_char_vec("bx");
        assert_eq!(common_suffix_len(&a, &b), 1);
    }

    #[test]
    fn test_common_suffix_len_one_char_different() {
        let a = string_to_char_vec("xa");
        let b = string_to_char_vec("xb");
        assert_eq!(common_suffix_len(&a, &b), 0);
    }

    #[test]
    fn test_find_middle_snake_empty_strings() {
        let (a, b, mut vf, mut vr) = init_test_find_middle_snake("", "");
        assert_eq!(find_middle_snake(&a, &b, &mut vf, &mut vr), None);
    }

    #[test]
    fn test_find_middle_snake_identical_strings() {
        let (a, b, mut vf, mut vr) = init_test_find_middle_snake("abcde", "abcde");
        assert_eq!(
            find_middle_snake(&a, &b, &mut vf, &mut vr),
            Some((0, 0, 5, 5))
        );
    }

    #[test]
    fn test_find_middle_snake_long_identical_strings() {
        let (a, b, mut vf, mut vr) =
            init_test_find_middle_snake("longidenticalstring", "longidenticalstring");
        assert_eq!(
            find_middle_snake(&a, &b, &mut vf, &mut vr),
            Some((0, 0, 19, 19))
        );
    }

    #[test]
    fn test_find_middle_snake_completely_different() {
        let (a, b, mut vf, mut vr) = init_test_find_middle_snake("aaaaa", "bbbbb");
        assert_eq!(find_middle_snake(&a, &b, &mut vf, &mut vr), None);
    }

    #[test]
    fn test_find_middle_snake_simple_edit() {
        let (a, b, mut vf, mut vr) = init_test_find_middle_snake("BC", "XBC");
        assert_eq!(
            find_middle_snake(&a, &b, &mut vf, &mut vr),
            Some((0, 1, 2, 3))
        );
    }

    #[test]
    fn test_find_middle_snake() {
        let (a, b, mut vf, mut vr) = init_test_find_middle_snake("ABCABBA", "CBABAC");
        assert_eq!(
            find_middle_snake(&a, &b, &mut vf, &mut vr),
            Some((3, 2, 5, 4))
        );
    }

    fn init_test_find_middle_snake(a: &str, b: &str) -> (Vec<char>, Vec<char>, V, V) {
        let a = string_to_char_vec(a);
        let b = string_to_char_vec(b);
        let max_d = max_d(a.len(), b.len());
        let vf = V::new(max_d);
        let vb = V::new(max_d);
        (a, b, vf, vb)
    }

    fn string_to_char_vec(s: &str) -> Vec<char> {
        s.chars().collect()
    }
}

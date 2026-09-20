//! Colinear chaining -- `modules/colinear_solver.py`.
//!
//! Algorithm 15.1 from Mäkinen et al., *Genome scale algorithmic design*. The
//! reference ships two implementations and picks between them on size:
//! `read_coverage` (O(n^2), used when a chromosome has fewer than 90 mems) and
//! `n_logn_read_coverage` (range-max trees). Both must be ported: the choice is
//! made per chromosome per read, so real data exercises both.
//!
//! THE TIE-BREAKS ARE THE CONTRACT. Three of them, none obvious, all reproduced
//! here and each the kind of thing that silently changes one read in a
//! thousand:
//!
//!   * `max(reversed(T_values), key=...)` -- iterating REVERSED means that
//!     among equal scores Python returns the one it meets first in reverse
//!     order, i.e. the LARGEST j_prime. A forward `max` would pick the
//!     smallest and be wrong.
//!   * `max_both([C_a, C_b])` returns `(index, value)` and ties go to index 0,
//!     so C_a wins a draw against C_b.
//!   * `argmax` is `max(enumerate(...), key=x[1])`, which returns the FIRST
//!     maximal index.

#[derive(Debug, Clone, PartialEq)]
pub struct Mem {
    pub x: i64,
    pub y: i64,
    pub c: i64,
    pub d: i64,
    pub val: i64,
    pub j: i64,
    pub exon_part_id: String,
}

/// `max(enumerate(it), key=x[1])` -- first maximal index.
fn argmax(v: &[i64]) -> usize {
    let mut best = 0usize;
    for (i, x) in v.iter().enumerate() {
        if *x > v[best] {
            best = i;
        }
    }
    best
}

/// `all_solutions_c_max_indicies`: every index whose value equals the max.
fn all_max_indices(c: &[i64], c_max: i64) -> Vec<usize> {
    c.iter()
        .enumerate()
        .filter(|(_, v)| **v == c_max)
        .map(|(i, _)| i)
        .collect()
}

/// `reconstruct_all_solutions` in the non-mam mode: walk the traceback from
/// each maximal index, emitting mems[index-1] because the vectors are shifted
/// by one.
fn reconstruct_all(
    mems: &[Mem],
    max_indices: &[usize],
    trace: &[usize],
) -> Vec<Vec<Mem>> {
    let mut out = Vec::with_capacity(max_indices.len());
    for &start in max_indices {
        let mut idx = start;
        let mut sol = Vec::new();
        while idx > 0 {
            sol.push(mems[idx - 1].clone());
            idx = trace[idx];
        }
        sol.reverse();
        out.push(sol);
    }
    out
}

/// `read_coverage(mems, max_intron)`.
///
/// `mems` must already be sorted by `y`, which `get_mems_from_input` guarantees.
pub fn read_coverage(mems: &[Mem], max_intron: i64) -> (Vec<Vec<Mem>>, i64) {
    let n = mems.len();
    let mut c = vec![0i64; n + 1];
    let mut trace = vec![0usize; n + 1];

    for j in 0..n {
        let v = &mems[j];

        // T: predecessors that end strictly before v starts on the read.
        // Python: max(reversed(T_values)) -> largest j_prime among ties.
        let mut t_idx: i64 = -1;
        let mut t_best: i64 = 0;
        {
            let mut found = false;
            for j_prime in 0..j {
                let m = &mems[j_prime];
                if m.d < v.c && v.y - m.y < max_intron {
                    let val = c[j_prime + 1];
                    // reversed iteration + `max` keeps the FIRST maximum seen
                    // in reverse, i.e. the largest j_prime. Scanning forward,
                    // that means >= rather than >.
                    if !found || val >= t_best {
                        t_best = val;
                        t_idx = j_prime as i64;
                        found = true;
                    }
                }
            }
            if !found {
                t_best = 0;
                t_idx = -1;
            }
        }

        // I: predecessors overlapping v on the read, scored with the chord
        // difference (v.d - m.d) added.
        let mut i_idx: i64 = -1;
        let mut i_best: i64 = 0;
        let c_b;
        {
            let mut found = false;
            for j_prime in 0..j {
                let m = &mems[j_prime];
                if v.c <= m.d && m.d <= v.d && v.y - m.y < max_intron {
                    let val = c[j_prime + 1] + (v.d - m.d);
                    if !found || val >= i_best {
                        i_best = val;
                        i_idx = j_prime as i64;
                        found = true;
                    }
                }
            }
            if found {
                c_b = i_best;
            } else {
                i_idx = -1;
                c_b = 0;
            }
        }

        let c_a = (v.d - v.c + 1) + t_best;

        // max_both([C_a, C_b]) -- ties go to C_a
        let (which, value) = if c_a >= c_b { (0, c_a) } else { (1, c_b) };
        c[j + 1] = value;
        let j_prime = if which == 0 { t_idx } else { i_idx };
        trace[j + 1] = if j_prime < 0 { 0 } else { (j_prime + 1) as usize };
    }

    let c_max = c[argmax(&c)];
    let idxs = all_max_indices(&c, c_max);
    let solutions = reconstruct_all(mems, &idxs, &trace);
    (solutions, c_max)
}

// ---------------------------------------------------------------------------
// The n log n path: `n_logn_read_coverage` + `range_query_max_search_tree.py`.
//
// Chosen when a chromosome has 90 or more mems. RARE in practice -- measured on
// 20 000 Drosophila reads it is a tiny fraction of chaining calls, and on SIRV
// it never fires at all (PORTING.md Finding 29). Rare is precisely why it must
// be ported against a recorded oracle rather than by eye: a bug here would show
// up on a handful of reads in a large run and nowhere else.
// ---------------------------------------------------------------------------

const NEG: i64 = -(1i64 << 32);

#[derive(Clone, Debug)]
struct TNode {
    d: i64,
    j: i64,
    cj: i64,
    j_max: i64,
}

/// `make_leafs_power_of_2`: one leaf per mem, plus a sentinel at d = -1, padded
/// with negative-index dummies up to a power of two, then sorted by `d`.
///
/// Python's `sorted` is stable, so leaves with equal `d` keep the order they
/// were appended in: sentinel first, then mems by index, then the padding.
fn make_leafs(mems: &[Mem]) -> Vec<TNode> {
    let mut nodes = vec![TNode { d: -1, j: -1, cj: NEG, j_max: -1 }];
    for m in mems {
        nodes.push(TNode { d: m.d, j: m.j, cj: NEG, j_max: m.j });
    }
    // grow to the next power of two
    let mut target = 1usize;
    while target < nodes.len() {
        target <<= 1;
    }
    let remainder = target - nodes.len();
    for i in 0..remainder {
        let neg = -(i as i64) - 2;
        nodes.push(TNode { d: -1, j: neg, cj: NEG, j_max: neg });
    }
    nodes.sort_by_key(|n| n.d); // stable, matching Python's sorted
    nodes
}

struct Rmq {
    tree: Vec<TNode>,
    n: usize,
}

impl Rmq {
    fn new(leafs: &[TNode]) -> Self {
        let n = leafs.len();
        let blank = TNode { d: 0, j: 0, cj: 0, j_max: 0 };
        let mut tree = vec![blank; 2 * n];
        for i in 0..n {
            tree[n + i] = leafs[i].clone();
        }
        for i in (1..n).rev() {
            // `max(..., key=d)`; Python's max keeps the FIRST on a tie, which
            // is the left child.
            let l = tree[2 * i].clone();
            let r = tree[2 * i + 1].clone();
            let pick = if r.d > l.d { r } else { l };
            tree[i] = pick;
        }
        Rmq { tree, n }
    }

    /// `update(tree, leaf_pos, value, n)`.
    fn update(&mut self, leaf_pos: usize, value: i64) {
        let mut pos = leaf_pos + self.n;
        self.tree[pos].cj = value;
        while pos > 1 {
            pos >>= 1;
            // max(sorted([l, r], key=j_max, desc), key=Cj): among equal Cj the
            // larger j_max wins.
            let l = self.tree[2 * pos].clone();
            let r = self.tree[2 * pos + 1].clone();
            let (first, second) = if r.j_max > l.j_max { (r, l) } else { (l, r) };
            let best = if second.cj > first.cj { second } else { first };
            self.tree[pos].cj = best.cj;
            self.tree[pos].j_max = best.j_max;
        }
    }

    /// `range_query(tree, l, r, n)` -> (C_max, j_max, node_pos)
    fn range_query(&self, l: i64, r: i64) -> (i64, i64, usize) {
        debug_assert!(l <= r);
        let (mut l_pos, mut r_pos) = (1usize, 1usize);
        let mut v = 1usize;
        let mut v_prime: Vec<usize> = Vec::new();
        let mut v_biss: Vec<usize> = Vec::new();

        loop {
            let l_l = 2 * l_pos;
            let l_r = 2 * l_pos + 1;
            let r_l = 2 * r_pos;
            let r_r = 2 * r_pos + 1;

            r_pos = if r >= self.tree[r_l].d { r_r } else { r_l };
            l_pos = if l > self.tree[l_l].d { l_r } else { l_l };

            if l_pos == r_pos {
                v = l_pos;
            }
            if l_pos != r_pos && v != (l_pos >> 1) {
                if r_pos == r_r {
                    push_unique(&mut v_biss, r_l);
                }
                if l_pos == l_l {
                    push_unique(&mut v_prime, l_r);
                }
            }
            if l_pos >= self.n && r_pos >= self.n {
                break;
            }
        }

        if self.tree[r_pos].d == r {
            push_unique(&mut v_biss, r_pos);
        }
        if self.tree[l_pos].d >= l {
            push_unique(&mut v_prime, l_pos);
        }

        let pick = |set: &Vec<usize>| -> Option<usize> {
            if set.is_empty() {
                return None;
            }
            // sorted by j_max descending, then max by Cj -> largest j_max among
            // equal Cj. V_prime/V_biss are Python SETS of ints, but the sort
            // makes the result independent of their iteration order.
            let mut s: Vec<usize> = set.clone();
            s.sort_by(|a, b| self.tree[*b].j_max.cmp(&self.tree[*a].j_max));
            let mut best = s[0];
            for &x in &s[1..] {
                if self.tree[x].cj > self.tree[best].cj {
                    best = x;
                }
            }
            Some(best)
        };

        let vl = pick(&v_prime);
        let vr = pick(&v_biss);
        let v_max_pos = match (vl, vr) {
            (Some(a), Some(b)) => {
                let (first, second) = if self.tree[b].j_max > self.tree[a].j_max { (b, a) } else { (a, b) };
                if self.tree[second].cj > self.tree[first].cj { second } else { first }
            }
            (Some(a), None) => a,
            (None, Some(b)) => b,
            // the reference prints "BUG" and exits here
            (None, None) => return (NEG, -1, 0),
        };
        (self.tree[v_max_pos].cj, self.tree[v_max_pos].j_max, v_max_pos)
    }
}

fn push_unique(v: &mut Vec<usize>, x: usize) {
    if !v.contains(&x) {
        v.push(x);
    }
}

/// `n_logn_read_coverage(mems)`.
pub fn n_logn_read_coverage(mems: &[Mem]) -> (Vec<Vec<Mem>>, i64) {
    let t_leafs = make_leafs(mems);
    let n = t_leafs.len();
    let mut t = Rmq::new(&t_leafs);
    let mut i_tree = Rmq::new(&make_leafs(mems));

    let mut mem_to_leaf: std::collections::HashMap<i64, usize> = Default::default();
    for (i, l) in t_leafs.iter().enumerate() {
        mem_to_leaf.insert(l.j, i);
    }

    let mut c = vec![0i64; mems.len() + 1];
    let mut trace = vec![0usize; mems.len() + 1];

    t.update(0, 0);
    i_tree.update(0, 0);
    let _ = n;

    for (j, mem) in mems.iter().enumerate() {
        let leaf = match mem_to_leaf.get(&(j as i64)) {
            Some(&x) => x,
            None => continue,
        };
        let (t_max, j_prime_a, _) = t.range_query(-1, mem.c - 1);
        let c_a = t_max + mem.d - mem.c + 1;
        let (i_max, j_prime_b, _) = i_tree.range_query(mem.c, mem.d);
        let c_b = i_max + mem.d;

        let (which, value) = if c_a >= c_b { (0, c_a) } else { (1, c_b) };
        c[j + 1] = value;
        let j_prime = if which == 0 { j_prime_a } else { j_prime_b };

        trace[j + 1] = if j_prime < 0 || value == 0 { 0 } else { (j_prime + 1) as usize };

        t.update(leaf, value);
        i_tree.update(leaf, value - mem.d);
    }

    let c_max = c[argmax(&c)];
    let idxs = all_max_indices(&c, c_max);
    (reconstruct_all(mems, &idxs, &trace), c_max)
}

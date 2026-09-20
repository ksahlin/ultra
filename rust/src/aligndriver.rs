//! `align.annotate_guaranteed_optimal_bound` and `align.find_exons`.

use crate::colinear::{Mam, Mem};
use crate::index::Key;
use std::collections::{BTreeMap, BTreeSet};

/// `annotate_guaranteed_optimal_bound`.
///
/// Splits each chromosome's mems wherever the genomic gap exceeds the allowed
/// intron, re-indexes `j` within each piece, and annotates the piece with the
/// best coverage it could possibly reach -- the union of the read intervals the
/// mems cover. `align_single` uses that bound to skip whole chromosomes.
///
/// Returns (chr_id, instance) -> (tot_cov, is_rc, mems).
pub fn annotate_guaranteed_optimal_bound(
    mems: &BTreeMap<u64, Vec<Mem>>,
    is_rc: bool,
    max_intron_chr: &BTreeMap<u64, i64>,
    max_global_intron: i64,
) -> BTreeMap<(u64, u64), (i64, bool, Vec<Mem>)> {
    let mut instances: Vec<((u64, u64), Vec<Mem>)> = Vec::new();

    for (chr_id, all_mems) in mems.iter() {
        let max_allowed_intron = std::cmp::min(
            max_intron_chr.get(chr_id).copied().unwrap_or(0) + 20000,
            max_global_intron,
        );
        let mut inst_idx: u64 = 0;
        if all_mems.len() > 1 {
            let mut j = 0i64;
            let mut cur: Vec<Mem> = Vec::new();
            let mut m1 = all_mems[0].clone();
            m1.j = j;
            cur.push(m1);
            for w in all_mems.windows(2) {
                let (a, b) = (&w[0], &w[1]);
                if b.x - a.y > max_allowed_intron {
                    instances.push(((*chr_id, inst_idx), std::mem::take(&mut cur)));
                    inst_idx += 1;
                    j = 0;
                    let mut m2 = b.clone();
                    m2.j = j;
                    cur.push(m2);
                } else {
                    j += 1;
                    let mut m2 = b.clone();
                    m2.j = j;
                    cur.push(m2);
                }
            }
            instances.push(((*chr_id, inst_idx), cur));
        } else {
            instances.push(((*chr_id, inst_idx), all_mems.clone()));
        }
    }

    let mut out = BTreeMap::new();
    for ((chr_id, inst), ms) in instances {
        // Sweep starts and stops. The reference builds `starts + stops` and
        // sorts by position with a STABLE sort, so a start and a stop at the
        // same coordinate keep starts first -- which closes intervals one base
        // later than a stops-first order would.
        let mut events: Vec<(u8, i64)> = Vec::with_capacity(ms.len() * 2);
        for m in &ms {
            events.push((0, m.c)); // 0 = start, listed first
        }
        for m in &ms {
            events.push((1, m.d)); // 1 = stop
        }
        events.sort_by_key(|e| e.1); // stable: preserves starts-before-stops
        if events.is_empty() {
            out.insert((chr_id, inst), (0, is_rc, ms));
            continue;
        }
        let mut active_start = events[0].1;
        let mut nr_active = 1i64;
        let mut tot_cov = 0i64;
        for &(site, pos) in &events[1..] {
            if nr_active == 0 {
                active_start = pos;
            }
            if site == 1 {
                nr_active -= 1;
            } else {
                nr_active += 1;
            }
            if nr_active == 0 {
                tot_cov += pos - active_start + 1; // mem coordinates are inclusive
            }
        }
        out.insert((chr_id, inst), (tot_cov, is_rc, ms));
    }
    out
}

/// `help_functions.find_all_paths`, depth-first over a coverage graph.
///
/// The reference iterates `set(graph[start]).difference(path)`, so the order
/// candidates are pushed depends on CPython set iteration, and `find_exons`
/// reads only `paths[0]`. Finding 33's family: the port iterates in sorted
/// order instead.
fn find_all_paths(graph: &BTreeMap<i64, Vec<i64>>, start: i64, end: i64) -> Vec<Vec<i64>> {
    let mut paths = Vec::new();
    let mut queue: Vec<(i64, Vec<i64>)> = vec![(start, Vec::new())];
    while let Some((node, path)) = queue.pop() {
        let mut path = path;
        path.push(node);
        if node == end {
            paths.push(path.clone());
        }
        if let Some(next) = graph.get(&node) {
            let seen: BTreeSet<i64> = path.iter().cloned().collect();
            let mut cands: Vec<i64> = next.iter().cloned().collect::<BTreeSet<i64>>()
                .into_iter().filter(|n| !seen.contains(n)).collect();
            // The reference pushes onto a LIFO stack fed from `set(...)`, so
            // which path reaches `end` first -- and `find_exons` reads only
            // paths[0] -- depends on CPython set iteration. Finding 33's
            // family. Ascending was measured to match more often than
            // descending (597 vs 595 of 600 recorded calls), so ascending it is;
            // the residual is budgeted in tests/aligndriver_oracle.rs.
            cands.sort_unstable();
            for n in cands {
                queue.push((n, path.clone()));
            }
        }
    }
    paths
}

pub struct Exons {
    pub exons: Vec<(i64, i64, i64, i64, i64)>,
    pub created_ref_seq: String,
    pub predicted_exons: Vec<(i64, i64)>,
    pub predicted_splices: Vec<(i64, i64)>,
    pub covered: i64,
}

/// `find_exons`.
pub fn find_exons(
    chr_id: u64,
    mam_solution: &[Mam],
    ref_exon_sequences: &BTreeMap<Key, String>,
    ref_segment_sequences: &BTreeMap<Key, String>,
    ref_flank_sequences: &BTreeMap<Key, String>,
    splice_pairs: &BTreeSet<(i64, i64)>,
) -> Exons {
    // split at annotated intron sites
    let mut parts: Vec<Vec<Mam>> = Vec::new();
    if mam_solution.len() > 1 {
        let mut prev = 0usize;
        for i in 0..mam_solution.len() - 1 {
            let (m1, m2) = (&mam_solution[i], &mam_solution[i + 1]);
            if splice_pairs.contains(&(m1.y, m2.x)) {
                parts.push(mam_solution[prev..=i].to_vec());
                prev = i + 1;
            }
        }
        if prev < mam_solution.len() {
            parts.push(mam_solution[prev..].to_vec());
        }
    } else {
        parts.push(mam_solution.to_vec());
    }

    // split again at gaps too large to be deletions
    let mut parts2: Vec<Vec<Mam>> = Vec::new();
    for part in parts {
        if part.len() > 1 {
            let mut prev = 0usize;
            for i in 0..part.len() - 1 {
                if part[i + 1].x - part[i].y > 10 {
                    parts2.push(part[prev..=i].to_vec());
                    prev = i + 1;
                }
            }
            if prev < part.len() {
                parts2.push(part[prev..].to_vec());
            }
        } else {
            parts2.push(part);
        }
    }

    let mut exons: Vec<(i64, i64, i64, i64, i64)> = Vec::new();
    for part in &parts2 {
        let adjacent = part.len() == 1
            || part.windows(2).all(|w| w[0].y == w[1].x);
        if adjacent {
            for m in part {
                exons.push((m.x, m.y, m.c, m.d, m.ref_chr_id));
            }
            continue;
        }
        // look for a chain of annotated exons / segments spanning the region
        let mut all_points: Vec<i64> = Vec::new();
        let mut segm: BTreeMap<(u64, i64, i64), &Mam> = BTreeMap::new();
        let mut start_points: BTreeMap<i64, &Mam> = BTreeMap::new();
        let mut end_points: BTreeMap<i64, &Mam> = BTreeMap::new();
        for m in part {
            all_points.push(m.x);
            all_points.push(m.y);
            segm.insert((chr_id, m.x, m.y), m);
            start_points.insert(m.x, m);
            end_points.insert(m.y, m);
        }
        let mut sorted_pts = all_points.clone();
        sorted_pts.sort_unstable();
        let mut cover: BTreeMap<i64, Vec<i64>> = BTreeMap::new();
        for &p1 in &sorted_pts {
            let e = cover.entry(p1).or_default();
            for &p2 in &sorted_pts {
                let key: Key = (chr_id, p1 as u64, p2 as u64);
                if ref_exon_sequences.contains_key(&key) || segm.contains_key(&(chr_id, p1, p2)) {
                    e.push(p2);
                }
            }
        }
        let paths = find_all_paths(&cover, part[0].x, part[part.len() - 1].y);
        if let Some(path) = paths.first() {
            for w in path.windows(2) {
                let (p1, p2) = (w[0], w[1]);
                if let Some(m) = segm.get(&(chr_id, p1, p2)) {
                    exons.push((m.x, m.y, m.c, m.d, m.ref_chr_id));
                } else {
                    let (c, x) = match start_points.get(&p1) {
                        Some(m) => (m.c, m.x),
                        None => { let m = end_points[&p1]; (m.d, m.y) }
                    };
                    let (d, y) = match start_points.get(&p2) {
                        Some(m) => (m.c, m.x),
                        None => { let m = end_points[&p2]; (m.d, m.y) }
                    };
                    exons.push((x, y, c, d, chr_id as i64));
                }
            }
        } else {
            for m in part {
                exons.push((m.x, m.y, m.c, m.d, m.ref_chr_id));
            }
        }
    }

    let mut chained = String::new();
    let mut predicted_exons: Vec<(i64, i64)> = Vec::new();
    let mut prev_y = -1i64;
    let mut covered = 0i64;
    for &(x, y, c, d, seq_id) in &exons {
        let key: Key = (seq_id as u64, x as u64, y as u64);
        let seq = ref_exon_sequences
            .get(&key)
            .or_else(|| ref_segment_sequences.get(&key))
            .or_else(|| ref_flank_sequences.get(&key));
        if let Some(s) = seq {
            covered += d - c + 1;
            chained.push_str(s);
        }
        if prev_y >= x {
            // adjacent segments: extend the previous exon rather than splitting
            let last = predicted_exons.len() - 1;
            predicted_exons[last].1 = y;
        } else {
            predicted_exons.push((x, y));
        }
        prev_y = y;
    }
    let predicted_splices: Vec<(i64, i64)> = predicted_exons
        .windows(2)
        .map(|w| (w[0].1, w[1].0))
        .collect();

    Exons { exons, created_ref_seq: chained, predicted_exons, predicted_splices, covered }
}

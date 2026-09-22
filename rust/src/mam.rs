//! The MAM layer -- `modules/classify_read_with_mams.py`.
//!
//! DIVERGENCE (Finding 33): the reference derives part of its ordering from
//! CPython's `set` iteration, which the port replaces with a total order.
//! Measured cost: 6 reads of 10 000 on SIRV, 0 of 2 000 on Drosophila.
//! Decision taken by the author on 2026-09-20.

use crate::colinear::Mam;
use crate::edlib;

/// `cigar_to_accuracy`: matches divided by total aligned length.
///
/// Note the denominator counts every operation, `=`, `X`, `I` and `D` alike,
/// so it is not identity-over-aligned-columns in the usual sense.
pub fn cigar_to_accuracy(cigar: &str) -> f64 {
    let (mut matches, mut aln_len) = (0f64, 0f64);
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
        } else {
            let n: f64 = num.parse().unwrap_or(0.0);
            num.clear();
            if ch == '=' {
                matches += n;
            }
            aln_len += n;
        }
    }
    if aln_len == 0.0 {
        return 0.0;
    }
    matches / aln_len
}

/// `classify_read_with_mams.edlib_alignment` with `task='path'`.
/// Returns (locations, edit_distance, accuracy).
pub fn edlib_path(query: &[u8], target: &[u8], k: i32) -> (Vec<(i64, i64)>, i32, f64) {
    let r = edlib::align(query, target, edlib::MODE_HW, edlib::TASK_PATH, k);
    if r.edit_distance == -1 {
        // the reference returns the literal [0, 0] here, not an empty list
        return (vec![(0, 0)], -1, 0.0);
    }
    let acc = r.cigar.as_deref().map(cigar_to_accuracy).unwrap_or(0.0);
    (r.locations, r.edit_distance, acc)
}

/// `add_segment_to_mam`, appending to `out`.
///
/// Three independent blocks, all reproduced including their quirks:
///
///   * segments of 9 bp or more are aligned with `k = 0.4*min(len)`, and each
///     distinct start position contributes one MAM -- but the "keep the best
///     per start" bookkeeping is buggy in a way that is contract:
///     `max_score` is overwritten with the LAST accepted score rather than the
///     maximum, so a later lower-scoring hit at a repeated start can still
///     displace nothing while a higher one is admitted. Faithful.
///   * segments of 5..8 bp require an exact match (`k = 0`).
///   * independently, a segment at least 0.8x the read length is aligned the
///     other way round (read against segment), and on success the MAM spans
///     the WHOLE read, `start, stop = 0, len(read)-1`.
#[allow(clippy::too_many_arguments)]
pub fn add_segment_to_mam(
    read_seq: &[u8],
    ref_chr_id: i64,
    exon_seq: &[u8],
    e_start: i64,
    e_stop: i64,
    segm_id: &str,
    min_acc: f64,
    annot_label: &str,
    out: &mut Vec<Mam>,
) {
    let mam_id = format!("{segm_id}{annot_label}");

    if e_stop - e_start >= 5 {
        if e_stop - e_start >= 9 {
            let k = (0.4 * std::cmp::min(read_seq.len(), exon_seq.len()) as f64) as i32;
            let (locations, ed, acc) = edlib_path(exon_seq, read_seq, k);
            if ed >= 0 && acc > min_acc {
                let mut considered: Vec<i64> = Vec::new();
                let mut max_score = 0.0f64;
                for &(start, stop) in &locations {
                    let msl = (stop - start + 1) as f64;
                    let acc_approx = f64::max((msl - ed as f64) / msl, acc);
                    let score = acc_approx * msl;
                    if considered.contains(&start) {
                        if score <= max_score {
                            continue;
                        }
                    } else {
                        max_score = score;
                    }
                    out.push(Mam {
                        x: e_start, y: e_stop, c: start, d: stop,
                        val: score, j: 0, min_segment_length: msl as i64,
                        mam_id: mam_id.clone(), ref_chr_id,
                    });
                    considered.push(start);
                    max_score = score; // yes, unconditionally -- see above
                }
            }
        } else {
            // 5..8 bp: exact match only, or the noise swamps everything
            let (locations, ed, _acc) = edlib_path(exon_seq, read_seq, 0);
            if ed == 0 {
                let score = exon_seq.len() as f64;
                let mut considered: Vec<i64> = Vec::new();
                for &(start, stop) in &locations {
                    if considered.contains(&start) {
                        continue;
                    }
                    out.push(Mam {
                        x: e_start, y: e_stop, c: start, d: stop,
                        val: score, j: score as i64, min_segment_length: 0,
                        mam_id: mam_id.clone(), ref_chr_id,
                    });
                    considered.push(start);
                }
            }
        }
    }

    // Independent of the above: the read may be contained within the segment.
    if (e_stop - e_start) as f64 >= 0.8 * read_seq.len() as f64 {
        let k = (0.4 * std::cmp::min(read_seq.len(), exon_seq.len()) as f64) as i32;
        let (locations, ed, acc) = edlib_path(read_seq, exon_seq, k);
        if ed >= 0 {
            let (start, stop) = locations[0];
            let msl = (stop - start + 1) as f64;
            let score = acc * msl;
            if acc > min_acc {
                out.push(Mam {
                    x: e_start, y: e_stop,
                    c: 0, d: read_seq.len() as i64 - 1, // the WHOLE read
                    val: score, j: 0, min_segment_length: msl as i64,
                    mam_id: mam_id.clone(), ref_chr_id,
                });
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The plumbing: get_unique_exon_and_flank_locations and
// get_unique_segment_and_flank_choordinates.
// ---------------------------------------------------------------------------

use crate::colinear::Mem;
use crate::index::Key;
use std::collections::{BTreeMap, BTreeSet};

/// Hit coordinates for one read, as `get_unique_exon_and_flank_locations`
/// produces them.
#[derive(Debug, Default)]
pub struct HitLocations {
    /// DIVERGENCE (Finding 33): the reference produces this as
    /// `list(set(...))` sorted by start ALONE, so entries tied on start come
    /// out in CPython set-iteration order. The port sorts by the full
    /// (start, chr, stop) triple instead. The SET is identical either way;
    /// only the order of ties differs.
    pub segment_hit_locations: Vec<(u64, u64, u64)>,
    pub flank_hit_locations: Vec<(u64, u64, u64)>,
    /// (chr, start, stop) -> [ref_start, ref_stop, read_start, read_stop]
    pub partial_segment_hit_locations: BTreeMap<(u64, u64, u64), [i64; 4]>,
    pub partial_flank_hit_locations: BTreeMap<(u64, u64, u64), [i64; 4]>,
    pub choord_to_exon_id: BTreeMap<(u64, u64, u64), Key>,
    pub first_part_stop: u64,
    pub last_part_start: u64,
}

fn is_overlapping(a1: i64, a2: i64, b1: i64, b2: i64) -> bool {
    (a1 <= b1 && b1 <= a2) || (a1 <= b2 && b2 <= a2) || (b1 <= a1 && a1 <= b2) || (b1 <= a2 && a2 <= b2)
}

/// `get_unique_exon_and_flank_locations(solution, parts_to_segments)`.
pub fn get_unique_exon_and_flank_locations(
    solution: &[Mem],
    parts_to_segments: &BTreeMap<Key, Vec<Key>>,
) -> HitLocations {
    let mut h = HitLocations {
        first_part_stop: u64::from(1u32 << 31) * 2, // 2**32, as the reference writes it
        last_part_start: 0,
        ..Default::default()
    };
    let mut seen_segments: BTreeSet<(u64, u64, u64)> = BTreeSet::new();

    for mem in solution {
        let mut it = mem.exon_part_id.split('^');
        let (c, a, b) = match (it.next(), it.next(), it.next()) {
            (Some(c), Some(a), Some(b)) => (c, a, b),
            _ => continue,
        };
        let (chr, rs, re): (u64, u64, u64) = match (c.parse(), a.parse(), b.parse()) {
            (Ok(x), Ok(y), Ok(z)) => (x, y, z),
            _ => continue,
        };
        let key: Key = (chr, rs, re);

        match parts_to_segments.get(&key) {
            None => {
                // not a part -> it is a flank
                h.flank_hit_locations.push((chr, rs, re));
                h.partial_flank_hit_locations
                    .entry((chr, rs, re))
                    .and_modify(|v| {
                        v[1] = mem.y;
                        v[3] = mem.d;
                    })
                    .or_insert([mem.x, mem.y, mem.c, mem.d]);
            }
            Some(segs) => {
                if re <= h.first_part_stop {
                    h.first_part_stop = re;
                }
                if rs >= h.last_part_start {
                    h.last_part_start = rs;
                }
                for s in segs {
                    let (_schr, ss, se) = *s;
                    if is_overlapping(ss as i64, se as i64, mem.x, mem.y) {
                        h.choord_to_exon_id.insert((chr, ss, se), (s.0, ss, se));
                        seen_segments.insert((chr, ss, se));
                        h.partial_segment_hit_locations
                            .entry((chr, ss, se))
                            .and_modify(|v| {
                                v[1] = mem.y;
                                v[3] = mem.d;
                            })
                            .or_insert([mem.x, mem.y, mem.c, mem.d]);
                    }
                }
            }
        }
    }

    // The reference: list(set(...)) then sort by x[1]. We sort by the full
    // triple -- Finding 33.
    let mut v: Vec<(u64, u64, u64)> = seen_segments.into_iter().collect();
    v.sort_by_key(|t| (t.1, t.0, t.2));
    h.segment_hit_locations = v;
    h
}

/// `get_unique_segment_and_flank_choordinates`.
pub struct UniqueChoords {
    /// (chr, start, stop) -> the segment ids mapping to it
    pub segments: BTreeMap<(u64, u64, u64), BTreeSet<Key>>,
    pub segments_partial: BTreeMap<(u64, u64, u64), (u64, i64, i64, Key)>,
    pub flanks: BTreeSet<(u64, u64, u64)>,
    pub flanks_partial: BTreeMap<(u64, u64, u64), (u64, i64, i64)>,
}

pub fn get_unique_segment_and_flank_choordinates(
    h: &HitLocations,
    segment_to_gene: &BTreeMap<Key, BTreeSet<String>>,
    gene_to_small_segments: &BTreeMap<String, Vec<Key>>,
) -> UniqueChoords {
    let mut u = UniqueChoords {
        segments: BTreeMap::new(),
        segments_partial: BTreeMap::new(),
        flanks: BTreeSet::new(),
        flanks_partial: BTreeMap::new(),
    };

    if !h.segment_hit_locations.is_empty() {
        // every small segment of every gene touched by a hit, added up front
        let mut genes: BTreeSet<&String> = BTreeSet::new();
        for loc in &h.segment_hit_locations {
            if let Some(eid) = h.choord_to_exon_id.get(loc) {
                if let Some(gs) = segment_to_gene.get(eid) {
                    genes.extend(gs.iter());
                }
            }
        }
        for g in genes {
            if let Some(smalls) = gene_to_small_segments.get(g) {
                for s in smalls {
                    u.segments.entry(*s).or_default().insert(*s);
                }
            }
        }
    }

    for loc in &h.segment_hit_locations {
        if let Some(eid) = h.choord_to_exon_id.get(loc) {
            u.segments.entry(*loc).or_default().insert(*eid);
            if let Some(p) = h.partial_segment_hit_locations.get(loc) {
                u.segments_partial.insert(*loc, (loc.0, p[0], p[1], *eid));
            }
        }
    }

    for loc in &h.flank_hit_locations {
        u.flanks.insert(*loc);
        if let Some(p) = h.partial_flank_hit_locations.get(loc) {
            let (_c, rs, re) = *loc;
            let span = (re - rs) as f64;
            // a hit that starts well inside the flank, or ends well before it
            if (p[0] - rs as i64) as f64 > 0.05 * span {
                u.flanks_partial.insert(*loc, (loc.0, p[0], p[1]));
            }
            if (re as i64 - p[1]) as f64 > 0.05 * span {
                u.flanks_partial.insert(*loc, (loc.0, p[0], p[1]));
            }
        }
    }
    u
}

/// `classify_read_with_mams.main` -> (non_covered_regions, value, mam_solution)
///
/// Assembles the pieces above in the reference's order: full segments, then
/// full flanks, then partial segment hits, then partial flank hits; filter the
/// partial hits to the outermost ones; sort by `y` and re-index; chain.
#[allow(clippy::too_many_arguments)]
pub fn classify_read(
    solution: &[Mem],
    ref_segment_sequences: &BTreeMap<Key, String>,
    ref_flank_sequences: &BTreeMap<Key, String>,
    parts_to_segments: &BTreeMap<Key, Vec<Key>>,
    segment_to_gene: &BTreeMap<Key, BTreeSet<String>>,
    gene_to_small_segments: &BTreeMap<String, Vec<Key>>,
    read_seq: &[u8],
    min_acc: f64,
) -> (Vec<i64>, f64, Vec<Mam>) {
    let h = get_unique_exon_and_flank_locations(solution, parts_to_segments);
    let u = get_unique_segment_and_flank_choordinates(&h, segment_to_gene, gene_to_small_segments);

    let mut mams: Vec<Mam> = Vec::new();

    // full segments, in start order
    let mut seg_keys: Vec<&(u64, u64, u64)> = u.segments.keys().collect();
    seg_keys.sort_by_key(|t| (t.1, t.0, t.2));
    for loc in seg_keys {
        let key: Key = (loc.0, loc.1, loc.2);
        let seq = match ref_segment_sequences.get(&key) {
            Some(s) => s,
            None => continue,
        };
        // The reference does `all_segm_ids.pop()` -- an arbitrary element of a
        // SET. Finding 33: we take the smallest, deterministically.
        let segm_id = match u.segments.get(loc).and_then(|s| s.iter().next()) {
            Some(k) => format!("{:?}", k),
            None => continue,
        };
        add_segment_to_mam(read_seq, loc.0 as i64, seq.as_bytes(), loc.1 as i64, loc.2 as i64,
                           &segm_id, min_acc, "_full_segment", &mut mams);
    }

    // full flanks, in start order
    let mut flank_keys: Vec<&(u64, u64, u64)> = u.flanks.iter().collect();
    flank_keys.sort_by_key(|t| (t.1, t.0, t.2));
    for loc in flank_keys {
        let key: Key = (loc.0, loc.1, loc.2);
        let seq = match ref_flank_sequences.get(&key) {
            Some(s) => s,
            None => continue,
        };
        let fid = format!("flank_{}_{}", loc.1, loc.2);
        add_segment_to_mam(read_seq, loc.0 as i64, seq.as_bytes(), loc.1 as i64, loc.2 as i64,
                           &fid, min_acc, "_full_flank", &mut mams);
    }

    // Partial hits are only allowed at the read's own ends, and only outside
    // the span already covered by a valid full hit.
    let (first_valid_stop, last_valid_start) = if mams.is_empty() {
        (-1i64, 1i64 << 32)
    } else {
        (mams.iter().map(|m| m.y).min().unwrap(), mams.iter().map(|m| m.x).max().unwrap())
    };
    let final_first_stop = std::cmp::max(h.first_part_stop as i64, first_valid_stop);
    let final_last_start = std::cmp::min(h.last_part_start as i64, last_valid_start);

    let mut tried: BTreeSet<String> = BTreeSet::new();
    for (loc, (_c, s_start, _s_stop)) in u.segments_partial.iter().map(|(k, v)| (k, (v.0, v.1, v.2))) {
        let key: Key = (loc.0, loc.1, loc.2);
        let seq = match ref_segment_sequences.get(&key) {
            Some(s) => s.as_bytes(),
            None => continue,
        };
        let (e_start, e_stop) = (loc.1 as i64, loc.2 as i64);
        let eid = u.segments_partial.get(loc).map(|v| format!("{:?}", v.3)).unwrap_or_default();
        // Python slicing CLAMPS: seq[n:] with n past the end is empty, and
        // seq[:n] with n past the end is the whole string. Rejecting those
        // cases instead -- as an early version of this port did -- silently
        // drops partial hits the reference keeps.
        if e_stop <= final_first_stop {
            let off = ((s_start - e_start).max(0) as usize).min(seq.len());
            let part = &seq[off..];
            let ps = String::from_utf8_lossy(part).into_owned();
            if part.len() > 5 && !tried.contains(&ps) {
                add_segment_to_mam(read_seq, loc.0 as i64, part, e_start, e_stop, &eid,
                                   min_acc, "_partial_segment_start", &mut mams);
                tried.insert(ps);
            }
        } else if final_last_start <= e_start {
            let s_stop = u.segments_partial.get(loc).map(|v| v.2).unwrap_or(0);
            let cut = seq.len() as i64 - (e_stop - (s_stop + 1));
            let cut = cut.max(0).min(seq.len() as i64) as usize;
            let part = &seq[..cut];
            let ps = String::from_utf8_lossy(part).into_owned();
            if part.len() > 5 && !tried.contains(&ps) {
                add_segment_to_mam(read_seq, loc.0 as i64, part, e_start, e_stop, &eid,
                                   min_acc, "_partial_segment_end", &mut mams);
                tried.insert(ps);
            }
        }
    }

    // Keep only the OUTERMOST partial hits. Note the reference's `elif`: an id
    // containing both "_start" and "_end" only ever updates the start side.
    let mut outmost_start = 1i64 << 32;
    let mut min_segment_end = 1i64 << 32;
    let mut outmost_end = 0i64;
    let mut max_segment_start = 0i64;
    for m in &mams {
        if m.mam_id.contains("_start") && m.c < outmost_start {
            outmost_start = m.c;
        } else if m.mam_id.contains("_end") && m.d > outmost_end {
            outmost_end = m.d;
        }
        if m.d < min_segment_end { min_segment_end = m.d; }
        if m.c > max_segment_start { max_segment_start = m.c; }
    }
    let outmost_start = std::cmp::min(outmost_start, min_segment_end);
    let outmost_end = std::cmp::max(outmost_end, max_segment_start);
    mams.retain(|m| {
        !(m.mam_id.contains("_start") && m.c > outmost_start)
            && !(m.mam_id.contains("_end") && m.d < outmost_end)
    });

    mams.sort_by_key(|m| m.y);
    for (j, m) in mams.iter_mut().enumerate() {
        m.j = j as i64;
    }


    if mams.is_empty() {
        return (Vec::new(), -1.0, Vec::new());
    }
    let (solution, value, _unique) = if mams.len() > 200 {
        crate::colinear::n_logn_read_coverage_mams(&mams, 5)
    } else {
        crate::colinear::read_coverage_mam_score(&mams, 20)
    };

    let mut non_covered = Vec::new();
    if !solution.is_empty() {
        non_covered.push(solution[0].c);
        for w in solution.windows(2) {
            non_covered.push(w[1].c - w[0].d - 1);
        }
        non_covered.push(read_seq.len() as i64 - solution[solution.len() - 1].d);
    }
    (non_covered, value, solution)
}

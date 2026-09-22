//! Aligning one read -- the Rust equivalent of the body of
//! `modules/align.py:align_single`.
//!
//! DELIBERATELY NOT A TRANSLATION OF ITS SHAPE. In Python `align_single` is the
//! entry point of a spawned worker process: it calls `import_data()` to load
//! the entire index into that process, then pulls batches off a
//! `Manager().Queue()` and pushes SAM records back. Two consequences the port
//! does not inherit:
//!
//!   * the index is loaded once PER WORKER (PORTING.md Finding 4), so memory
//!     scales with `--t`. That is what makes issue #15 OOM at `--t 48`. Here
//!     one `&Index` is shared by every thread.
//!   * the output order is worker-completion order and therefore not
//!     reproducible (Finding 32), so there is nothing to imitate: this returns
//!     records for one read and the caller keeps them in input order, which is
//!     deterministic.
//!
//! So what is ported is the LOGIC of the loop body, not its plumbing.

use crate::colinear::{Mam, Mem};
use crate::index::Index;
use crate::{aligndriver, colinear, edlib, mam, reads, samout};
use std::collections::BTreeMap;

pub struct Params {
    pub max_intron: i64,
    pub min_acc: f64,
    pub dropoff: f64,
    pub max_loc: f64,
    pub alignment_threshold: f64,
    pub non_covered_cutoff: i64,
    pub reduce_read_polya: usize,
}

/// One candidate alignment, before the primary/secondary decision.
struct Candidate {
    score: f64,
    genome_start: i64,
    genome_stop: i64,
    chr_id: u64,
    classification: &'static str,
    predicted_exons: Vec<(i64, i64)>,
    read_aln: String,
    ref_aln: String,
    annotated_to: String,
    is_rc: bool,
}

/// `get_exact_alignment`: parasail below 20 kb, edlib above.
fn exact_alignment(read_seq: &[u8], created_ref: &[u8], mam_sol_exons_length: i64) -> (String, String, f64) {
    let long = created_ref.len() > 20000
        || read_seq.len() > 20000
        || (read_seq.len() > 1000 && (read_seq.len() as i64) < mam_sol_exons_length / 10);
    if long {
        let r = edlib::align(read_seq, created_ref, edlib::MODE_HW, edlib::TASK_PATH, -1);
        // the reference rebuilds the gapped strings from the cigar; the score
        // is matches*2 - 2*edit_distance
        let (ra, rfa) = crate::samout::expand_cigar(
            r.cigar.as_deref().unwrap_or(""), read_seq, created_ref, r.locations.first().copied());
        let matches = ra.bytes().zip(rfa.bytes()).filter(|(a, b)| a == b).count() as f64 * 2.0;
        (ra, rfa, matches - 2.0 * r.edit_distance as f64)
    } else {
        let p = crate::parasail::sg_trace(read_seq, created_ref, 2, -2, 3, 1);
        let (ra, rfa) = crate::samout::expand_cigar_full(&p.cigar, read_seq, created_ref);
        (ra, rfa, p.score as f64)
    }
}

/// Align one read and return its SAM records, primary first.
pub fn align_read(
    ix: &Index,
    read_acc: &str,
    seq: &str,
    qual: Option<&str>,
    hits: &[String],
    hits_rc: &[String],
    p: &Params,
) -> Vec<String> {
    // the reference compresses sequence AND quality together, with to_len = 1
    // here (not 5 as in the namfinder preprocessing)
    let (seq_mod, qual_mod) = reads::remove_read_polya_ends_q(seq, qual, p.reduce_read_polya, 1);

    let mems = to_mems(hits);
    let mems_rc = to_mems(hits_rc);
    let max_intron_chr: BTreeMap<u64, i64> =
        ix.max_intron_chr.iter().map(|(k, v)| (*k, *v as i64)).collect();

    let ub = aligndriver::annotate_guaranteed_optimal_bound(&mems, false, &max_intron_chr, p.max_intron);
    let ub_rc = aligndriver::annotate_guaranteed_optimal_bound(&mems_rc, true, &max_intron_chr, p.max_intron);

    // sorted by the bound, descending -- the reference's `sorted(..., reverse=True)`
    let mut all: Vec<((u64, u64), (i64, bool, Vec<Mem>))> =
        ub.into_iter().chain(ub_rc).collect();
    all.sort_by(|a, b| b.1 .0.cmp(&a.1 .0));

    let mut chainings: Vec<(u64, Vec<Mem>, i64, bool)> = Vec::new();
    let mut best_value = 0i64;
    for ((chr_id, _), (bound, is_rc, chr_mems)) in all {
        if (bound as f64) < best_value as f64 * p.dropoff {
            break;
        }
        let max_allowed = std::cmp::min(
            max_intron_chr.get(&chr_id).copied().unwrap_or(0) + 20000, p.max_intron);
        let (sols, value) = if chr_mems.len() < 90 {
            colinear::read_coverage(&chr_mems, max_allowed)
        } else {
            colinear::n_logn_read_coverage(&chr_mems)
        };
        if value > best_value {
            best_value = value;
        }
        for s in sols {
            chainings.push((chr_id, s, value, is_rc));
        }
    }

    if chainings.is_empty() {
        return vec![samout::sam_record(read_acc, &seq_mod, qual_mod.as_deref(), "*",
            "unaligned", &[], "*", "*", "*", false, false, 0)];
    }

    // key: (-score, span) -- the reference's sort
    chainings.sort_by(|a, b| {
        let sa = -(a.2);
        let sb = -(b.2);
        sa.cmp(&sb).then_with(|| {
            let spa = a.1.last().map(|m| m.y).unwrap_or(0) - a.1.first().map(|m| m.x).unwrap_or(0);
            let spb = b.1.last().map(|m| m.y).unwrap_or(0) - b.1.first().map(|m| m.x).unwrap_or(0);
            spa.cmp(&spb)
        })
    });
    let best_chaining_score = chainings[0].2 as f64;

    let mut candidates: Vec<Candidate> = Vec::new();
    let mut seen_mam_solutions: Vec<Vec<Mam>> = Vec::new();

    for (i, (chr_id, mem_solution, chaining_score, is_rc)) in chainings.iter().enumerate() {
        if (*chaining_score as f64) / best_chaining_score < p.dropoff || i as f64 >= p.max_loc {
            break;
        }
        let read_seq: String = if *is_rc { crate::prefilter::revcomp(&seq_mod) } else { seq_mod.clone() };

        let (non_covered, mam_value, mam_solution) = mam::classify_read(
            mem_solution, &ix.ref_segment_sequences, &ix.ref_flank_sequences,
            &ix.parts_to_segments, &ix.segment_to_gene, &ix.gene_to_small_segments,
            read_seq.as_bytes(), p.min_acc);

        // the reference skips a chaining that yields a MAM solution already seen
        if seen_mam_solutions.contains(&mam_solution) {
            continue;
        }
        seen_mam_solutions.push(mam_solution.clone());
        if mam_value <= 0.0 {
            continue;
        }
        let mam_sol_exons_length: i64 = mam_solution.iter().map(|m| m.y - m.x).sum();

        let empty = Default::default();
        let pairs: std::collections::BTreeSet<(i64, i64)> = ix
            .all_splice_pairs_annotations.get(chr_id).unwrap_or(&empty)
            .keys().map(|(a, b)| (*a as i64, *b as i64)).collect();
        let ex = aligndriver::find_exons(*chr_id, &mam_solution, &ix.ref_exon_sequences,
            &ix.ref_segment_sequences, &ix.ref_flank_sequences, &pairs);


        let (classification, annotated_to) = classify_for(ix, *chr_id, &ex.predicted_splices);

        let largest_intron = mam_solution.windows(2)
            .map(|w| w[1].x - w[0].y).max().unwrap_or(0);
        let max_allowed = std::cmp::min(
            max_intron_chr.get(chr_id).copied().unwrap_or(0) + 20000, p.max_intron);
        if largest_intron > max_allowed && classification != "FSM" {
            continue;
        }

        let (read_aln, ref_aln, score) =
            exact_alignment(read_seq.as_bytes(), ex.created_ref_seq.as_bytes(), mam_sol_exons_length);
        if score < 2.0 * p.alignment_threshold * read_seq.len() as f64 && classification != "FSM" {
            continue;
        }

        let mut classification = classification;
        if non_covered.len() >= 3 {
            let internal_max = non_covered[1..non_covered.len() - 1].iter().copied().max().unwrap_or(0);
            if internal_max > p.non_covered_cutoff {
                classification = "Insufficient_junction_coverage_unclassified";
            }
        }

        candidates.push(Candidate {
            score,
            genome_start: mam_solution[0].x,
            genome_stop: mam_solution[mam_solution.len() - 1].y,
            chr_id: *chr_id,
            classification,
            predicted_exons: ex.predicted_exons,
            read_aln, ref_aln,
            annotated_to,
            is_rc: *is_rc,
        });
    }

    if candidates.is_empty() {
        return vec![samout::sam_record(read_acc, &seq_mod, qual_mod.as_deref(), "*",
            "unaligned", &[], "*", "*", "*", false, false, 0)];
    }

    // sorted by (-score, span, classification)
    candidates.sort_by(|a, b| {
        b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| (a.genome_stop - a.genome_start).cmp(&(b.genome_stop - b.genome_start)))
            .then_with(|| a.classification.cmp(b.classification))
    });
    let best = candidates[0].score;
    let more_than_one = candidates.len() > 1;

    let mut out = Vec::with_capacity(candidates.len());
    for (i, c) in candidates.iter().enumerate() {
        let (is_secondary, map_score) = if i == 0 {
            let ms = if more_than_one && candidates[1].score == best { 0 } else { 60 };
            (false, ms)
        } else {
            (true, 0)
        };
        let read_seq: String = if c.is_rc { crate::prefilter::revcomp(&seq_mod) } else { seq_mod.clone() };
        let read_qual = qual_mod.as_ref().map(|q| if c.is_rc {
            q.chars().rev().collect::<String>()
        } else {
            q.clone()
        });
        let chr_name = ix.id_to_chr.get(&c.chr_id).cloned().unwrap_or_else(|| "*".into());
        out.push(samout::sam_record(read_acc, &read_seq, read_qual.as_deref(), &chr_name,
            c.classification, &c.predicted_exons, &c.read_aln, &c.ref_aln,
            &c.annotated_to, c.is_rc, is_secondary, map_score));
    }
    out
}

fn classify_for(ix: &Index, chr_id: u64, predicted_splices: &[(i64, i64)]) -> (&'static str, String) {
    let empty_pairs = Default::default();
    let pairs_raw = ix.all_splice_pairs_annotations.get(&chr_id).unwrap_or(&empty_pairs);
    let pairs: BTreeMap<(i64, i64), std::collections::BTreeSet<String>> = pairs_raw
        .iter().map(|((a, b), t)| ((*a as i64, *b as i64), t.clone())).collect();
    let empty_sites = Default::default();
    let sites: std::collections::BTreeSet<i64> = ix
        .all_splice_sites_annotations.get(&chr_id).unwrap_or(&empty_sites)
        .iter().map(|x| *x as i64).collect();
    let empty_tx = Default::default();
    let t2s_raw = ix.transcripts_to_splices.get(&chr_id).unwrap_or(&empty_tx);
    let mut t2s: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();
    let mut s2t: BTreeMap<Vec<(i64, i64)>, std::collections::BTreeSet<String>> = BTreeMap::new();
    for (tid, sp) in t2s_raw {
        let v: Vec<(i64, i64)> = sp.iter().map(|(a, b)| (*a as i64, *b as i64)).collect();
        s2t.entry(v.clone()).or_default().insert(tid.clone());
        t2s.insert(tid.clone(), v);
    }
    samout::classify_alignment(predicted_splices, &s2t, &t2s, &pairs, &sites)
}

fn to_mems(hits: &[String]) -> BTreeMap<u64, Vec<Mem>> {
    let parsed = reads::get_mems_from_input(hits);
    // `get_mems_from_input` sorts each chromosome's mems by inclusive end and
    // the reference then assigns j by position in that sorted list.
    parsed.into_iter().map(|(c, v)| {
        (c, v.into_iter().enumerate().map(|(j, m)| Mem {
            x: m.x as i64, y: m.y as i64, c: m.c as i64, d: m.d as i64,
            val: m.val as i64, j: j as i64, exon_part_id: m.exon_part_id,
        }).collect())
    }).collect()
}

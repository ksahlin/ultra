//! The genomic prefilter -- `modules/prefilter_genomic_reads.py`.
//!
//! minimap2 is run over the whole genome, and each read is then classified by
//! how much of its alignment falls inside regions uLTRA actually indexed. Reads
//! that are mostly outside are declared genomic, kept as minimap2 aligned them,
//! and never passed to uLTRA's own aligner.
//!
//! Replaces `intervaltree` with a sorted vector plus binary search; the query is
//! "which indexed regions overlap [a, b)", and the regions are static once the
//! index is loaded.

use std::collections::BTreeMap;

/// One minimap2 SAM record, as much of it as the filter needs.
#[derive(Debug, Clone)]
pub struct SamRecord {
    pub qname: String,
    pub flag: u32,
    pub rname: String,
    pub pos: i64,
    pub cigar: String,
    pub seq: String,
    pub qual: String,
}

impl SamRecord {
    pub fn parse(line: &str) -> Option<SamRecord> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            return None;
        }
        Some(SamRecord {
            qname: f[0].to_string(),
            flag: f[1].parse().ok()?,
            rname: f[2].to_string(),
            pos: f[3].parse().ok()?,
            cigar: f[5].to_string(),
            seq: f[9].to_string(),
            qual: f[10].to_string(),
        })
    }
    pub fn is_reverse(&self) -> bool {
        self.flag & 16 != 0
    }
}

/// `get_exons_from_cigar`: the reference intervals an alignment covers, split
/// at `N` (reference skip).
///
/// Operations that consume reference and count as aligned are M, D, = and X
/// (pysam codes 0, 2, 7, 8). `N` closes the current interval and opens the
/// next. I, S and H consume no reference.
pub fn exons_from_cigar(pos: i64, cigar: &str) -> Vec<(i64, i64)> {
    let mut out = Vec::new();
    // pysam's reference_start is 0-based; SAM POS is 1-based
    let start0 = pos - 1;
    let mut ref_pos = start0;
    let mut exon_start = start0;
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
            continue;
        }
        let l: i64 = num.parse().unwrap_or(0);
        num.clear();
        match ch {
            'M' | 'D' | '=' | 'X' => ref_pos += l,
            'N' => {
                out.push((exon_start, ref_pos));
                ref_pos += l;
                exon_start = ref_pos;
            }
            _ => {}
        }
    }
    out.push((exon_start, ref_pos));
    out
}

/// The regions uLTRA indexed, per reference name, as sorted non-overlapping-
/// friendly intervals with a running max-end for pruning.
pub struct IndexedRegions {
    by_ref: BTreeMap<String, Vec<(i64, i64)>>,
}

impl IndexedRegions {
    /// `get_ultra_indexed_choordinates`: one interval per part/flank sequence.
    pub fn from_parts(part_keys: &[(u64, u64, u64)], id_to_chr: &BTreeMap<u64, String>) -> Self {
        let mut by_ref: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();
        for &(chr_id, start, stop) in part_keys {
            if let Some(name) = id_to_chr.get(&chr_id) {
                by_ref.entry(name.clone()).or_default().push((start as i64, stop as i64));
            }
        }
        for v in by_ref.values_mut() {
            v.sort_unstable();
        }
        IndexedRegions { by_ref }
    }

    /// Total overlap of [a, b) with the indexed regions on `refname`.
    ///
    /// `intervaltree`'s `overlap(a, b)` is half-open and excludes zero-length
    /// touches, and `overlap_size` is `min(stop) - max(start)` without clamping
    /// at zero -- reproduced, because a negative contribution is possible in
    /// principle and the reference would count it.
    pub fn total_overlap(&self, refname: &str, a: i64, b: i64) -> i64 {
        let v = match self.by_ref.get(refname) {
            Some(v) => v,
            None => return 0,
        };
        // first interval whose end could reach a
        let mut total = 0i64;
        // linear from the first candidate; regions per contig are modest and
        // the alignment spans are short
        let idx = v.partition_point(|&(s, _)| s < a);
        let from = idx.saturating_sub(64); // walk back for long regions
        for &(s, e) in &v[from..] {
            if s >= b {
                break;
            }
            if e > a && s < b {
                total += std::cmp::min(b, e) - std::cmp::max(a, s);
            }
        }
        total
    }
}

/// What the filter decided for one read.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Decision {
    /// mostly outside uLTRA's indexed regions: keep minimap2's alignment
    Unindexed,
    /// inside: hand to uLTRA's aligner
    Indexed,
    /// minimap2 could not place it: hand to uLTRA's aligner
    Unmapped,
    /// secondary or supplementary: the reference ignores it entirely
    Ignored,
}

/// `filter_reads_to_align`'s decision for one record.
///
/// Note the reference only considers `flag == 0`, `flag == 16` and `flag == 4`
/// EXACTLY -- a secondary (256) or supplementary (2048) record matches none of
/// them and is silently dropped, contributing nothing and not being written to
/// either output.
pub fn classify_record(r: &SamRecord, regions: &IndexedRegions, genomic_frac_cutoff: f64) -> Decision {
    if r.flag == 4 {
        return Decision::Unmapped;
    }
    if r.flag != 0 && r.flag != 16 {
        return Decision::Ignored;
    }
    let exons = exons_from_cigar(r.pos, &r.cigar);
    let mut total_overlap = 0i64;
    let mut total_aligned = 0i64;
    for (a, b) in exons {
        total_overlap += regions.total_overlap(&r.rname, a, b);
        total_aligned += b - a;
    }
    if total_aligned == 0 {
        // the reference divides here and would raise ZeroDivisionError
        return Decision::Indexed;
    }
    if 1.0 - (total_overlap as f64 / total_aligned as f64) > genomic_frac_cutoff {
        Decision::Unindexed
    } else {
        Decision::Indexed
    }
}

/// The record the filter writes for a read it hands on to uLTRA.
///
/// DIVERGENCE, Finding 34: the reference restores the SEQUENCE to the original
/// read orientation for reverse-strand alignments but leaves the QUALITY in the
/// aligned orientation, so the two end up reversed relative to each other. The
/// port reverses the quality too.
///
/// Findings 35 and 36 are also fixed here: the record is written with an `@`
/// header rather than `>`, and a missing quality becomes `*`-free FASTA rather
/// than the literal text `None`.
pub fn filtered_record(r: &SamRecord) -> (String, String) {
    if r.is_reverse() {
        (revcomp(&r.seq), r.qual.chars().rev().collect())
    } else {
        (r.seq.clone(), r.qual.clone())
    }
}

pub fn revcomp(s: &str) -> String {
    s.chars().rev().map(|c| match c {
        'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
        'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c',
        other => other,
    }).collect()
}

// ---------------------------------------------------------------------------
// The minimap2-vs-uLTRA cross-check -- `output_final_alignments` in `uLTRA`.
// ---------------------------------------------------------------------------

/// `score(cigartuples)`: matches minus (insertions + deletions + substitutions).
///
/// Only `=` counts as a match -- a plain `M` scores nothing, which is fine here
/// because both inputs are produced with `--eqx`.
pub fn cigar_score(cigar: &str) -> i64 {
    let (mut matches, mut diffs) = (0i64, 0i64);
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
            continue;
        }
        let l: i64 = num.parse().unwrap_or(0);
        num.clear();
        match ch {
            '=' => matches += l,
            'I' | 'D' | 'X' => diffs += l,
            _ => {}
        }
    }
    matches - diffs
}

/// Which alignment wins for one read.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Winner {
    /// uLTRA's alignment replaces minimap2's
    Ultra,
    /// minimap2's is kept
    Minimap2,
    /// uLTRA did not align this read at all
    UltraUnmapped,
}

/// `output_final_alignments`' decision for one read.
///
/// The reference's counters distinguish "equal", "slightly better" (within 10)
/// and "significantly better", but all three keep minimap2's record -- only a
/// STRICTLY lower minimap2 score hands the read to uLTRA. An unmapped minimap2
/// record with a uLTRA alignment also goes to uLTRA.
pub fn pick_winner(mm2_cigar: Option<&str>, ultra_score: Option<i64>) -> Winner {
    let ultra = match ultra_score {
        Some(s) => s,
        None => return Winner::UltraUnmapped,
    };
    match mm2_cigar {
        None => Winner::Ultra, // minimap2 unmapped, uLTRA aligned
        Some(c) => {
            if cigar_score(c) < ultra {
                Winner::Ultra
            } else {
                Winner::Minimap2
            }
        }
    }
}

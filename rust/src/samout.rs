//! SAM record construction -- `modules/sam_output.py`.
//!
//! This is where the port's bytes finally become `reads.sam`, so everything
//! here is contract down to the field separators.
//!
//! Note the header uLTRA writes elsewhere carries `@SQ` lines ONLY -- no `@PG`,
//! no `@HD` (PORTING.md Part 5), which is what makes SAM a byte-identity target
//! at all.

/// `get_type`: the CIGAR operation for one aligned column.
fn op_for(a: u8, b: u8) -> u8 {
    if a == b { b'=' } else if a == b'-' { b'D' } else if b == b'-' { b'I' } else { b'X' }
}

/// `get_segments`: cut the gapped alignment at the predicted exon boundaries.
///
/// Boundaries are counted in REFERENCE bases consumed, so insertions in the
/// read do not advance them. Consecutive breakpoints created by an insertion
/// sitting exactly at a junction are collapsed, keeping the first.
fn get_segments(read_aln: &[u8], ref_aln: &[u8], predicted_exons: &[(i64, i64)]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut breaks: std::collections::BTreeSet<i64> = Default::default();
    let mut prev = 0i64;
    for (p1, p2) in predicted_exons {
        breaks.insert(p2 - p1 + prev);
        prev += p2 - p1;
    }

    let mut aln_breaks: Vec<usize> = Vec::new();
    let mut cur_ref_pos = 0i64;
    let mut last_i = 0usize;
    for (i, &n) in ref_aln.iter().enumerate() {
        if breaks.contains(&cur_ref_pos) {
            aln_breaks.push(i);
        }
        if n != b'-' {
            cur_ref_pos += 1;
        }
        last_i = i;
    }
    // the reference's trailing check, using the loop's final n and position
    if !ref_aln.is_empty() && ref_aln[last_i] != b'-' && breaks.contains(&cur_ref_pos) {
        aln_breaks.push(last_i);
    }

    // drop consecutive breakpoints, keeping the first of each run
    let mut collapsed: Vec<usize> = Vec::new();
    for (k, &b) in aln_breaks.iter().enumerate() {
        if k == 0 || b != aln_breaks[k - 1] + 1 {
            collapsed.push(b);
        }
    }

    let mut segments = Vec::new();
    let mut e_start = 0usize;
    for (i, &e_stop) in collapsed.iter().enumerate() {
        if i == collapsed.len() - 1 {
            segments.push((read_aln[e_start..].to_vec(), ref_aln[e_start..].to_vec()));
        } else {
            segments.push((read_aln[e_start..e_stop].to_vec(), ref_aln[e_start..e_stop].to_vec()));
        }
        e_start = e_stop;
    }
    segments
}

/// `get_cigars`: run-length encode each segment, trimming the leading and
/// trailing gaps of the whole alignment.
///
/// A leading `D` becomes the reference start offset and is dropped; a leading
/// `I` becomes a soft clip. Same at the end.
fn get_cigars(segments: &[(Vec<u8>, Vec<u8>)]) -> (Vec<String>, i64) {
    let mut out = Vec::new();
    let mut start_offset = 0i64;
    let last = segments.len().saturating_sub(1);
    for (i, (read, rf)) in segments.iter().enumerate() {
        if read.is_empty() {
            out.push(String::new());
            continue;
        }
        let mut types: Vec<u8> = Vec::new();
        let mut lens: Vec<i64> = Vec::new();
        let mut prev = op_for(read[0], rf[0]);
        let mut len = 1i64;
        for k in 1..read.len() {
            let cur = op_for(read[k], rf[k]);
            if cur == prev {
                len += 1;
            } else {
                types.push(prev);
                lens.push(len);
                len = 1;
                prev = cur;
            }
        }
        types.push(prev);
        lens.push(len);

        if i == 0 {
            if types[0] == b'D' {
                start_offset = lens[0];
                lens.remove(0);
                types.remove(0);
            } else if types[0] == b'I' {
                types[0] = b'S';
            }
        }
        if i == last {
            if let Some(&t) = types.last() {
                if t == b'D' {
                    lens.pop();
                    types.pop();
                } else if t == b'I' {
                    let n = types.len() - 1;
                    types[n] = b'S';
                }
            }
        }
        let mut s = String::new();
        for (l, t) in lens.iter().zip(types.iter()) {
            s.push_str(&l.to_string());
            s.push(*t as char);
        }
        out.push(s);
    }
    (out, start_offset)
}

/// `get_genomic_cigar`: join the per-exon CIGARs with `N` runs for the introns.
fn genomic_cigar(read_aln: &[u8], ref_aln: &[u8], predicted_exons: &[(i64, i64)]) -> (String, i64) {
    let segments = get_segments(read_aln, ref_aln, predicted_exons);
    let (cigars, start_offset) = get_cigars(&segments);
    let intron_lengths: Vec<i64> = predicted_exons
        .windows(2)
        .map(|w| w[1].0 - w[0].1)
        .collect();
    let mut s = String::new();
    for (i, c) in cigars.iter().enumerate() {
        s.push_str(c);
        if i < intron_lengths.len() {
            s.push_str(&format!("{}N", intron_lengths[i]));
        }
    }
    (s, start_offset)
}

/// `edit_distance`: X, I and D lengths summed. N and = do not count.
fn edit_distance(cigar: &str) -> i64 {
    let mut ed = 0i64;
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
        } else {
            let n: i64 = num.parse().unwrap_or(0);
            num.clear();
            if ch == 'X' || ch == 'I' || ch == 'D' {
                ed += n;
            }
        }
    }
    ed
}

/// `sam_output.main`.
#[allow(clippy::too_many_arguments)]
pub fn sam_record(
    read_id: &str,
    read_seq: &str,
    read_qual: Option<&str>,
    ref_id: &str,
    classification: &str,
    predicted_exons: &[(i64, i64)],
    read_aln: &str,
    ref_aln: &str,
    annotated_to_transcript_id: &str,
    is_rc: bool,
    is_secondary: bool,
    map_score: i64,
) -> String {
    let aligned = classification != "unaligned";
    let (cigar, flag, reference_start, mapping_quality) = if aligned {
        let (gc, off) = genomic_cigar(read_aln.as_bytes(), ref_aln.as_bytes(), predicted_exons);
        let flag = match (is_secondary, is_rc) {
            (true, true) => 256 + 16,
            (true, false) => 256,
            (false, true) => 16,
            (false, false) => 0,
        };
        // SAM is 1-based, hence the +1
        let start = predicted_exons[0].0 + off + 1;
        (gc, flag, start, map_score)
    } else {
        ("*".to_string(), 4, 0, 0)
    };

    let qual = match read_qual {
        Some(q) if q.len() == read_seq.len() => q,
        _ => "*",
    };

    let mut s = String::with_capacity(read_seq.len() * 2 + 128);
    s.push_str(read_id);
    s.push('\t');
    s.push_str(&flag.to_string());
    s.push('\t');
    s.push_str(ref_id);
    s.push('\t');
    s.push_str(&reference_start.to_string());
    s.push('\t');
    s.push_str(&mapping_quality.to_string());
    s.push('\t');
    s.push_str(&cigar);
    s.push_str("\t*\t0\t0\t");
    s.push_str(read_seq);
    s.push('\t');
    s.push_str(qual);
    s.push('\t');
    if aligned {
        s.push_str(&format!("XA:Z:{annotated_to_transcript_id}\t"));
        s.push_str(&format!("XC:Z:{classification}\t"));
        s.push_str(&format!("NM:i:{}", edit_distance(&cigar)));
    }
    s.push('\n');
    s
}

// ---------------------------------------------------------------------------
// `modules/classify_alignment2.py`
// ---------------------------------------------------------------------------

use std::collections::{BTreeMap, BTreeSet};

/// `contains(sub, pri)`: is `sub` a contiguous subsequence of `pri`?
fn contains_subseq(sub: &[(i64, i64)], pri: &[(i64, i64)]) -> bool {
    if sub.is_empty() || sub.len() > pri.len() {
        return false;
    }
    pri.windows(sub.len()).any(|w| w == sub)
}

/// `classify_alignment2.main` -> (classification, annotated_to_transcript_id)
///
/// NOTE on the FSM branch: the reference joins a **set** of transcript ids with
/// commas, so when several transcripts share a splice pattern the order of the
/// `XA:Z` tag depends on CPython's string-set iteration -- Finding 33's family,
/// and the only member of it that reaches the SAM text directly. The port sorts.
/// Measured on SIRV: 0 of 1 000 FSM results name more than one transcript, so
/// it does not bite there; `tests/sam_oracle.rs` counts the multi-transcript
/// cases it sees so a corpus that does produce them cannot pass silently.
pub fn classify_alignment(
    predicted_splices: &[(i64, i64)],
    splices_to_transcripts: &BTreeMap<Vec<(i64, i64)>, BTreeSet<String>>,
    transcripts_to_splices: &BTreeMap<String, Vec<(i64, i64)>>,
    splice_pairs: &BTreeMap<(i64, i64), BTreeSet<String>>,
    splice_sites: &BTreeSet<i64>,
) -> (&'static str, String) {
    if predicted_splices.is_empty() {
        return ("NO_SPLICE", String::new());
    }

    if let Some(trs) = splices_to_transcripts.get(predicted_splices) {
        // sorted, not set-iteration order -- see above
        let joined = trs.iter().cloned().collect::<Vec<_>>().join(",");
        return ("FSM", joined);
    }

    // NIC: every donor and acceptor is annotated somewhere
    let is_nic = predicted_splices
        .iter()
        .all(|(a, b)| splice_sites.contains(a) && splice_sites.contains(b));
    if is_nic {
        let known_combination = predicted_splices.iter().all(|p| splice_pairs.contains_key(p));
        return if known_combination {
            ("ISM/NIC_known", String::new())
        } else {
            ("NIC_novel", String::new())
        };
    }

    // ISM: every predicted splice is annotated, and they form a contiguous run
    // of some transcript's splices
    let hits: Vec<&BTreeSet<String>> = predicted_splices
        .iter()
        .filter_map(|p| splice_pairs.get(p))
        .collect();
    let in_all: BTreeSet<String> = if hits.is_empty() {
        BTreeSet::new()
    } else {
        let mut acc = hits[0].clone();
        for h in &hits[1..] {
            acc = acc.intersection(h).cloned().collect();
        }
        acc
    };
    for tid in &in_all {
        if let Some(tsp) = transcripts_to_splices.get(tid) {
            if contains_subseq(predicted_splices, tsp) {
                return ("ISM", tid.clone());
            }
        }
    }

    ("NNC", String::new())
}

/// `help_functions.cigar_to_seq`: rebuild the two gapped strings from a CIGAR.
///
/// `=`/`X` consume both, `I` consumes the query and emits a gap in the
/// reference, `D` the reverse. `S` is skipped in the reference's version.
pub fn expand_cigar_full(cigar: &str, query: &[u8], rf: &[u8]) -> (String, String) {
    let (mut q, mut r) = (String::new(), String::new());
    let (mut qi, mut ri) = (0usize, 0usize);
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
            continue;
        }
        let l: usize = num.parse().unwrap_or(0);
        num.clear();
        match ch {
            '=' | 'X' | 'M' => {
                let qe = (qi + l).min(query.len());
                let re = (ri + l).min(rf.len());
                q.push_str(&String::from_utf8_lossy(&query[qi..qe]));
                r.push_str(&String::from_utf8_lossy(&rf[ri..re]));
                qi = qe; ri = re;
            }
            'I' => {
                let qe = (qi + l).min(query.len());
                q.push_str(&String::from_utf8_lossy(&query[qi..qe]));
                r.push_str(&"-".repeat(qe - qi));
                qi = qe;
            }
            'D' => {
                let re = (ri + l).min(rf.len());
                q.push_str(&"-".repeat(re - ri));
                r.push_str(&String::from_utf8_lossy(&rf[ri..re]));
                ri = re;
            }
            _ => {}
        }
    }
    (q, r)
}

/// `help_functions.edlib_alignment`: like the above but for an HW alignment,
/// where the CIGAR covers only `ref[start..=stop]` and the flanks are padded
/// with gaps on the query side.
pub fn expand_cigar(cigar: &str, query: &[u8], rf: &[u8], loc: Option<(i64, i64)>) -> (String, String) {
    let (start, stop) = match loc {
        Some((a, b)) => (a.max(0) as usize, b.max(0) as usize),
        None => (0, rf.len().saturating_sub(1)),
    };
    let end = (stop + 1).min(rf.len());
    let (q, r) = expand_cigar_full(cigar, query, &rf[start.min(rf.len())..end]);
    let tail = rf.len().saturating_sub(stop + 1);
    let qa = format!("{}{}{}", "-".repeat(start), q, "-".repeat(tail));
    let ra = format!("{}{}{}",
        String::from_utf8_lossy(&rf[..start.min(rf.len())]), r,
        String::from_utf8_lossy(&rf[end.min(rf.len())..]));
    (qa, ra)
}

//! Read parsing and preprocessing, matching `modules/help_functions.py`.
//!
//! `readfq` is Heng Li's generator, and the port reproduces its quirks rather
//! than improving on them, because they are observable:
//!
//!   * the accession is `last[1:].split()[0]` -- the first whitespace token.
//!     Unlike isONclust's copy, this one does NOT substitute spaces.
//!   * every line is taken as `l[:-1]`, i.e. the final character is dropped
//!     unconditionally. On a file whose last line has no trailing newline the
//!     last base (or quality value) is silently lost.
//!   * quality is accumulated until `leng >= len(seq)`; if EOF arrives first
//!     the record degrades to a fasta record with no quality.

use std::io::BufRead;

#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub name: String,
    pub seq: String,
    pub qual: Option<String>,
}

#[derive(Debug)]
pub enum ReadError {
    Io(std::io::Error),
    /// Finding 6: the reference cannot see this and silently drops reads.
    QualityLengthMismatch { name: String, seq_len: usize, qual_len: usize },
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadError::Io(e) => write!(f, "{e}"),
            ReadError::QualityLengthMismatch { name, seq_len, qual_len } => write!(
                f,
                "malformed fastq: read '{name}' has {seq_len} bases but {qual_len} quality values"
            ),
        }
    }
}

/// Parse a fasta/fastq stream the way `readfq` does.
///
/// `strict` turns Finding 6 into an error instead of a silent truncation. The
/// reference has no such mode: a quality line of the wrong length desynchronises
/// its parser, it yields fewer records, and `pc.file_IO`'s `zip()` then drops
/// every read after the first bad one without a word. See PORTING.md Finding 6.
pub fn readfq<R: BufRead>(r: R, strict: bool) -> Result<Vec<Record>, ReadError> {
    let mut out = Vec::new();
    let mut last: Option<String> = None;
    // `l[:-1]` drops the final character of every line, newline or not.
    let lines: Vec<String> = {
        let mut v = Vec::new();
        for l in r.lines() {
            v.push(l.map_err(ReadError::Io)?);
        }
        v
    };
    let mut i = 0usize;

    loop {
        if last.is_none() {
            while i < lines.len() {
                let l = &lines[i];
                i += 1;
                if l.starts_with('>') || l.starts_with('@') {
                    last = Some(l.clone());
                    break;
                }
            }
        }
        let header = match last.take() {
            Some(h) => h,
            None => break,
        };
        let name = header[1..].split_whitespace().next().unwrap_or("").to_string();

        let mut seqs: Vec<String> = Vec::new();
        let mut sep: Option<String> = None;
        while i < lines.len() {
            let l = &lines[i];
            i += 1;
            if l.starts_with('@') || l.starts_with('+') || l.starts_with('>') {
                sep = Some(l.clone());
                break;
            }
            seqs.push(l.clone());
        }

        let seq = seqs.concat();
        match sep {
            Some(ref sl) if sl.starts_with('+') => {
                // fastq: accumulate quality until it is at least as long as seq
                let mut quals: Vec<String> = Vec::new();
                let mut leng = 0usize;
                let mut done = false;
                while i < lines.len() {
                    let l = &lines[i];
                    i += 1;
                    quals.push(l.clone());
                    leng += l.len();
                    if leng >= seq.len() {
                        last = None;
                        done = true;
                        break;
                    }
                }
                let qual = quals.concat();
                if done {
                    if strict && qual.len() != seq.len() {
                        return Err(ReadError::QualityLengthMismatch {
                            name,
                            seq_len: seq.len(),
                            qual_len: qual.len(),
                        });
                    }
                    out.push(Record { name, seq, qual: Some(qual) });
                } else {
                    // EOF before enough quality: the reference degrades to fasta
                    if strict {
                        return Err(ReadError::QualityLengthMismatch {
                            name,
                            seq_len: seq.len(),
                            qual_len: qual.len(),
                        });
                    }
                    out.push(Record { name, seq, qual: None });
                    break;
                }
            }
            other => {
                out.push(Record { name, seq, qual: None });
                last = other;
                if last.is_none() {
                    break;
                }
            }
        }
    }
    Ok(out)
}

/// `help_functions.remove_read_polyA_ends`.
///
/// Compresses homopolymer runs of A or T longer than `threshold_len` to
/// `to_len`, but only within the last `min(len(seq)//2, 100)` bases.
///
/// EDGE CASE reproduced deliberately: when that window is 0 -- which happens
/// for a sequence shorter than 2 bases -- the reference computes `seq[:-0]`,
/// and in Python `-0 == 0`, so the prefix is `seq[:0]`, the empty string. The
/// whole read is dropped. Faithful, and noted here so it is not mistaken for a
/// porting slip.
pub fn remove_read_polya_ends(seq: &str, threshold_len: usize, to_len: usize) -> String {
    remove_read_polya_ends_q(seq, None, threshold_len, to_len).0
}

/// The same, carrying the quality string along.
///
/// The reference compresses BOTH together -- `qual_list.extend(qualgroup[:to_len])`
/// keeps the FIRST `to_len` quality values of a collapsed homopolymer run, so
/// the quality stays the same length as the sequence. Dropping this is what
/// made the port emit `*` for the quality on every read with a polyA tail.
pub fn remove_read_polya_ends_q(
    seq: &str, qual: Option<&str>, threshold_len: usize, to_len: usize,
) -> (String, Option<String>) {
    let n = seq.len();
    let window = std::cmp::min(n / 2, 100);
    if window == 0 {
        return (String::new(), qual.map(|_| String::new()));
    }
    let head = &seq[..n - window];
    let tail = seq[n - window..].as_bytes();
    let qtail: Option<&[u8]> = qual.map(|q| {
        let qn = q.len();
        if qn >= window { &q.as_bytes()[qn - window..] } else { q.as_bytes() }
    });

    let mut out = String::with_capacity(n);
    out.push_str(head);
    let mut qout: Option<String> = qual.map(|q| {
        let qn = q.len();
        if qn >= window { q[..qn - window].to_string() } else { String::new() }
    });

    let mut i = 0usize;
    while i < tail.len() {
        let ch = tail[i];
        let mut j = i;
        while j < tail.len() && tail[j] == ch { j += 1; }
        let run = j - i;
        let keep = if run > threshold_len && (ch == b'A' || ch == b'T') { to_len } else { run };
        for _ in 0..keep { out.push(ch as char); }
        if let (Some(qo), Some(qt)) = (qout.as_mut(), qtail) {
            // the FIRST `keep` quality values of this run
            for k in 0..keep.min(run) {
                if i + k < qt.len() { qo.push(qt[i + k] as char); }
            }
        }
        i = j;
    }
    (out, qout)
}

#[allow(dead_code)]
fn remove_read_polya_ends_seq_only(seq: &str, threshold_len: usize, to_len: usize) -> String {
    let n = seq.len();
    let window = std::cmp::min(n / 2, 100);
    if window == 0 {
        return String::new();
    }
    let head = &seq[..n - window];
    let tail = &seq[n - window..];

    let mut out = String::with_capacity(n);
    out.push_str(head);
    let bytes = tail.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        let ch = bytes[i];
        let mut j = i;
        while j < bytes.len() && bytes[j] == ch {
            j += 1;
        }
        let run = j - i;
        let keep = if run > threshold_len && (ch == b'A' || ch == b'T') { to_len } else { run };
        for _ in 0..keep {
            out.push(ch as char);
        }
        i = j;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    fn parse(s: &str) -> Vec<Record> {
        readfq(BufReader::new(s.as_bytes()), false).unwrap()
    }

    #[test]
    fn fasta_and_fastq() {
        let r = parse(">a desc\nACGT\n>b\nTTTT\n");
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].name, "a");          // first token only
        assert_eq!(r[0].seq, "ACGT");
        assert_eq!(r[0].qual, None);
        let q = parse("@x\nACGT\n+\nIIII\n");
        assert_eq!(q[0].seq, "ACGT");
        assert_eq!(q[0].qual.as_deref(), Some("IIII"));
    }

    #[test]
    fn strict_rejects_what_the_reference_silently_drops() {
        // Finding 6: three quality values for four bases.
        let e = readfq(BufReader::new("@x\nACGT\n+\nIII\n".as_bytes()), true);
        assert!(matches!(e, Err(ReadError::QualityLengthMismatch { .. })));
    }

    #[test]
    fn polya_window_and_threshold() {
        // 120 bases: window is min(60,100)=60, so only the last 60 are touched.
        let s = format!("{}{}", "C".repeat(60), "A".repeat(60));
        let got = remove_read_polya_ends(&s, 8, 5);
        assert_eq!(got, format!("{}{}", "C".repeat(60), "A".repeat(5)));
        // a run of exactly threshold_len is NOT compressed (strict >)
        let s2 = format!("{}{}{}", "C".repeat(100), "A".repeat(8), "C".repeat(92));
        assert_eq!(remove_read_polya_ends(&s2, 8, 5), s2);
    }

    #[test]
    fn polya_edge_case_drops_short_reads() {
        // window == 0 -> seq[:-0] == seq[:0] == "" in Python. Faithful.
        assert_eq!(remove_read_polya_ends("A", 8, 5), "");
    }
}

// ---------------------------------------------------------------------------
// namfinder output parsing: `seed_wrapper.read_seeds` + `align.get_mems_from_input`.
// ---------------------------------------------------------------------------

/// One NAM as uLTRA interprets it. Coordinates are converted exactly as
/// `get_mems_from_input` does: the genome start is the part's own start plus
/// namfinder's 1-based offset minus one, and the END IS INCLUSIVE, because MEM
/// solvers report it that way and uLTRA keeps that convention.
#[derive(Debug, Clone, PartialEq)]
pub struct Mem {
    pub chr_id: u64,
    /// inclusive genome interval
    pub x: u64,
    pub y: u64,
    /// inclusive read interval, 0-based
    pub c: u64,
    pub d: u64,
    pub val: u64,
    /// the `chr^start^stop` token, kept verbatim as the reference does
    pub exon_part_id: String,
}

/// A read's NAMs, forward and reverse-complement.
#[derive(Debug, Clone)]
pub struct SeedRecord {
    pub acc: String,
    pub hits: Vec<String>,
    pub acc_rev: String,
    pub hits_rc: Vec<String>,
}

/// `read_seeds`: split namfinder's output into per-read forward/reverse blocks.
///
/// FAITHFUL QUIRK: the reference only yields when BOTH `curr_acc` and
/// `curr_acc_rev` are set, and it resets both after each yield. A read whose
/// Reverse block never appears is therefore dropped silently -- which is one
/// half of Finding 6's truncation, the other half being `zip()`.
pub fn read_seeds<R: BufRead>(r: R) -> Result<Vec<SeedRecord>, ReadError> {
    let mut out = Vec::new();
    let (mut acc, mut acc_rev) = (String::new(), String::new());
    let (mut hits, mut hits_rc): (Vec<String>, Vec<String>) = (Vec::new(), Vec::new());
    let mut is_rc = false;

    for line in r.lines() {
        let line = line.map_err(ReadError::Io)?;
        if line.starts_with('>') {
            if !acc.is_empty() && !acc_rev.is_empty() {
                out.push(SeedRecord {
                    acc: std::mem::take(&mut acc),
                    hits: std::mem::take(&mut hits),
                    acc_rev: std::mem::take(&mut acc_rev),
                    hits_rc: std::mem::take(&mut hits_rc),
                });
                is_rc = false;
            }
            let name = line[1..].trim().to_string();
            if name.contains("Reverse") {
                acc_rev = name;
                is_rc = true;
            } else {
                acc = name;
                is_rc = false;
            }
        } else if is_rc {
            hits_rc.push(line);
        } else {
            hits.push(line);
        }
    }
    if !acc.is_empty() && !acc_rev.is_empty() {
        out.push(SeedRecord { acc, hits, acc_rev, hits_rc });
    }
    Ok(out)
}

/// `get_mems_from_input`: parse hit lines into Mems, grouped by chr and sorted
/// by the inclusive genome END (`x[1]`), which is the reference's sort key.
///
/// The sort is STABLE and keyed on the end alone, so entries tied on it keep
/// their input order -- that is the mechanism behind Finding 27, and it is
/// reproduced rather than tidied.
pub fn get_mems_from_input(hits: &[String]) -> std::collections::BTreeMap<u64, Vec<Mem>> {
    let mut by_chr: std::collections::BTreeMap<u64, Vec<Mem>> = Default::default();
    for line in hits {
        let v: Vec<&str> = line.split_whitespace().collect();
        if v.len() < 4 {
            continue;
        }
        let part = v[0];
        let mut it = part.split('^');
        let (c, ps) = match (it.next(), it.next()) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };
        let chr_id: u64 = match c.parse() { Ok(x) => x, Err(_) => continue };
        let part_start: u64 = match ps.parse() { Ok(x) => x, Err(_) => continue };
        let ref_off: u64 = match v[1].parse::<u64>() { Ok(x) => x - 1, Err(_) => continue };
        let read_start: u64 = match v[2].parse::<u64>() { Ok(x) => x - 1, Err(_) => continue };
        let len: u64 = match v[3].parse() { Ok(x) => x, Err(_) => continue };
        let g = part_start + ref_off;
        by_chr.entry(chr_id).or_default().push(Mem {
            chr_id,
            x: g,
            y: g + len - 1,
            c: read_start,
            d: read_start + len - 1,
            val: len,
            exon_part_id: part.to_string(),
        });
    }
    for v in by_chr.values_mut() {
        v.sort_by_key(|m| m.y); // stable, end only -- see above
    }
    by_chr
}

#[cfg(test)]
mod seed_tests {
    use super::*;
    use std::io::BufReader;

    const SAMPLE: &str = "> read1\n  6^2285^2620 93 151 125\n  6^2740^2828 4 398 81\n> read1 Reverse\n> read4\n> read4 Reverse\n  3^1532^1764 144 305 89\n";

    #[test]
    fn splits_forward_and_reverse_blocks() {
        let r = read_seeds(BufReader::new(SAMPLE.as_bytes())).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].acc, "read1");
        assert_eq!(r[0].hits.len(), 2);
        assert_eq!(r[0].hits_rc.len(), 0);
        assert_eq!(r[1].acc, "read4");
        assert_eq!(r[1].hits.len(), 0);
        assert_eq!(r[1].hits_rc.len(), 1);
    }

    #[test]
    fn mem_coordinates_match_the_reference_arithmetic() {
        // "6^2285^2620 93 151 125": part starts at 2285, namfinder's ref offset
        // is 1-based 93, read start 1-based 151, length 125.
        let m = get_mems_from_input(&["  6^2285^2620 93 151 125".to_string()]);
        let v = &m[&6];
        assert_eq!(v[0].x, 2285 + 92);
        assert_eq!(v[0].y, 2285 + 92 + 125 - 1); // END IS INCLUSIVE
        assert_eq!(v[0].c, 150);
        assert_eq!(v[0].d, 150 + 125 - 1);
        assert_eq!(v[0].val, 125);
        assert_eq!(v[0].exon_part_id, "6^2285^2620");
    }

    #[test]
    fn sorted_by_inclusive_end_per_chromosome() {
        let hits: Vec<String> = SAMPLE
            .lines()
            .filter(|l| !l.starts_with('>'))
            .map(|s| s.to_string())
            .collect();
        let m = get_mems_from_input(&hits);
        for v in m.values() {
            assert!(v.windows(2).all(|w| w[0].y <= w[1].y));
        }
    }
}

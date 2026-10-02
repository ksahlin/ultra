//! The minimap2 pre-pass and the final merge -- `prefilter_genomic_reads.py`
//! and `output_final_alignments` in the `uLTRA` script.
//!
//! The decision logic (what counts as genomic, which alignment wins) lives in
//! `prefilter.rs` and is pinned by `tests/prefilter_oracle.rs`. This module is
//! the I/O around it: running minimap2, splitting its SAM, and merging.

use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::prefilter::{self, Decision, IndexedRegions, SamRecord};
use crate::samfmt::format_g;

type R<T> = std::io::Result<T>;

fn err(msg: String) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::Other, msg)
}

/// `align_with_minimap2`: `minimap2 -ax splice --eqx -k K -t T ref reads`.
///
/// stdout becomes `minimap2.sam`, stderr `minimap2_errors.1`, exactly as the
/// reference names them.
pub fn run_minimap2(
    ref_path: &str, reads_path: &str, outdir: &Path, nr_cores: i64, k_size: i64,
) -> R<PathBuf> {
    let sam = outdir.join("minimap2.sam");
    let errs = outdir.join("minimap2_errors.1");
    let out_f = std::fs::File::create(&sam)?;
    let err_f = std::fs::File::create(&errs)?;

    let status = std::process::Command::new("minimap2")
        .args(["-ax", "splice", "--eqx", "-k"])
        .arg(k_size.to_string())
        .arg("-t")
        .arg(nr_cores.to_string())
        .arg(ref_path)
        .arg(reads_path)
        .stdout(std::process::Stdio::from(out_f))
        .stderr(std::process::Stdio::from(err_f))
        .status()
        .map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                err("minimap2 was not found on PATH. uLTRA uses it to detect reads \
                     aligning outside the annotated regions. Install it, or pass \
                     --disable_mm2 to skip the genomic prefilter."
                    .to_string())
            } else {
                e
            }
        })?;
    if !status.success() {
        return Err(err(format!(
            "minimap2 failed with status {status}; see {}",
            errs.display()
        )));
    }
    Ok(sam)
}

/// What one read looked like to the prefilter, for the FASTQ it writes.
///
/// DIVERGENCE, Findings 35 and 36: the reference writes every record as
/// `>acc\nseq\n+\nqual\n` -- a `>` header over a four-line body, and the
/// literal text `None` when the input had no qualities. The port writes a
/// well-formed record instead: FASTQ when there are qualities, FASTA when
/// there are not. See Finding 39 for what the malformed file costs.
fn write_read(w: &mut impl Write, acc: &str, seq: &str, qual: Option<&str>) -> R<()> {
    match qual {
        Some(q) => writeln!(w, "@{acc}\n{seq}\n+\n{q}"),
        None => writeln!(w, ">{acc}\n{seq}"),
    }
}

/// Rewrite a minimap2 record the way pysam does when it copies one: every
/// `f`-typed tag is reformatted. Everything else is passed through untouched
/// -- verified field by field against 10 000 round-tripped records.
fn pysam_rewrite(line: &str) -> String {
    if !line.contains(":f:") {
        return line.to_string();
    }
    line.split('\t')
        .map(|f| {
            let b = f.as_bytes();
            if b.len() > 5 && &f[2..5] == ":f:" {
                if let Ok(v) = f[5..].parse::<f64>() {
                    return format!("{}:f:{}", &f[..2], format_g(v));
                }
            }
            f.to_string()
        })
        .collect::<Vec<_>>()
        .join("\t")
}

pub struct FilterOutcome {
    pub nr_unindexed: usize,
    pub reads_path: PathBuf,
    pub indexed_path: PathBuf,
    pub unindexed_path: PathBuf,
}

/// `filter_reads_to_align`: split minimap2's SAM into the reads uLTRA should
/// align and the ones it should not touch.
pub fn filter_reads_to_align(
    mm2_sam: &Path, regions: &IndexedRegions, outdir: &Path, genomic_frac: f64,
) -> R<FilterOutcome> {
    let reads_path = outdir.join("reads_after_genomic_filtering.fastq");
    let indexed_path = outdir.join("indexed.sam");
    let unindexed_path = outdir.join("unindexed.sam");

    let f = std::fs::File::open(mm2_sam)?;
    let rdr = std::io::BufReader::new(f);
    let mut reads_w = BufWriter::new(std::fs::File::create(&reads_path)?);
    let mut idx_w = BufWriter::new(std::fs::File::create(&indexed_path)?);
    let mut unidx_w = BufWriter::new(std::fs::File::create(&unindexed_path)?);

    let mut in_header = true;
    let mut nr_unindexed = 0usize;

    for line in rdr.lines() {
        let line = line?;
        if in_header {
            if line.starts_with('@') {
                // both split files are written with minimap2's header, as
                // pysam's `template=SAM_file` does
                writeln!(idx_w, "{line}")?;
                writeln!(unidx_w, "{line}")?;
                continue;
            }
            in_header = false;
        }
        let rec = match SamRecord::parse(&line) {
            Some(r) => r,
            None => continue,
        };
        let qual = if rec.qual == "*" { None } else { Some(rec.qual.as_str()) };

        match prefilter::classify_record(&rec, regions, genomic_frac) {
            Decision::Unindexed => {
                // pysam set_tag appends a tag that is not already present
                writeln!(unidx_w, "{}\tXA:Z:\tXC:Z:uLTRA_unindexed", pysam_rewrite(&line))?;
                nr_unindexed += 1;
            }
            Decision::Indexed => {
                let (seq, q) = prefilter::filtered_record(&rec);
                let q = qual.map(|_| q);
                write_read(&mut reads_w, &rec.qname, &seq, q.as_deref())?;
                writeln!(idx_w, "{}", pysam_rewrite(&line))?;
            }
            Decision::Unmapped => {
                // minimap2 could not place it; hand the read on untouched
                write_read(&mut reads_w, &rec.qname, &rec.seq, qual)?;
            }
            Decision::Ignored => {}
        }
    }
    reads_w.flush()?;
    idx_w.flush()?;
    unidx_w.flush()?;
    Ok(FilterOutcome { nr_unindexed, reads_path, indexed_path, unindexed_path })
}

// ---------------------------------------------------------------------------
// The merge
// ---------------------------------------------------------------------------

#[derive(Default, Debug)]
pub struct MergeStats {
    pub ultra_better: usize,
    pub equal_score: usize,
    pub slightly_worse: usize,
    pub worse: usize,
    pub ultra_unmapped: usize,
    pub not_attempted: usize,
    /// Finding 39: READS minimap2 never produced a record for (not records --
    /// uLTRA can emit more than one alignment per read). The reference drops
    /// these; the port keeps uLTRA's alignments.
    pub mm2_absent: usize,
}

fn qname_of(line: &str) -> &str {
    line.split('\t').next().unwrap_or("")
}

fn field(line: &str, n: usize) -> &str {
    line.split('\t').nth(n).unwrap_or("")
}

fn flag_of(line: &str) -> u32 {
    field(line, 1).parse().unwrap_or(0)
}

fn is_secondary(line: &str) -> bool {
    flag_of(line) & 256 != 0
}

fn is_unmapped(line: &str) -> bool {
    flag_of(line) & 4 != 0
}

/// Run `f` over every non-header line of a SAM file.
fn for_each_record(path: &Path, mut f: impl FnMut(&str) -> R<()>) -> R<()> {
    let fh = std::fs::File::open(path)?;
    for line in std::io::BufReader::new(fh).lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }
        if line.is_empty() {
            continue;
        }
        f(&line)?;
    }
    Ok(())
}

/// `output_final_alignments`: keep minimap2's record for each read unless
/// uLTRA's scores strictly better.
///
/// Follows the reference's multi-pass shape rather than holding both SAMs in
/// memory -- the reference's own comment says it chose passes over RAM, and
/// that trade is right for the port too.
pub fn output_final_alignments(
    ultra_sam: &Path, indexed: &Path, unindexed: &Path, header: &[String],
) -> R<MergeStats> {
    use std::collections::{HashMap, HashSet};
    let mut st = MergeStats::default();

    // Pass 1: uLTRA's score per read, primary and mapped only.
    let mut ultra_scores: HashMap<String, i64> = HashMap::new();
    for_each_record(ultra_sam, |l| {
        if !is_secondary(l) && !is_unmapped(l) {
            ultra_scores.insert(qname_of(l).to_string(), prefilter::cigar_score(field(l, 5)));
        }
        Ok(())
    })?;

    // Pass 2: which reads uLTRA wins, and which reads minimap2 covered at all.
    let mut ultra_better: HashSet<String> = HashSet::new();
    let mut mm2_seen: HashSet<String> = HashSet::new();
    for_each_record(indexed, |l| {
        let q = qname_of(l);
        mm2_seen.insert(q.to_string());
        let ultra = match ultra_scores.get(q) {
            Some(s) => *s,
            None => {
                st.ultra_unmapped += 1;
                return Ok(());
            }
        };
        if !is_secondary(l) {
            let mm2 = prefilter::cigar_score(field(l, 5));
            if mm2 < ultra {
                ultra_better.insert(q.to_string());
            } else if mm2 == ultra {
                st.equal_score += 1;
            } else if mm2 <= ultra + 10 {
                st.slightly_worse += 1;
            } else {
                st.worse += 1;
            }
            ultra_scores.remove(q);
        }
        Ok(())
    })?;
    drop(ultra_scores);
    for_each_record(unindexed, |l| {
        mm2_seen.insert(qname_of(l).to_string());
        Ok(())
    })?;

    // Pass 3: the full uLTRA records for the reads it wins.
    let mut winners: HashMap<String, String> = HashMap::new();
    for_each_record(ultra_sam, |l| {
        let q = qname_of(l);
        if !is_secondary(l) && ultra_better.contains(q) {
            winners.insert(q.to_string(), l.to_string());
        }
        Ok(())
    })?;
    st.ultra_better = winners.len();

    // Pass 4: write the merged file -- uLTRA's header, minimap2's record order.
    let tmp = ultra_sam.with_extension("samtmp");
    {
        let mut w = BufWriter::new(std::fs::File::create(&tmp)?);
        for h in header {
            writeln!(w, "{h}")?;
        }
        for_each_record(indexed, |l| {
            if !is_secondary(l) {
                if let Some(u) = winners.get(qname_of(l)) {
                    return writeln!(w, "{u}");
                }
            }
            writeln!(w, "{l}")
        })?;
        for_each_record(unindexed, |l| {
            if !is_secondary(l) {
                st.not_attempted += 1;
            }
            writeln!(w, "{l}")
        })?;

        // Pass 5 -- DIVERGENCE, Finding 39. A read minimap2 produced no record
        // for appears in neither split file, so the reference's merge has no
        // slot for it and uLTRA's alignment is discarded. The author wrote the
        // branch for this case (`replaced_unaligned_cntr`) but it can never
        // fire. The port emits uLTRA's records instead of dropping the read.
        let mut kept: HashSet<String> = HashSet::new();
        for_each_record(ultra_sam, |l| {
            let q = qname_of(l);
            if !mm2_seen.contains(q) {
                kept.insert(q.to_string());
                return writeln!(w, "{l}");
            }
            Ok(())
        })?;
        st.mm2_absent = kept.len();
        w.flush()?;
    }
    std::fs::rename(&tmp, ultra_sam)?;
    Ok(st)
}

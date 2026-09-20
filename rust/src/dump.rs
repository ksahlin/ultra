//! Render an in-memory `Index` in exactly the canonical text form that
//! `bench/dump_reference.py` produces from the reference's pickles.
//!
//! This is the index stage's contract. The two renderers must agree
//! character for character, so the formatting rules live in one place in each
//! and are stated identically:
//!
//!   * one file per structure, `<name>.txt`
//!   * mappings: `key<TAB>value`, sorted by the rendered key (byte-wise)
//!   * sets: one entry per line, sorted
//!   * a (chr, start, stop) key renders as `chr/start-stop`
//!   * a list of such keys renders comma-separated, SORTED (Finding 21: the
//!     reference's order is hash-seed noise and both consumers de-duplicate)
//!   * every file ends with a trailing newline if it has any lines

use crate::index::{Index, Key};
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::io::Write as _;

fn k(key: &Key) -> String {
    format!("{}/{}-{}", key.0, key.1, key.2)
}

fn join_keys_sorted(v: &[Key]) -> String {
    let mut s: Vec<&Key> = v.iter().collect();
    s.sort();
    s.iter().map(|x| k(x)).collect::<Vec<_>>().join(",")
}

fn write_lines(dir: &std::path::Path, name: &str, lines: &[String]) -> std::io::Result<()> {
    let mut f = std::io::BufWriter::new(std::fs::File::create(dir.join(format!("{name}.txt")))?);
    for l in lines {
        f.write_all(l.as_bytes())?;
        f.write_all(b"\n")?;
    }
    f.flush()
}

/// Sort rendered lines the way Python's `sort()` does on `str`: by Unicode
/// code point, which for our ASCII keys is byte order.
fn sorted_lines(mut v: Vec<String>) -> Vec<String> {
    v.sort();
    v
}

pub fn dump(ix: &Index, out: &std::path::Path) -> std::io::Result<Vec<(String, usize)>> {
    std::fs::create_dir_all(out)?;
    let mut written = Vec::new();

    macro_rules! emit {
        ($name:expr, $lines:expr) => {{
            let lines: Vec<String> = $lines;
            write_lines(out, $name, &lines)?;
            written.push(($name.to_string(), lines.len()));
        }};
    }

    emit!("chr_to_id", sorted_lines(ix.chr_to_id.iter().map(|(n, i)| format!("{n}\t{i}")).collect()));
    emit!("id_to_chr", sorted_lines(ix.id_to_chr.iter().map(|(i, n)| format!("{i}\t{n}")).collect()));
    emit!("refs_lengths", sorted_lines(ix.refs_lengths.iter().map(|(n, l)| format!("{n}\t{l}")).collect()));
    emit!("refs_id_lengths", sorted_lines(ix.refs_id_lengths.iter().map(|(i, l)| format!("{i}\t{l}")).collect()));
    emit!("max_intron_chr", sorted_lines(ix.max_intron_chr.iter().map(|(c, m)| format!("{c}\t{m}")).collect()));

    emit!("segment_id_to_choordinates", sorted_lines(
        ix.segment_id_to_choordinates.iter().map(|(key, (a, b))| format!("{}\t{}/{}", k(key), a, b)).collect()));
    emit!("segment_to_ref", sorted_lines(
        ix.segment_to_ref.iter().map(|(key, c)| format!("{}\t{}", k(key), c)).collect()));
    emit!("segment_to_gene", sorted_lines(
        ix.segment_to_gene.iter().map(|(key, g)| {
            format!("{}\t{}", k(key), g.iter().cloned().collect::<Vec<_>>().join(","))
        }).collect()));
    emit!("parts_to_segments", sorted_lines(
        ix.parts_to_segments.iter().map(|(key, v)| format!("{}\t{}", k(key), join_keys_sorted(v))).collect()));
    emit!("gene_to_small_segments", sorted_lines(
        ix.gene_to_small_segments.iter().map(|(g, v)| format!("{}\t{}", g, join_keys_sorted(v))).collect()));

    emit!("exon_choordinates_to_id", sorted_lines(ix.exon_choordinates_to_id.iter().map(k).collect()));
    emit!("flank_choordinates", sorted_lines(ix.flank_choordinates.iter().map(k).collect()));

    // nested: chr -> {(site1,site2): {transcript,...}}
    emit!("all_splice_pairs_annotations", sorted_lines(
        ix.all_splice_pairs_annotations.iter().map(|(c, m)| {
            let mut parts: Vec<String> = m.iter().map(|((a, b), ts)| {
                format!("({a},{b})={}", ts.iter().cloned().collect::<Vec<_>>().join(","))
            }).collect();
            parts.sort();
            format!("{c}\t{}", parts.join(";"))
        }).collect()));

    emit!("all_splice_sites_annotations", sorted_lines(
        ix.all_splice_sites_annotations.iter().map(|(c, s)| {
            let mut v: Vec<String> = s.iter().map(|x| x.to_string()).collect();
            v.sort();
            format!("{c}\t{}", v.join(","))
        }).collect()));

    // chr -> {(splice,...): {transcript,...}}
    emit!("splices_to_transcripts", sorted_lines(
        ix.splices_to_transcripts.iter().map(|(c, m)| {
            let mut parts: Vec<String> = m.iter().map(|(sp, ts)| {
                let key = render_splice_tuple(sp);
                format!("{key}={}", ts.iter().cloned().collect::<Vec<_>>().join(","))
            }).collect();
            parts.sort();
            format!("{c}\t{}", parts.join(";"))
        }).collect()));

    emit!("transcripts_to_splices", sorted_lines(
        ix.transcripts_to_splices.iter().map(|(c, m)| {
            let mut parts: Vec<String> = m.iter().map(|(tid, sp)| {
                format!("{tid}={}", render_splice_tuple_value(sp))
            }).collect();
            parts.sort();
            format!("{c}\t{}", parts.join(";"))
        }).collect()));

    emit!("ref_part_sequences", seq_lines(&ix.ref_part_sequences));
    emit!("ref_segment_sequences", seq_lines(&ix.ref_segment_sequences));
    emit!("ref_exon_sequences", seq_lines(&ix.ref_exon_sequences));
    emit!("ref_flank_sequences", seq_lines(&ix.ref_flank_sequences));

    Ok(written)
}

fn seq_lines(m: &std::collections::BTreeMap<Key, String>) -> Vec<String> {
    sorted_lines(m.iter().map(|(key, s)| format!("{}\t{}", k(key), s)).collect())
}

/// A tuple-of-pairs used as a dict KEY renders via `_atom`, i.e.
/// "(a,b)" joined by "," inside an outer "(...)": `((1,2),(3,4))`.
/// An EMPTY tuple renders as "()".
fn render_splice_tuple(sp: &[(u64, u64)]) -> String {
    let mut s = String::from("(");
    let mut first = true;
    for (a, b) in sp {
        if !first {
            s.push(',');
        }
        let _ = write!(s, "({a},{b})");
        first = false;
    }
    s.push(')');
    s
}

/// A tuple-of-pairs used as a dict VALUE goes through `render_value`, which
/// joins with "/" rather than wrapping in parens.
fn render_splice_tuple_value(sp: &[(u64, u64)]) -> String {
    sp.iter().map(|(a, b)| format!("({a},{b})")).collect::<Vec<_>>().join("/")
}

#[allow(dead_code)]
fn unused(_: &BTreeSet<u64>) {}

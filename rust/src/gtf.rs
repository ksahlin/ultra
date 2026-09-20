//! GTF parsing, replacing gffutils and its sqlite database.
//!
//! The reference builds a gffutils database -- 140.6 s and 3 719 MB for GENCODE
//! v47, 240.3 s for Ensembl Drosophila on the default path -- and then issues
//! exactly three queries against it:
//!
//!     db.features_of_type('exon', order_by='seqid')          <- 'seqid,start' after Finding 13
//!     db.features_of_type('transcript', order_by='seqid')
//!     db.children(transcript, featuretype='exon', order_by='start')
//!
//! All three are satisfied by parsing the file once and sorting. See PORTING.md
//! Part 3, Hypothesis 1.
//!
//! ORDERING. `order_by='seqid,start'` is SQL `ORDER BY seqid, start` where seqid
//! is TEXT under SQLite's default BINARY collation and start is INTEGER. Rust's
//! `str` Ord is byte-wise and `u64` Ord is numeric, so the comparison matches;
//! ties are resolved by gffutils in rowid (file) order, which a STABLE sort over
//! records held in file order reproduces.

use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Exon {
    /// Identity token only. gffutils gives every exon LINE its own id and does
    /// not deduplicate by coordinate (measured: 339 ids for 186 distinct
    /// coordinate pairs on SIRV, max multiplicity 9), so a running counter is
    /// equivalent.
    pub id: u32,
    pub seqid: String,
    /// 0-based, i.e. the reference's `exon.start - 1`.
    pub start: u64,
    /// 1-based inclusive end == 0-based exclusive end, i.e. `exon.stop`.
    pub stop: u64,
    pub gene_ids: Vec<String>,
    pub transcript_id: Option<String>,
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub seqid: String,
    /// (start0, stop) of each exon of this transcript, sorted by start.
    pub exons: Vec<(u64, u64)>,
}

pub struct Gtf {
    pub exons: Vec<Exon>,
    pub transcripts: Vec<Transcript>,
}

#[derive(Debug)]
pub enum GtfError {
    Io(std::io::Error),
    /// end < start, which is what issue #20 reports from a downstream tool.
    BadCoordinates { line_no: usize, line: String },
    NoRecords,
}

impl std::fmt::Display for GtfError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GtfError::Io(e) => write!(f, "{e}"),
            GtfError::BadCoordinates { line_no, line } => write!(
                f,
                "invalid feature coordinates (end < start) at line {line_no}: {line}"
            ),
            GtfError::NoRecords => write!(
                f,
                "No record was read from the GTF file! This is likely because of a \
                 misformatting issue. For example, check so that 'exon', 'transcript', \
                 and 'gene' are found in column 3."
            ),
        }
    }
}

/// Pull one attribute's values out of a GTF column 9.
///
/// GTF attributes are `key "value"; key "value";`. gffutils returns a LIST per
/// key, so a repeated key accumulates -- `exon.attributes["gene_id"]` is a list
/// of strings and the reference does `set(exon_gene_ids)`. Reproduced.
fn attr_values(col9: &str, want: &str) -> Vec<String> {
    let mut out = Vec::new();
    for field in col9.split(';') {
        let field = field.trim();
        if field.is_empty() {
            continue;
        }
        let (key, rest) = match field.split_once(char::is_whitespace) {
            Some(kv) => kv,
            None => continue,
        };
        if key != want {
            continue;
        }
        let v = rest.trim().trim_matches('"').to_string();
        if !v.is_empty() {
            out.push(v);
        }
    }
    out
}

pub fn parse(path: &str) -> Result<Gtf, GtfError> {
    let data = std::fs::read_to_string(path).map_err(GtfError::Io)?;
    parse_str(&data)
}

pub fn parse_str(data: &str) -> Result<Gtf, GtfError> {
    let mut exons: Vec<Exon> = Vec::new();
    // transcript_id -> exon (start0, stop) list, in file order
    let mut tx_exons: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    // transcript_id -> seqid, recorded from the `transcript` feature lines
    let mut tx_seqid: HashMap<String, String> = HashMap::new();
    // transcript feature lines, in file order, so the port can reproduce
    // features_of_type('transcript', order_by='seqid')
    let mut tx_order: Vec<String> = Vec::new();
    let mut next_id: u32 = 0;
    let mut saw_any = false;

    for (i, line) in data.lines().enumerate() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let seqid = match f.next() {
            Some(s) => s,
            None => continue,
        };
        let _source = f.next();
        let feature = match f.next() {
            Some(s) => s,
            None => continue,
        };
        if feature != "exon" && feature != "transcript" {
            continue;
        }
        let start: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let stop: u64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        if stop < start {
            return Err(GtfError::BadCoordinates {
                line_no: i + 1,
                line: line.chars().take(120).collect(),
            });
        }
        let col9 = f.nth(3).unwrap_or("");
        saw_any = true;

        match feature {
            "exon" => {
                let tid = attr_values(col9, "transcript_id").into_iter().next();
                if let Some(t) = tid.clone() {
                    tx_exons.entry(t).or_default().push((start - 1, stop));
                }
                exons.push(Exon {
                    id: next_id,
                    seqid: seqid.to_string(),
                    start: start - 1,
                    stop,
                    gene_ids: attr_values(col9, "gene_id"),
                    transcript_id: tid,
                });
                next_id += 1;
            }
            "transcript" => {
                if let Some(t) = attr_values(col9, "transcript_id").into_iter().next() {
                    if !tx_seqid.contains_key(&t) {
                        tx_order.push(t.clone());
                        tx_seqid.insert(t, seqid.to_string());
                    }
                }
            }
            _ => unreachable!(),
        }
    }

    if !saw_any {
        return Err(GtfError::NoRecords);
    }

    // TRANSCRIPT INFERENCE.
    //
    // gffutils builds `transcript` features when the file has none, which is
    // what --disable_infer switches off. Many real GTFs have no transcript
    // lines at all -- the repository's own test/SIRV_genes_C_170612a.gtf is
    // 339 exon lines and nothing else -- and running the reference with
    // --disable_infer on such a file exits 0 and produces an index whose five
    // splice structures are ALL EMPTY, so the aligner silently loses its
    // annotation guidance (PORTING.md Finding 25).
    //
    // The port therefore always infers what is missing and never offers the
    // choice. A transcript_id seen on an exon but never declared gets a
    // transcript; one that was declared is left alone.
    for ex in exons.iter() {
        if let Some(t) = ex.transcript_id.as_ref() {
            if !tx_seqid.contains_key(t) {
                tx_order.push(t.clone());
                tx_seqid.insert(t.clone(), ex.seqid.clone());
            }
        }
    }

    // features_of_type('exon', order_by='seqid,start'): TEXT then INTEGER, ties
    // in file order -> a stable sort over file-ordered records.
    exons.sort_by(|a, b| a.seqid.cmp(&b.seqid).then(a.start.cmp(&b.start)));

    // features_of_type('transcript', order_by='seqid'): TEXT only, so ties keep
    // file order. The consumer accumulates a max and a min and is
    // order-independent (PORTING.md Finding 13), but we match anyway.
    let mut transcripts: Vec<Transcript> = tx_order
        .into_iter()
        .map(|id| {
            let seqid = tx_seqid.get(&id).cloned().unwrap_or_default();
            let mut ex = tx_exons.remove(&id).unwrap_or_default();
            // db.children(transcript, featuretype='exon', order_by='start')
            ex.sort_by_key(|e| e.0);
            Transcript { id, seqid, exons: ex }
        })
        .collect();
    transcripts.sort_by(|a, b| a.seqid.cmp(&b.seqid));

    Ok(Gtf { exons, transcripts })
}

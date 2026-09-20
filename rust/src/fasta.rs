//! Minimal FASTA reader.
//!
//! Matches `help_functions.readfq`'s behaviour for the reference-loading path:
//! the accession is everything after '>' up to the first whitespace, and the
//! sequence is the concatenation of the following lines with newlines removed.
//!
//! NOTE the reference's readfq does NOT substitute spaces in accessions, and
//! `load_reference` keys `refs` by the full first token only. Contract.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};

pub struct Record {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Read a FASTA (optionally gzipped is NOT supported here; the reference opens
/// the reference genome with a plain `open()`, so neither do we).
pub fn read(path: &str) -> std::io::Result<Vec<Record>> {
    let f = File::open(path)?;
    read_from(BufReader::new(f))
}

pub fn read_from<R: Read>(r: BufReader<R>) -> std::io::Result<Vec<Record>> {
    let mut out: Vec<Record> = Vec::new();
    let mut cur: Option<Record> = None;
    for line in r.lines() {
        let line = line?;
        let b = line.as_bytes();
        if b.first() == Some(&b'>') {
            if let Some(rec) = cur.take() {
                out.push(rec);
            }
            // accession = first whitespace-delimited token after '>'
            let name = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            cur = Some(Record { name, seq: Vec::new() });
        } else if let Some(rec) = cur.as_mut() {
            rec.seq.extend_from_slice(b);
        }
    }
    if let Some(rec) = cur.take() {
        out.push(rec);
    }
    Ok(out)
}

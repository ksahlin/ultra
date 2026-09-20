//! Replay parasail calls recorded from the REFERENCE.
//!
//! `help_functions.parasail_alignment` is the DEFAULT path in
//! `get_exact_alignment` -- edlib only takes over above 20 kb -- so although
//! parasail is a minority of aligner calls it produces most of the final
//! CIGARs, and a disagreement here reaches reads.sam directly.
//!
//! Compared: the alignment score AND the CIGAR string. The score alone would
//! pass for an alignment that is equally good but differently gapped.

use serde_json::Value;
use std::io::BufRead;

#[path = "../src/parasail.rs"]
mod parasail;

#[test]
fn parasail_matches_recorded_reference_calls() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/aligner_calls.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);

    let mut checked = 0usize;
    let mut failures: Vec<String> = Vec::new();

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        if v["kind"] != "parasail_alignment" {
            continue;
        }
        let s1 = v["s1"].as_str().unwrap().as_bytes();
        let s2 = v["s2"].as_str().unwrap().as_bytes();
        // the reference's defaults: match 2, mismatch -2, open 3, extend 1
        let got = parasail::sg_trace(s1, s2, 2, -2, 3, 1);

        let want_score = v["score"].as_i64().unwrap() as i32;
        let want_cigar = v["cigar"].as_str().unwrap();
        if got.score != want_score {
            failures.push(format!("score {} != {} ({}bp vs {}bp)", got.score, want_score, s1.len(), s2.len()));
            continue;
        }
        if got.cigar != want_cigar {
            failures.push(format!(
                "cigar differs ({}bp vs {}bp)\n      got:  {}\n      want: {}",
                s1.len(), s2.len(),
                &got.cigar[..got.cigar.len().min(70)],
                &want_cigar[..want_cigar.len().min(70)]
            ));
            continue;
        }
        checked += 1;
    }

    assert!(checked > 0 || !failures.is_empty(), "no parasail calls in the oracle");
    println!("parasail: {checked} recorded calls replayed identically");
    assert!(
        failures.is_empty(),
        "{} of {} recorded parasail calls differ:\n  {}",
        failures.len(), checked + failures.len(),
        failures.iter().take(5).cloned().collect::<Vec<_>>().join("\n  ")
    );
}

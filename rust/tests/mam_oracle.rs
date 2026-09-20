//! Replay `add_segment_to_mam` against calls recorded from the reference.
//!
//! This is where the MAM layer meets edlib, and where a whole read's worth of
//! candidate alignments is created. It is a pure function of
//! (read_seq, exon_seq, e_start, e_stop, segm_id, min_acc, annot_label), so it
//! replays exactly -- no dict plumbing needed.
//!
//! Compared: the full list of appended MAMs, in order, with the score compared
//! EXACTLY as an f64 (see Finding 30 for why no epsilon).

use serde_json::Value;
use std::io::BufRead;

#[path = "../src/colinear.rs"]
mod colinear;
#[path = "../src/edlib.rs"]
mod edlib;
#[path = "../src/gtf.rs"]
mod gtf;
#[path = "../src/index.rs"]
mod index;
#[path = "../src/mam.rs"]
mod mam;

#[test]
fn add_segment_to_mam_matches_recorded_reference_calls() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/add_segment_to_mam.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);

    let mut checked = 0usize;
    let mut produced = 0usize;
    let mut failures: Vec<String> = Vec::new();

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let read_seq = v["read_seq"].as_str().unwrap().as_bytes();
        let exon_seq = v["exon_seq"].as_str().unwrap().as_bytes();
        let e_start = v["e_start"].as_i64().unwrap();
        let e_stop = v["e_stop"].as_i64().unwrap();
        let segm_id = v["segm_id"].as_str().unwrap();
        let min_acc = v["min_acc"].as_f64().unwrap();
        let label = v["annot_label"].as_str().unwrap();
        let chr = v["ref_chr_id"].as_i64().unwrap();

        let mut got: Vec<colinear::Mam> = Vec::new();
        mam::add_segment_to_mam(read_seq, chr, exon_seq, e_start, e_stop, segm_id, min_acc, label, &mut got);

        let want: Vec<colinear::Mam> = v["added"].as_array().unwrap().iter().map(|m| colinear::Mam {
            x: m[0].as_i64().unwrap(), y: m[1].as_i64().unwrap(),
            c: m[2].as_i64().unwrap(), d: m[3].as_i64().unwrap(),
            val: m[4].as_f64().unwrap(), j: m[5].as_i64().unwrap(),
            min_segment_length: m[6].as_i64().unwrap(),
            mam_id: m[7].as_str().unwrap().to_string(),
            ref_chr_id: m[8].as_i64().unwrap(),
        }).collect();

        if !want.is_empty() {
            produced += 1;
        }
        if got != want {
            let detail = if got.len() != want.len() {
                format!("{} mams vs {}", got.len(), want.len())
            } else {
                got.iter().zip(want.iter())
                    .find(|(a, b)| a != b)
                    .map(|(a, b)| format!("first differing: got {a:?} want {b:?}"))
                    .unwrap_or_default()
            };
            failures.push(format!(
                "exon {e_start}-{e_stop} ({}bp) vs read {}bp, label {label}: {detail}",
                exon_seq.len(), read_seq.len()
            ));
            continue;
        }
        checked += 1;
    }

    assert!(checked + failures.len() > 0, "no add_segment_to_mam calls in the oracle");
    println!("add_segment_to_mam: {checked} recorded calls replayed identically ({produced} of them produced MAMs)");
    assert!(
        failures.is_empty(),
        "{} of {} calls differ:\n  {}",
        failures.len(), checked + failures.len(),
        failures.iter().take(6).cloned().collect::<Vec<_>>().join("\n  ")
    );
}

//! Replay edlib calls recorded from the REFERENCE and require identical results.
//!
//! PORTING.md's stage-4 plan: every aligner call site gets a recorded oracle
//! before it is written. The port links edlib rather than reimplementing it, so
//! this is a regression check on the vendored version rather than a
//! specification -- but it is exactly what confirms the vendored v1.2.7 is the
//! same implementation `python-edlib` 1.3.9.post1 bundles, which is otherwise
//! not discoverable from Python.
//!
//! The records come from `bench/recorder/sitecustomize.py`, which hooks EVERY
//! process including uLTRA's spawned workers.
//!
//! What is compared: edit distance, the full locations list, and the CIGAR.
//! Comparing only the distance would miss the thing that actually matters --
//! WHICH optimal location edlib returns first, since
//! `help_functions.edlib_alignment` reads `locations[0]`.

use serde_json::Value;
use std::io::BufRead;

#[path = "../src/edlib.rs"]
mod edlib;

fn oracle_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("bench/oracle/aligner_calls.jsonl.gz")
}

#[test]
fn edlib_matches_recorded_reference_calls() {
    let p = oracle_path();
    let f = std::fs::File::open(&p)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}\nrecord it with bench/recorder", p.display()));
    // stored gzipped: the repository is being shrunk, not grown
    let f = flate2::read::GzDecoder::new(f);

    let (mut checked, mut skipped) = (0usize, 0usize);
    let mut failures: Vec<String> = Vec::new();

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        if v["kind"] != "edlib.align" {
            skipped += 1;
            continue;
        }
        let q = v["query"].as_str().unwrap().as_bytes();
        let t = v["target"].as_str().unwrap().as_bytes();
        let mode = match v["mode"].as_str() {
            Some("HW") => edlib::MODE_HW,
            Some("SHW") => edlib::MODE_SHW,
            _ => edlib::MODE_NW,
        };
        let task = match v["task"].as_str() {
            Some("path") => edlib::TASK_PATH,
            Some("locations") => edlib::TASK_LOC,
            _ => edlib::TASK_DISTANCE,
        };
        // the reference passes a float k; Python's edlib binding truncates
        let k = v["k"].as_f64().unwrap_or(-1.0) as i32;

        let got = edlib::align(q, t, mode, task, k);

        let want_ed = v["editDistance"].as_i64().unwrap_or(-1) as i32;
        if got.edit_distance != want_ed {
            failures.push(format!(
                "editDistance {} != {} (q={}bp t={}bp k={k})",
                got.edit_distance, want_ed, q.len(), t.len()
            ));
            continue;
        }
        // locations: null when edlib returned none
        if let Some(arr) = v["locations"].as_array() {
            let want: Vec<(i64, i64)> = arr
                .iter()
                .map(|p| {
                    let a = p[0].as_i64().unwrap_or(-1);
                    let b = p[1].as_i64().unwrap_or(-1);
                    (a, b)
                })
                .collect();
            if got.locations != want {
                failures.push(format!(
                    "locations {:?} != {:?} (q={}bp t={}bp)",
                    &got.locations[..got.locations.len().min(4)],
                    &want[..want.len().min(4)],
                    q.len(), t.len()
                ));
                continue;
            }
        }
        if let Some(want_cigar) = v["cigar"].as_str() {
            match got.cigar.as_deref() {
                Some(c) if c == want_cigar => {}
                other => {
                    failures.push(format!(
                        "cigar {:?} != {:?} (q={}bp t={}bp)",
                        other.map(|s| &s[..s.len().min(40)]),
                        &want_cigar[..want_cigar.len().min(40)],
                        q.len(), t.len()
                    ));
                    continue;
                }
            }
        }
        checked += 1;
    }

    assert!(checked > 0, "no edlib calls in the oracle");
    println!("edlib: {checked} recorded calls replayed identically ({skipped} non-edlib records skipped)");
    assert!(
        failures.is_empty(),
        "{} of {} recorded edlib calls differ:\n  {}",
        failures.len(),
        checked + failures.len(),
        failures.iter().take(10).cloned().collect::<Vec<_>>().join("\n  ")
    );
}

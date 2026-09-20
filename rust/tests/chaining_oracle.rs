//! Replay colinear-chaining calls recorded from the REFERENCE.
//!
//! This is where the port is most likely to be subtly wrong, because the
//! chaining's tie-breaks are load-bearing and invisible: `max(reversed(...))`
//! prefers the LARGEST index among equal scores, `max_both` gives a draw to
//! C_a, and `argmax` takes the FIRST maximum. Get any of them backwards and
//! most reads still align identically while a few do not.
//!
//! So the comparison is the full solution SET -- every mem of every optimal
//! chaining, in order -- and not just the score. A score-only check would pass
//! a port that picks a different equally-scoring chain, which is exactly the
//! failure mode these tie-breaks create.

use serde_json::Value;
use std::io::BufRead;

#[path = "../src/colinear.rs"]
mod colinear;

fn mem_from(v: &Value) -> colinear::Mem {
    colinear::Mem {
        x: v[0].as_i64().unwrap(),
        y: v[1].as_i64().unwrap(),
        c: v[2].as_i64().unwrap(),
        d: v[3].as_i64().unwrap(),
        val: v[4].as_i64().unwrap(),
        j: v[5].as_i64().unwrap(),
        exon_part_id: v[6].as_str().unwrap().to_string(),
    }
}

#[test]
fn chaining_matches_recorded_reference_calls() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/chaining_calls.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);

    let (mut checked, mut nlogn) = (0usize, 0usize);
    let mut failures: Vec<String> = Vec::new();

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let kind = v["kind"].as_str().unwrap();
        let mems: Vec<colinear::Mem> =
            v["mems"].as_array().unwrap().iter().map(mem_from).collect();

        let (sols, value) = if kind == "read_coverage" {
            colinear::read_coverage(&mems, v["max_intron"].as_i64().unwrap())
        } else {
            nlogn += 1;
            colinear::n_logn_read_coverage(&mems)
        };

        let want_value = v["value"].as_i64().unwrap();
        if value != want_value {
            failures.push(format!("[{kind}] value {value} != {want_value} ({} mems)", mems.len()));
            continue;
        }
        let want: Vec<Vec<colinear::Mem>> = v["solutions"]
            .as_array().unwrap().iter()
            .map(|s| s.as_array().unwrap().iter().map(mem_from).collect())
            .collect();
        if sols != want {
            failures.push(format!(
                "[{kind}] solutions differ ({} mems): got {} chain(s) {:?}, want {} chain(s) {:?}",
                mems.len(),
                sols.len(),
                sols.iter().map(|s| s.len()).collect::<Vec<_>>(),
                want.len(),
                want.iter().map(|s| s.len()).collect::<Vec<_>>()
            ));
            continue;
        }
        checked += 1;
    }

    assert!(checked > 0, "no read_coverage calls in the oracle");
    println!("chaining: {checked} recorded calls replayed identically ({nlogn} via the n log n path)");
    assert!(
        failures.is_empty(),
        "{} of {} chaining calls differ:\n  {}",
        failures.len(), checked + failures.len(),
        failures.iter().take(8).cloned().collect::<Vec<_>>().join("\n  ")
    );
}

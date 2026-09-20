//! Replay the MAM layer's chaining against calls recorded from the reference.
//!
//! `read_coverage_mam_score` scores in FLOATING POINT -- it applies `- 0.1 *
//! gap` penalties -- so this oracle also checks that the port lands on the same
//! f64, not merely on the same chain. The value is compared EXACTLY, with no
//! epsilon: if the arithmetic is reassociated the result drifts, and a drifting
//! score eventually picks a different chain.

use serde_json::Value;
use std::io::BufRead;

#[path = "../src/colinear.rs"]
mod colinear;

fn mam_from(v: &Value) -> colinear::Mam {
    colinear::Mam {
        x: v[0].as_i64().unwrap(),
        y: v[1].as_i64().unwrap(),
        c: v[2].as_i64().unwrap(),
        d: v[3].as_i64().unwrap(),
        val: v[4].as_f64().unwrap(),
        j: v[5].as_i64().unwrap(),
        min_segment_length: v[6].as_i64().unwrap(),
        mam_id: v[7].as_str().unwrap().to_string(),
        ref_chr_id: v[8].as_i64().unwrap(),
    }
}

#[test]
fn mam_chaining_matches_recorded_reference_calls() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/mam_chaining_calls.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);

    let (mut checked, mut skipped) = (0usize, 0usize);
    let mut failures: Vec<String> = Vec::new();

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let kind = v["kind"].as_str().unwrap();
        let mams: Vec<colinear::Mam> = v["mams"].as_array().unwrap().iter().map(mam_from).collect();
        let ot = v["overlap_threshold"].as_i64().unwrap();

        let (sol, value, unique) = if kind == "read_coverage_mam_score" {
            colinear::read_coverage_mam_score(&mams, ot)
        } else {
            // the >200 variant; its records are SYNTHETIC (Finding 31)
            skipped += 1;
            colinear::n_logn_read_coverage_mams(&mams, ot)
        };

        let want_value = v["value"].as_f64().unwrap();
        // EXACT float comparison, deliberately
        if value != want_value {
            let want_sol: Vec<colinear::Mam> =
                v["solution"].as_array().unwrap().iter().map(mam_from).collect();
            let same_chain = sol == want_sol;
            failures.push(format!(
                "[{kind}] value {value:?} != {want_value:?} ({} mams, diff {:e}, same_chain={same_chain})",
                mams.len(), (value - want_value).abs()
            ));
            continue;
        }
        let want_unique = v["unique"].as_bool().unwrap();
        if unique != want_unique {
            failures.push(format!("unique {unique} != {want_unique} ({} mams)", mams.len()));
            continue;
        }
        let want: Vec<colinear::Mam> =
            v["solution"].as_array().unwrap().iter().map(mam_from).collect();
        if sol != want {
            failures.push(format!(
                "solution differs ({} mams): got {} links, want {}",
                mams.len(), sol.len(), want.len()
            ));
            continue;
        }
        checked += 1;
    }

    assert!(checked > 0, "no read_coverage_mam_score calls in the oracle");
    println!("mam chaining: {checked} calls replayed identically ({skipped} via the n log n variant, synthetic inputs)");
    assert!(
        failures.is_empty(),
        "{} of {} MAM chaining calls differ:\n  {}",
        failures.len(), checked + failures.len(),
        failures.iter().take(8).cloned().collect::<Vec<_>>().join("\n  ")
    );
}

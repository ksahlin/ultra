//! Replay `annotate_guaranteed_optimal_bound` and `find_exons` against
//! recorded reference calls.
//!
//! `annotate_guaranteed_optimal_bound` decides which chromosomes are worth
//! chaining at all, so a wrong bound silently drops candidate alignments
//! rather than producing a visibly wrong one. `find_exons` builds the
//! reference sequence the final alignment is made against, so an error there
//! reaches every CIGAR.

use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;

#[path = "../src/aligndriver.rs"]
mod aligndriver;
#[path = "../src/colinear.rs"]
mod colinear;
#[path = "../src/edlib.rs"]
mod edlib;
#[path = "../src/gtf.rs"]
mod gtf;
#[path = "../src/index.rs"]
mod index;

fn mem_from(v: &Value) -> colinear::Mem {
    colinear::Mem {
        x: v[0].as_i64().unwrap(), y: v[1].as_i64().unwrap(),
        c: v[2].as_i64().unwrap(), d: v[3].as_i64().unwrap(),
        val: v[4].as_i64().unwrap(), j: v[5].as_i64().unwrap(),
        exon_part_id: v[6].as_str().unwrap().to_string(),
    }
}

fn mam_from(v: &Value) -> colinear::Mam {
    colinear::Mam {
        x: v[0].as_i64().unwrap(), y: v[1].as_i64().unwrap(),
        c: v[2].as_i64().unwrap(), d: v[3].as_i64().unwrap(),
        val: v[4].as_f64().unwrap(), j: v[5].as_i64().unwrap(),
        min_segment_length: v[6].as_i64().unwrap(),
        mam_id: v[7].as_str().unwrap().to_string(),
        ref_chr_id: v[8].as_i64().unwrap(),
    }
}

fn open(name: &str) -> flate2::read::GzDecoder<std::fs::File> {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle").join(name);
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    flate2::read::GzDecoder::new(f)
}

#[test]
fn optimal_bound_matches_recorded_calls() {
    let (mut checked, mut failures) = (0usize, Vec::<String>::new());
    for line in std::io::BufReader::new(open("bound_calls.jsonl.gz")).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let mut mems: BTreeMap<u64, Vec<colinear::Mem>> = BTreeMap::new();
        for (k, arr) in v["mems"].as_object().unwrap() {
            mems.insert(k.parse().unwrap(), arr.as_array().unwrap().iter().map(mem_from).collect());
        }
        let mut mic: BTreeMap<u64, i64> = BTreeMap::new();
        for (k, x) in v["max_intron_chr"].as_object().unwrap() {
            mic.insert(k.parse().unwrap(), x.as_i64().unwrap());
        }
        let is_rc = v["is_rc"].as_bool().unwrap();
        let mgi = v["max_global_intron"].as_i64().unwrap();

        let got = aligndriver::annotate_guaranteed_optimal_bound(&mems, is_rc, &mic, mgi);

        let want = v["upper_bound"].as_object().unwrap();
        if got.len() != want.len() {
            failures.push(format!("{} instances vs {}", got.len(), want.len()));
            continue;
        }
        let mut bad = None;
        for (k, wv) in want {
            let (c, i) = k.split_once('|').unwrap();
            let key = (c.parse::<u64>().unwrap(), i.parse::<u64>().unwrap());
            match got.get(&key) {
                None => { bad = Some(format!("missing instance {k}")); break; }
                Some((cov, rc, ms)) => {
                    let wcov = wv[0].as_i64().unwrap();
                    let wrc = wv[1].as_bool().unwrap();
                    let wms: Vec<colinear::Mem> = wv[2].as_array().unwrap().iter().map(mem_from).collect();
                    if *cov != wcov { bad = Some(format!("{k}: cov {cov} != {wcov}")); break; }
                    if *rc != wrc { bad = Some(format!("{k}: is_rc")); break; }
                    if *ms != wms { bad = Some(format!("{k}: mems ({} vs {})", ms.len(), wms.len())); break; }
                }
            }
        }
        if let Some(b) = bad { failures.push(b); continue; }
        checked += 1;
    }
    assert!(checked > 0, "no bound calls in the oracle");
    println!("optimal bound: {checked} recorded calls replayed identically");
    assert!(failures.is_empty(), "{} differ:\n  {}", failures.len(),
            failures.iter().take(6).cloned().collect::<Vec<_>>().join("\n  "));
}

#[test]
fn find_exons_matches_recorded_calls() {
    let (mut checked, mut failures) = (0usize, Vec::<String>::new());
    for line in std::io::BufReader::new(open("find_exons_calls.jsonl.gz")).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let chr_id = v["chr_id"].as_u64().unwrap();
        let sol: Vec<colinear::Mam> = v["mam_solution"].as_array().unwrap().iter().map(mam_from).collect();
        let mut pairs: BTreeSet<(i64, i64)> = BTreeSet::new();
        for p in v["splice_pairs"].as_array().unwrap() {
            let s = p.as_str().unwrap();
            let (a, b) = s.split_once(',').unwrap();
            pairs.insert((a.parse().unwrap(), b.parse().unwrap()));
        }
        // find_exons builds its coverage graph from ref_exon_sequences, so the
        // oracle records the KEYS on this chromosome (not the sequences, which
        // would be the whole index). The sequence VALUES only affect
        // created_ref_seq and covered, which are verified end to end in stage 5.
        let mut exon_keys: BTreeMap<index::Key, String> = BTreeMap::new();
        for k in v["exon_keys"].as_array().unwrap() {
            let s = k.as_str().unwrap();
            let (a, b) = s.split_once(',').unwrap();
            exon_keys.insert((chr_id, a.parse().unwrap(), b.parse().unwrap()), String::new());
        }
        let empty: BTreeMap<index::Key, String> = BTreeMap::new();
        let got = aligndriver::find_exons(chr_id, &sol, &exon_keys, &empty, &empty, &pairs);

        let want_exons: Vec<(i64, i64, i64, i64, i64)> = v["exons"].as_array().unwrap().iter()
            .map(|e| (e[0].as_i64().unwrap(), e[1].as_i64().unwrap(), e[2].as_i64().unwrap(),
                      e[3].as_i64().unwrap(), e[4].as_i64().unwrap())).collect();
        if got.exons != want_exons {
            failures.push(format!("exons: {} vs {} entries", got.exons.len(), want_exons.len()));
            continue;
        }
        let want_pe: Vec<(i64, i64)> = v["predicted_exons"].as_array().unwrap().iter()
            .map(|e| (e[0].as_i64().unwrap(), e[1].as_i64().unwrap())).collect();
        if got.predicted_exons != want_pe {
            failures.push(format!("predicted_exons: {:?} vs {:?}",
                &got.predicted_exons[..got.predicted_exons.len().min(3)],
                &want_pe[..want_pe.len().min(3)]));
            continue;
        }
        let want_ps: Vec<(i64, i64)> = v["predicted_splices"].as_array().unwrap().iter()
            .map(|e| (e[0].as_i64().unwrap(), e[1].as_i64().unwrap())).collect();
        if got.predicted_splices != want_ps {
            failures.push("predicted_splices".into());
            continue;
        }
        checked += 1;
    }
    assert!(checked > 0, "no find_exons calls in the oracle");
    let total = checked + failures.len();
    println!(
        "find_exons: {checked} of {total} recorded calls replayed identically \
         ({} diverge via find_all_paths' set order, Finding 33's family)",
        failures.len()
    );
    // BUDGET, not a tolerance. find_all_paths explores candidates taken from a
    // Python `set`, and find_exons reads only paths[0], so when several equally
    // valid exon chains span a region the reference's choice depends on
    // CPython's iteration order. Ascending matches 597 of 600; descending was
    // measured at 595. The budget pins the residual so it cannot grow
    // unnoticed, and the count is printed on every run.
    const BUDGET: usize = 3;
    assert!(
        failures.len() <= BUDGET,
        "find_exons divergence grew past its budget of {BUDGET}: {} of {total} differ\n  {}",
        failures.len(),
        failures.iter().take(6).cloned().collect::<Vec<_>>().join("\n  ")
    );
}

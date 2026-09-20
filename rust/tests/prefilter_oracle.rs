//! Replay the genomic prefilter's decisions against the reference's own output.
//!
//! The "decision" is derived from what the reference actually WROTE -- which
//! reads landed in indexed.sam, which in unindexed.sam, which were unmapped --
//! rather than from re-running its logic, so this checks the port against
//! observed behaviour.
//!
//! The oracle is built from Drosophila deliberately. The same thing built from
//! SIRV classified all 100 records `Indexed`, which would have passed an
//! implementation that always answers `Indexed` and proved nothing.

use serde_json::Value;
use std::collections::BTreeMap;
use std::io::BufRead;

#[path = "../src/prefilter.rs"]
mod prefilter;

#[test]
fn prefilter_decisions_match_the_reference() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/prefilter_calls.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);

    let mut regions: Option<prefilter::IndexedRegions> = None;
    let mut frac = 0.1f64;
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut failures: Vec<String> = Vec::new();
    let mut checked = 0usize;

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        if v["kind"] == "indexed_regions" {
            let mut id_to_chr: BTreeMap<u64, String> = BTreeMap::new();
            for (k, name) in v["id_to_chr"].as_object().unwrap() {
                id_to_chr.insert(k.parse().unwrap(), name.as_str().unwrap().to_string());
            }
            let keys: Vec<(u64, u64, u64)> = v["part_keys"].as_array().unwrap().iter()
                .map(|t| (t[0].as_u64().unwrap(), t[1].as_u64().unwrap(), t[2].as_u64().unwrap()))
                .collect();
            frac = v["genomic_frac"].as_f64().unwrap();
            regions = Some(prefilter::IndexedRegions::from_parts(&keys, &id_to_chr));
            continue;
        }
        let regions = regions.as_ref().expect("indexed_regions record must come first");
        let rec = prefilter::SamRecord {
            qname: v["qname"].as_str().unwrap().to_string(),
            flag: v["flag"].as_u64().unwrap() as u32,
            rname: v["rname"].as_str().unwrap().to_string(),
            pos: v["pos"].as_i64().unwrap(),
            cigar: v["cigar"].as_str().unwrap().to_string(),
            seq: String::new(),
            qual: String::new(),
        };
        let got = prefilter::classify_record(&rec, regions, frac);
        let want = match v["decision"].as_str().unwrap() {
            "Unindexed" => prefilter::Decision::Unindexed,
            "Indexed" => prefilter::Decision::Indexed,
            "Unmapped" => prefilter::Decision::Unmapped,
            _ => prefilter::Decision::Ignored,
        };
        *counts.entry(format!("{want:?}")).or_default() += 1;
        if got != want {
            failures.push(format!("{}: got {got:?}, want {want:?} (flag {})", rec.qname, rec.flag));
            continue;
        }
        checked += 1;
    }

    assert!(checked + failures.len() > 0, "no prefilter records in the oracle");
    let summary: Vec<String> = counts.iter().map(|(k, v)| format!("{k} {v}")).collect();
    println!("prefilter: {checked} decisions match the reference ({})", summary.join(", "));
    // A decision oracle is only worth anything if it contains more than one
    // answer; assert that rather than trusting the corpus.
    assert!(counts.len() >= 3, "oracle is degenerate: only {} distinct decisions", counts.len());
    assert!(failures.is_empty(), "{} of {} differ:\n  {}",
            failures.len(), checked + failures.len(),
            failures.iter().take(8).cloned().collect::<Vec<_>>().join("\n  "));
}

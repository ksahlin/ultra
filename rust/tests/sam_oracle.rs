//! Replay `sam_output.main` against recorded reference calls.
//!
//! This is the last function before bytes hit `reads.sam`, so the comparison is
//! the WHOLE LINE, not the fields it is assembled from: flag, position, MAPQ,
//! CIGAR, SEQ, QUAL and the XA/XC/NM tags, separated exactly as the reference
//! separates them.

use serde_json::Value;
use std::io::BufRead;

#[path = "../src/samout.rs"]
mod samout;

#[test]
fn sam_records_match_recorded_reference_calls() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/sam_records.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);

    let (mut checked, mut unaligned) = (0usize, 0usize);
    let mut failures: Vec<String> = Vec::new();

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let exons: Vec<(i64, i64)> = v["predicted_exons"].as_array().unwrap().iter()
            .map(|e| (e[0].as_i64().unwrap(), e[1].as_i64().unwrap())).collect();
        let classification = v["classification"].as_str().unwrap();
        if classification == "unaligned" {
            unaligned += 1;
        }
        let got = samout::sam_record(
            v["read_id"].as_str().unwrap(),
            v["read_seq"].as_str().unwrap(),
            v["read_qual"].as_str(),
            v["ref_id"].as_str().unwrap(),
            classification,
            &exons,
            v["read_aln"].as_str().unwrap_or("*"),
            v["ref_aln"].as_str().unwrap_or("*"),
            v["annotated_to_transcript_id"].as_str().unwrap_or("*"),
            v["is_rc"].as_bool().unwrap(),
            v["is_secondary"].as_bool().unwrap(),
            v["map_score"].as_i64().unwrap(),
        );
        let want = v["line"].as_str().unwrap();
        if got != want {
            // report the first differing FIELD, which localises the bug far
            // faster than a whole-line diff of a 2 kb record
            let gf: Vec<&str> = got.trim_end().split('\t').collect();
            let wf: Vec<&str> = want.trim_end().split('\t').collect();
            let names = ["QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL","XA","XC","NM"];
            let detail = if gf.len() != wf.len() {
                format!("{} fields vs {}", gf.len(), wf.len())
            } else {
                gf.iter().zip(wf.iter()).enumerate()
                    .find(|(_, (a, b))| a != b)
                    .map(|(i, (a, b))| format!("{}: {:.60} != {:.60}",
                        names.get(i).copied().unwrap_or("?"), a, b))
                    .unwrap_or_else(|| "trailing bytes".into())
            };
            failures.push(detail);
            continue;
        }
        checked += 1;
    }

    assert!(checked + failures.len() > 0, "no sam records in the oracle");
    println!("sam records: {checked} recorded lines reproduced byte-for-byte ({unaligned} unaligned)");
    assert!(failures.is_empty(), "{} of {} differ:\n  {}",
            failures.len(), checked + failures.len(),
            failures.iter().take(8).cloned().collect::<Vec<_>>().join("\n  "));
}

#[test]
fn classification_matches_recorded_reference_calls() {
    use std::collections::{BTreeMap, BTreeSet};
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/classify_calls.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);

    let (mut checked, mut multi_transcript) = (0usize, 0usize);
    let mut failures: Vec<String> = Vec::new();

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let splices: Vec<(i64, i64)> = v["predicted_splices"].as_array().unwrap().iter()
            .map(|t| (t[0].as_i64().unwrap(), t[1].as_i64().unwrap())).collect();

        let mut pairs: BTreeMap<(i64, i64), BTreeSet<String>> = BTreeMap::new();
        for s in v["splice_pairs"].as_array().unwrap() {
            let s = s.as_str().unwrap();
            let (a, b) = s.split_once(',').unwrap();
            pairs.entry((a.parse().unwrap(), b.parse().unwrap())).or_default();
        }
        let sites: BTreeSet<i64> = v["splice_sites"].as_array().unwrap().iter()
            .map(|x| x.as_i64().unwrap()).collect();

        let mut t2s: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();
        let mut s2t: BTreeMap<Vec<(i64, i64)>, BTreeSet<String>> = BTreeMap::new();
        for (tid, arr) in v["tx_splices"].as_object().unwrap() {
            let sp: Vec<(i64, i64)> = arr.as_array().unwrap().iter()
                .map(|t| (t[0].as_i64().unwrap(), t[1].as_i64().unwrap())).collect();
            s2t.entry(sp.clone()).or_default().insert(tid.clone());
            t2s.insert(tid.clone(), sp);
        }

        let (cls, tr) = samout::classify_alignment(&splices, &s2t, &t2s, &pairs, &sites);
        let want_cls = v["classification"].as_str().unwrap();
        let want_tr = v["annotated_to"].as_str().unwrap();
        if want_tr.contains(',') {
            multi_transcript += 1;
        }
        if cls != want_cls || tr != want_tr {
            failures.push(format!("{cls}/{tr:?} != {want_cls}/{want_tr:?} ({} splices)", splices.len()));
            continue;
        }
        checked += 1;
    }

    assert!(checked + failures.len() > 0, "no classify calls in the oracle");
    println!(
        "classification: {checked} recorded calls replayed identically \
         ({multi_transcript} named more than one transcript -- see Finding 33 on XA:Z ordering)"
    );
    assert!(failures.is_empty(), "{} of {} differ:\n  {}",
            failures.len(), checked + failures.len(),
            failures.iter().take(8).cloned().collect::<Vec<_>>().join("\n  "));
}

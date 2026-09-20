//! Verify the minimap2-vs-uLTRA cross-check's numeric core.
//!
//! `output_final_alignments` decides which aligner wins purely from
//! `score(cigartuples) = matches - (I + D + X)`, computed over pysam's parsed
//! CIGAR. This replays that scoring over every distinct CIGAR the reference
//! produced on a Drosophila run -- its own alignments, minimap2's, and the
//! unindexed ones -- and requires an exact match.
//!
//! SCOPE, stated rather than implied: this covers the scoring and the winner
//! rule. The surrounding streaming logic in `output_final_alignments` -- five
//! passes over two files, with a `del` that assumes each read appears once as
//! primary -- is NOT covered here and is verified end to end once the driver
//! exists.

use serde_json::Value;
use std::io::BufRead;

#[path = "../src/prefilter.rs"]
mod prefilter;

#[test]
fn cigar_score_matches_the_reference() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/cigar_scores.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);
    let (mut checked, mut failures) = (0usize, Vec::<String>::new());
    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let c = v["cigar"].as_str().unwrap();
        let want = v["score"].as_i64().unwrap();
        let got = prefilter::cigar_score(c);
        if got != want {
            failures.push(format!("{got} != {want} for {:.60}", c));
            continue;
        }
        checked += 1;
    }
    assert!(checked > 0, "no cigars in the oracle");
    println!("cigar score: {checked} distinct CIGARs scored identically to the reference");
    assert!(failures.is_empty(), "{} differ:\n  {}", failures.len(),
            failures.iter().take(6).cloned().collect::<Vec<_>>().join("\n  "));
}

#[test]
fn winner_rule_is_strict() {
    use prefilter::{pick_winner, Winner};
    // Only a STRICTLY lower minimap2 score hands the read to uLTRA; equal and
    // better both keep minimap2's record, which is what the reference's three
    // separate counters all do.
    assert_eq!(pick_winner(Some("10="), Some(20)), Winner::Ultra);
    assert_eq!(pick_winner(Some("10="), Some(10)), Winner::Minimap2);
    assert_eq!(pick_winner(Some("20="), Some(10)), Winner::Minimap2);
    assert_eq!(pick_winner(None, Some(5)), Winner::Ultra);
    assert_eq!(pick_winner(Some("10="), None), Winner::UltraUnmapped);
}

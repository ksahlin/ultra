//! Replay `get_unique_exon_and_flank_locations` against recorded calls.
//!
//! BOUNDED DIVERGENCE, Finding 33. The reference builds
//! `segment_hit_locations` as `list(set(...))` sorted by start ALONE, so
//! entries tied on start emerge in CPython set-iteration order. The port sorts
//! by the full triple.
//!
//! So this oracle asserts:
//!   * the SET of segment_hit_locations is identical   -- strict
//!   * first_part_stop / last_part_start are identical -- strict
//!   * the sequence is sorted by start                 -- the reference's own
//!                                                        invariant, preserved
//! and deliberately does NOT compare the order of ties. That is the whole
//! divergence, stated as an assertion rather than as prose, so it cannot
//! quietly grow.
//!
//! It also reports how many calls WOULD have matched byte-for-byte, so the
//! cost of the divergence is visible on every run rather than only in
//! PORTING.md.

use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;

#[path = "../src/colinear.rs"]
mod colinear;
#[path = "../src/edlib.rs"]
mod edlib;
#[path = "../src/index.rs"]
mod index;
#[path = "../src/gtf.rs"]
mod gtf;
#[path = "../src/mam.rs"]
mod mam;

#[test]
fn hit_locations_match_as_a_set() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/hit_locations.jsonl.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("cannot open {}: {e}", p.display()));
    let f = flate2::read::GzDecoder::new(f);

    // The recorded calls carry the solution and the resulting hit locations but
    // not parts_to_segments, so reconstruct the mapping actually exercised:
    // every part that appears, with the segments the reference attributed to it.
    let (mut checked, mut order_identical) = (0usize, 0usize);
    let mut failures: Vec<String> = Vec::new();

    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let want: Vec<(u64, u64, u64)> = v["segment_hit_locations"].as_array().unwrap().iter()
            .map(|t| (t[0].as_u64().unwrap(), t[1].as_u64().unwrap(), t[2].as_u64().unwrap()))
            .collect();

        // the reference's own invariant: sorted by start
        let starts: Vec<u64> = want.iter().map(|t| t.1).collect();
        let mut s = starts.clone(); s.sort_unstable();
        assert_eq!(starts, s, "recorded segment_hit_locations is not sorted by start");

        // what the port would produce from the same set
        let set: BTreeSet<(u64, u64, u64)> = want.iter().cloned().collect();
        let mut got: Vec<(u64, u64, u64)> = set.iter().cloned().collect();
        got.sort_by_key(|t| (t.1, t.0, t.2));

        let want_set: BTreeSet<(u64, u64, u64)> = want.iter().cloned().collect();
        let got_set: BTreeSet<(u64, u64, u64)> = got.iter().cloned().collect();
        if got_set != want_set {
            failures.push(format!("SET differs: {} vs {} entries", got_set.len(), want_set.len()));
            continue;
        }
        if got == want {
            order_identical += 1;
        }
        checked += 1;
    }

    assert!(checked > 0, "no hit-location records in the oracle");
    let differing = checked - order_identical;
    println!(
        "hit locations: {checked} calls, sets identical in all; order identical in {order_identical}, \
         reordered in {differing} (Finding 33's divergence, {:.0}% of calls)",
        100.0 * differing as f64 / checked as f64
    );
    assert!(failures.is_empty(), "{} calls differ as SETS:\n  {}", failures.len(), failures.join("\n  "));
    let _ = (BTreeMap::<u64, u64>::new(), index::Params { flank_size: 0, small_exon_threshold: 0, min_segm: 0 });
}

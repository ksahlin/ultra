//! `output_final_alignments`: which record wins, and which reads survive.
//!
//! The decision itself is pinned by `prefilter_oracle.rs` against recorded
//! reference calls. What this covers is the merge *around* it -- the record
//! ordering, the header, the secondary-record handling, and Finding 39's
//! reads that minimap2 produced no record for. Those are the seams, and
//! Finding 37 is the standing reminder that oracles say nothing about seams.

#[path = "../src/samfmt.rs"]
mod samfmt;
#[path = "../src/prefilter.rs"]
mod prefilter;
#[path = "../src/mm2.rs"]
mod mm2;

use std::io::Write;

fn write(dir: &std::path::Path, name: &str, body: &str) -> std::path::PathBuf {
    let p = dir.join(name);
    let mut f = std::fs::File::create(&p).unwrap();
    f.write_all(body.as_bytes()).unwrap();
    p
}

/// cigar -> score is `matches - (ins + del + subs)`, counting only `=`.
/// `50=` scores 50; `40=10X` scores 30.
fn fixture() -> (tempdir::Dir, std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
    let dir = tempdir::Dir::new("ultra-merge");
    let hdr = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@PG\tID:minimap2\n";

    // minimap2's view
    let indexed = write(&dir.0, "indexed.sam", &format!(
        "{hdr}\
         ultra_wins\t0\tc1\t100\t60\t40=10X\t*\t0\t0\tAAAA\tIIII\tde:f:0.0450\n\
         mm2_wins\t0\tc1\t200\t60\t50=\t*\t0\t0\tCCCC\tIIII\n\
         tie\t0\tc1\t300\t60\t50=\t*\t0\t0\tGGGG\tIIII\n\
         secondary\t256\tc1\t400\t60\t10=\t*\t0\t0\tTTTT\tIIII\n"));
    let unindexed = write(&dir.0, "unindexed.sam", &format!(
        "{hdr}genomic\t0\tc1\t900\t60\t60=\t*\t0\t0\tTTTT\tIIII\tXA:Z:\tXC:Z:uLTRA_unindexed\n"));

    // uLTRA's view: scores 50, 30, 50 -- plus a read minimap2 never placed
    let ultra = write(&dir.0, "reads.sam", "\
@SQ\tSN:c1\tLN:1000\n\
ultra_wins\t0\tc1\t100\t60\t50=\t*\t0\t0\tAAAA\tIIII\tNM:i:0\n\
mm2_wins\t0\tc1\t200\t60\t40=10X\t*\t0\t0\tCCCC\tIIII\tNM:i:10\n\
tie\t0\tc1\t300\t60\t50=\t*\t0\t0\tGGGG\tIIII\tNM:i:0\n\
mm2_never_placed\t0\tc1\t500\t60\t70=\t*\t0\t0\tGGGG\tIIII\tNM:i:0\n");
    (dir, ultra, indexed, unindexed)
}

mod tempdir {
    pub struct Dir(pub std::path::PathBuf);
    impl Dir {
        pub fn new(tag: &str) -> Dir {
            let p = std::env::temp_dir()
                .join(format!("{tag}-{}-{:?}", std::process::id(), std::thread::current().id()));
            let _ = std::fs::remove_dir_all(&p);
            std::fs::create_dir_all(&p).unwrap();
            Dir(p)
        }
    }
    impl std::ops::Deref for Dir {
        type Target = std::path::Path;
        fn deref(&self) -> &std::path::Path { &self.0 }
    }
    impl Drop for Dir {
        fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); }
    }
}

fn run() -> (mm2::MergeStats, Vec<String>) {
    let (_d, ultra, indexed, unindexed) = fixture();
    let header = vec!["@SQ\tSN:c1\tLN:1000".to_string()];
    let st = mm2::output_final_alignments(&ultra, &indexed, &unindexed, &header).unwrap();
    let body = std::fs::read_to_string(&ultra).unwrap();
    (st, body.lines().map(|s| s.to_string()).collect())
}

#[test]
fn the_better_scoring_alignment_wins_and_a_tie_keeps_minimap2() {
    let (st, lines) = run();
    let field = |name: &str, n: usize| -> String {
        lines.iter().find(|l| l.starts_with(&format!("{name}\t")))
            .unwrap_or_else(|| panic!("{name} missing from the merged file"))
            .split('\t').nth(n).unwrap().to_string()
    };
    // uLTRA 50 vs minimap2 30 -> uLTRA, recognisable by its NM tag
    assert_eq!(field("ultra_wins", 5), "50=", "uLTRA scored better and must replace minimap2");
    // minimap2 50 vs uLTRA 30 -> minimap2
    assert_eq!(field("mm2_wins", 5), "50=");
    // equal scores keep minimap2: only a STRICTLY lower mm2 score hands it over
    assert_eq!(field("tie", 5), "50=");
    assert_eq!(st.equal_score, 1);
    assert_eq!(st.ultra_better, 1);
    assert_eq!(st.worse, 1, "mm2 50 vs uLTRA 30 is more than 10 better");
}

#[test]
fn unindexed_records_are_appended_and_counted() {
    let (st, lines) = run();
    let g = lines.iter().find(|l| l.starts_with("genomic\t")).expect("unindexed record dropped");
    assert!(g.contains("XC:Z:uLTRA_unindexed"), "the uLTRA_unindexed tag must survive the merge");
    assert_eq!(st.not_attempted, 1);
}

#[test]
fn a_read_minimap2_never_placed_keeps_its_ultra_alignment() {
    // Finding 39. The reference drops this read entirely.
    let (st, lines) = run();
    assert_eq!(st.mm2_absent, 1);
    assert!(lines.iter().any(|l| l.starts_with("mm2_never_placed\t")),
            "Finding 39: a read with no minimap2 record must keep uLTRA's alignment");
}

#[test]
fn every_read_appears_exactly_once() {
    let (_st, lines) = run();
    let mut names: Vec<&str> = lines.iter().filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').next().unwrap()).collect();
    names.sort();
    assert_eq!(names, vec!["genomic", "mm2_never_placed", "mm2_wins", "secondary", "tie", "ultra_wins"]);
}

#[test]
fn the_merged_file_carries_ultras_header_not_minimap2s() {
    let (_st, lines) = run();
    let hdr: Vec<&String> = lines.iter().take_while(|l| l.starts_with('@')).collect();
    assert_eq!(hdr, vec!["@SQ\tSN:c1\tLN:1000"],
               "the merged file must carry uLTRA's @SQ-only header, not minimap2's @HD/@PG");
}

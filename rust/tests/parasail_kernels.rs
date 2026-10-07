//! parasail offers three vectorisation strategies for the same objective.
//!
//! The port calls `sg_trace_scan_16` because `help_functions.parasail_alignment`
//! does. `striped` (Farrar) and `diag` compute the same optimal score by a
//! different traversal, so the choice is a pure performance question -- except
//! that a different traversal can report a different equally-optimal CIGAR,
//! which would reach reads.sam. Both are measured here.

use libparasail_sys as ps;
use serde_json::Value;
use std::ffi::CString;
use std::io::BufRead;
use std::os::raw::c_char;

#[path = "../src/parasail.rs"]
mod parasail;

unsafe fn run(f: &str, s1: &[u8], s2: &[u8]) -> (i32, String) {
    let alphabet = CString::new("ACGT").unwrap();
    let matrix = ps::parasail_matrix_create(alphabet.as_ptr(), 2, -2);
    let (p1, l1) = (s1.as_ptr() as *const c_char, s1.len() as i32);
    let (p2, l2) = (s2.as_ptr() as *const c_char, s2.len() as i32);
    let mut r = match f {
        "scan" => ps::parasail_sg_trace_scan_16(p1, l1, p2, l2, 3, 1, matrix),
        "striped" => ps::parasail_sg_trace_striped_16(p1, l1, p2, l2, 3, 1, matrix),
        _ => ps::parasail_sg_trace_diag_16(p1, l1, p2, l2, 3, 1, matrix),
    };
    if ps::parasail_result_is_saturated(r) != 0 {
        ps::parasail_result_free(r);
        r = match f {
            "scan" => ps::parasail_sg_trace_scan_32(p1, l1, p2, l2, 3, 1, matrix),
            "striped" => ps::parasail_sg_trace_striped_32(p1, l1, p2, l2, 3, 1, matrix),
            _ => ps::parasail_sg_trace_diag_32(p1, l1, p2, l2, 3, 1, matrix),
        };
    }
    let score = (*r).score;
    let cig = ps::parasail_result_get_cigar(r, p1, l1, p2, l2, matrix);
    let decoded = ps::parasail_cigar_decode(cig);
    let out = std::ffi::CStr::from_ptr(decoded).to_string_lossy().into_owned();
    ps::parasail_cigar_free(cig);
    ps::parasail_result_free(r);
    ps::parasail_matrix_free(matrix);
    (score, out)
}

#[test]
fn compare_parasail_vectorisation_strategies() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap().join("bench/oracle/aligner_calls.jsonl.gz");
    let f = flate2::read::GzDecoder::new(std::fs::File::open(&p).unwrap());
    let mut pairs: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    for line in std::io::BufReader::new(f).lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        if v["kind"] != "parasail_alignment" { continue; }
        let (a, b) = (v["s1"].as_str().unwrap(), v["s2"].as_str().unwrap());
        if !a.is_empty() && !b.is_empty() { pairs.push((a.as_bytes().to_vec(), b.as_bytes().to_vec())); }
    }
    println!("parasail strategies over {} recorded calls:", pairs.len());
    let base: Vec<(i32, String)> = pairs.iter().map(|(a, b)| unsafe { run("scan", a, b) }).collect();
    for k in ["scan", "striped", "diag"] {
        let t = std::time::Instant::now();
        let mut got = Vec::with_capacity(pairs.len());
        for _ in 0..5 {
            got.clear();
            for (a, b) in &pairs { got.push(unsafe { run(k, a, b) }); }
        }
        let el = t.elapsed();
        let same_s = got.iter().zip(&base).filter(|(a, b)| a.0 == b.0).count();
        let same_c = got.iter().zip(&base).filter(|(a, b)| a.1 == b.1).count();
        println!("   {k:8} {el:>12.2?}   same score {same_s}/{}   same CIGAR {same_c}/{}",
                 pairs.len(), pairs.len());
    }
}

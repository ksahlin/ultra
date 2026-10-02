//! parasail, via `libparasail-sys` (which bundles the C sources).
//!
//! uLTRA uses exactly one parasail entry point, `help_functions.parasail_alignment`:
//!
//!     matrix = parasail.matrix_create("ACGT", 2, -2)
//!     r = parasail.sg_trace_scan_16(s1, s2, 3, 1, matrix)
//!     if r.saturated: r = parasail.sg_trace_scan_32(...)
//!
//! It is the DEFAULT alignment path in `get_exact_alignment` -- edlib is only
//! used for sequences over 20 kb -- so although it is a minority of aligner
//! CALLS it produces most of the final CIGARs.

use libparasail_sys as ps;
use std::ffi::CString;
// c_char is signed on x86_64 but UNSIGNED on aarch64-linux, so a hardcoded
// `*const i8` builds on one and not the other.
use std::os::raw::c_char;

pub struct Result_ {
    pub score: i32,
    pub cigar: String,
}

/// `parasail_alignment(s1, s2)` with the reference's fixed parameters.
pub fn sg_trace(s1: &[u8], s2: &[u8], match_score: i32, mismatch_penalty: i32, open: i32, ext: i32) -> Result_ {
    unsafe {
        let alphabet = CString::new("ACGT").unwrap();
        let matrix = ps::parasail_matrix_create(alphabet.as_ptr(), match_score, mismatch_penalty);
        let mut r = ps::parasail_sg_trace_scan_16(
            s1.as_ptr() as *const c_char, s1.len() as i32,
            s2.as_ptr() as *const c_char, s2.len() as i32,
            open, ext, matrix,
        );
        // The reference tests `result.saturated`; parasail exposes it as a
        // predicate. The 16-bit kernel saturates on long or divergent pairs and
        // the reference then recomputes with the 32-bit one.
        if ps::parasail_result_is_saturated(r) != 0 {
            ps::parasail_result_free(r);
            r = ps::parasail_sg_trace_scan_32(
                s1.as_ptr() as *const c_char, s1.len() as i32,
                s2.as_ptr() as *const c_char, s2.len() as i32,
                open, ext, matrix,
            );
        }
        let score = ps::parasail_result_get_score(r);
        let cig = ps::parasail_result_get_cigar(r, s1.as_ptr() as *const c_char, s1.len() as i32,
                                                s2.as_ptr() as *const c_char, s2.len() as i32, matrix);
        let cs = ps::parasail_cigar_decode(cig);
        let cigar = std::ffi::CStr::from_ptr(cs).to_string_lossy().into_owned();
        ps::parasail_cigar_free(cig);
        ps::parasail_result_free(r);
        ps::parasail_matrix_free(matrix);
        Result_ { score, cigar }
    }
}

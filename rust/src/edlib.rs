//! The vendored edlib, linked in.
//!
//! uLTRA reaches edlib through two wrappers, and both call
//! `edlib.align(query, target, task=..., mode=..., k=...)`. Measured on
//! sirv-10k, every live call is `mode="HW", task="path"` -- 11 810 of 12 000
//! aligner calls in total.
//!
//! `help_functions.edlib_alignment` reads `result['locations'][0]`, and which
//! location edlib puts first when several are optimal is not specified
//! anywhere. Linking edlib's own implementation is what makes that a non-issue.

use std::os::raw::{c_char, c_int};

#[repr(C)]
#[derive(Clone, Copy)]
struct EdlibAlignConfig {
    k: c_int,
    mode: c_int,
    task: c_int,
    additional_equalities: *const EdlibEqualityPair,
    additional_equalities_length: c_int,
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct EdlibEqualityPair {
    pub first: c_char,
    pub second: c_char,
}

#[repr(C)]
struct EdlibAlignResult {
    status: c_int,
    edit_distance: c_int,
    end_locations: *mut c_int,
    start_locations: *mut c_int,
    num_locations: c_int,
    alignment: *mut u8,
    alignment_length: c_int,
    alphabet_length: c_int,
}

extern "C" {
    fn edlibAlign(
        query: *const c_char, query_length: c_int,
        target: *const c_char, target_length: c_int,
        config: EdlibAlignConfig,
    ) -> EdlibAlignResult;
    fn edlibFreeAlignResult(result: EdlibAlignResult);
    fn edlibAlignmentToCigar(alignment: *const u8, alignment_length: c_int, cigar_format: c_int) -> *mut c_char;
}

pub const MODE_NW: c_int = 0;
pub const MODE_SHW: c_int = 1;
pub const MODE_HW: c_int = 2;
pub const TASK_DISTANCE: c_int = 0;
pub const TASK_LOC: c_int = 1;
pub const TASK_PATH: c_int = 2;
/// edlib's EDLIB_CIGAR_EXTENDED, which is what python-edlib returns (`=`/`X`).
const CIGAR_EXTENDED: c_int = 1;

#[derive(Debug, Clone, PartialEq)]
pub struct Alignment {
    pub edit_distance: i32,
    /// (start, end) pairs. `start` is -1 when edlib did not compute starts.
    pub locations: Vec<(i64, i64)>,
    pub cigar: Option<String>,
}

/// `edlib.align(query, target, mode=, task=, k=)`.
///
/// NOTE on `k`: the reference passes a FLOAT here
/// (`k = 0.4*min(len(read_seq), len(exon_seq))`), which the Python binding
/// coerces. We take an i32 and the caller does the same truncation.
pub fn align(query: &[u8], target: &[u8], mode: c_int, task: c_int, k: i32) -> Alignment {
    let cfg = EdlibAlignConfig {
        k: k as c_int,
        mode,
        task,
        additional_equalities: std::ptr::null(),
        additional_equalities_length: 0,
    };
    unsafe {
        let r = edlibAlign(
            query.as_ptr() as *const c_char, query.len() as c_int,
            target.as_ptr() as *const c_char, target.len() as c_int,
            cfg,
        );
        let mut locations = Vec::new();
        for i in 0..r.num_locations as isize {
            let end = *r.end_locations.offset(i) as i64;
            let start = if r.start_locations.is_null() { -1 } else { *r.start_locations.offset(i) as i64 };
            locations.push((start, end));
        }
        let cigar = if !r.alignment.is_null() && r.alignment_length > 0 {
            let c = edlibAlignmentToCigar(r.alignment, r.alignment_length, CIGAR_EXTENDED);
            if c.is_null() {
                None
            } else {
                let s = std::ffi::CStr::from_ptr(c).to_string_lossy().into_owned();
                libc_free(c as *mut std::ffi::c_void);
                Some(s)
            }
        } else {
            None
        };
        let out = Alignment { edit_distance: r.edit_distance as i32, locations, cigar };
        edlibFreeAlignResult(r);
        out
    }
}

extern "C" {
    #[link_name = "free"]
    fn libc_free(p: *mut std::ffi::c_void);
}

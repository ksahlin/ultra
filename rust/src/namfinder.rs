//! The vendored namfinder, linked in.
//!
//! uLTRA's seeding step shells out to a `namfinder` binary and pipes its stdout
//! through gzip (`seed_wrapper.py:74`, via `os.system` with an unquoted shell
//! string -- PORTING.md Finding 8). The port calls the same code in-process.
//!
//! `run` takes the argument vector namfinder's own CLI takes, because it IS
//! namfinder's CLI: see vendor/namfinder/shim.cpp.

use std::ffi::CString;
use std::os::raw::{c_char, c_int};

extern "C" {
    fn nf_run(argc: c_int, argv: *mut *mut c_char) -> c_int;
}

/// Invoke namfinder with `args` (argv[0] included). Returns its exit status.
///
/// NOTE this writes NAMs to stdout and logging to stderr, exactly as the
/// binary does; the caller redirects. Reproducing uLTRA's behaviour means the
/// output must be gzipped at level 1 into `seeds.txt.gz`.
pub fn run(args: &[String]) -> i32 {
    let cstrings: Vec<CString> = args
        .iter()
        .map(|a| CString::new(a.as_str()).expect("argument contains a NUL byte"))
        .collect();
    let mut ptrs: Vec<*mut c_char> = cstrings.iter().map(|c| c.as_ptr() as *mut c_char).collect();
    ptrs.push(std::ptr::null_mut());
    unsafe { nf_run((ptrs.len() - 1) as c_int, ptrs.as_mut_ptr()) as i32 }
}

/// Build namfinder's argument vector exactly as `find_nams_namfinder` does,
/// including its thinning arithmetic. Reproduced from seed_wrapper.py:43-58.
///
///   thinning 0 -> s = strobe_size,     l = strobe_size,      u = strobe_size + 1
///   thinning 1 -> s = strobe_size - 2, l = (strobe_size+1)/3, u = (strobe_size+1)/3 + 1
///   thinning 2 -> s = strobe_size - 4, l = (strobe_size+1)/5, u = (strobe_size+1)/5 + 1
///
/// The divisions are Python's `//` on positive ints, i.e. truncating.
pub fn argv_for(refs: &str, reads: &str, threads: i64, strobe_size: i64, thinning: i64) -> Vec<String> {
    let (s, l, u) = match thinning {
        1 => (strobe_size - 2, (strobe_size + 1) / 3, (strobe_size + 1) / 3 + 1),
        2 => (strobe_size - 4, (strobe_size + 1) / 5, (strobe_size + 1) / 5 + 1),
        _ => (strobe_size, strobe_size, strobe_size + 1),
    };
    [
        "namfinder", "-k", &strobe_size.to_string(), "-s", &s.to_string(),
        "-l", &l.to_string(), "-u", &u.to_string(),
        "-C", "500", "-L", "1000", "-t", &threads.to_string(), "-S", refs, reads,
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

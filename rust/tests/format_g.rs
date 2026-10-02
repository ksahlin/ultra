//! `format_g` must agree with C's `%g`, which is how pysam spells a float tag.
//!
//! Compared against the platform's own `printf` rather than a table of
//! expected strings, so the test cannot drift from the thing it is imitating.

#[path = "../src/samfmt.rs"]
mod samfmt;

extern "C" {
    fn snprintf(s: *mut u8, n: usize, fmt: *const u8, ...) -> i32;
}

fn c_printf_g(v: f64) -> String {
    let mut buf = [0u8; 64];
    let n = unsafe { snprintf(buf.as_mut_ptr(), 64, b"%g\0".as_ptr(), v) };
    String::from_utf8_lossy(&buf[..n as usize]).into_owned()
}

#[test]
fn matches_c_printf_g() {
    // the shapes that actually decide a %g branch: the %e/%f switchover at
    // 1e-5 and 1e6, rounding that carries into a new exponent, and zero
    let mut vals: Vec<f64> = vec![
        0.0, -0.0, 0.045, 0.007, 1.0, 123456.0, 1234567.0, 1e-4, 1e-5, 1e-7,
        123456789.0, 0.1, 0.999999, 0.9999995, 3.14159265, 1e20, 5e-324,
        0.000123456789, 99999.5, 100000.5, 2.5, 1e6, -0.0450, -1234567.0,
    ];
    // plus a deterministic spread across the exponent range
    let mut seed: u64 = 0x2545F4914F6CDD1D;
    let mut next = || {
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        seed
    };
    for _ in 0..4000 {
        let m = (next() % 2_000_000) as f64 / 1_000_000.0 - 1.0;
        let e = (next() % 25) as i32 - 12;
        vals.push(m * 10f64.powi(e));
    }

    let bad: Vec<String> = vals.iter()
        .filter_map(|&v| {
            let (got, want) = (samfmt::format_g(v), c_printf_g(v));
            (got != want).then(|| format!("{v:e}: got {got:?}, C says {want:?}"))
        })
        .collect();
    assert!(bad.is_empty(), "{} of {} mismatched:\n{}",
            bad.len(), vals.len(), bad.iter().take(10).cloned().collect::<Vec<_>>().join("\n"));
    println!("format_g agrees with C's %g on {} values", vals.len());
}

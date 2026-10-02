//! How pysam spells a SAM tag when it rewrites a record.
//!
//! Kept separate from `mm2.rs` so `tests/format_g.rs` can include it with no
//! dependencies and check it against the platform's own `printf`.

/// C's `%g` with the default precision of 6, which is how pysam writes a
/// float tag.
///
/// This is not cosmetic. minimap2 emits `de:f:0.0450`; pysam parses the tag
/// and writes it back as `de:f:0.045`, so every SAM the reference produces by
/// round-tripping minimap2's records carries pysam's spelling, not minimap2's.
/// Copying the record verbatim -- the obvious thing for a port to do, and what
/// this one did first -- diverges on 704 of 10 000 records on sirv-10k.
pub fn format_g(v: f64) -> String {
    const P: i32 = 6;
    if !v.is_finite() {
        return format!("{v}");
    }
    // the decimal exponent AFTER rounding to P significant digits, which is
    // what C's %g keys off
    let e = format!("{:.*e}", (P - 1) as usize, v);
    let x: i32 = e[e.find('e').unwrap() + 1..].parse().unwrap_or(0);

    let out = if x < -4 || x >= P {
        let mant = &e[..e.find('e').unwrap()];
        let mant = strip_zeros(mant);
        format!("{}e{}{:02}", mant, if x < 0 { '-' } else { '+' }, x.abs())
    } else {
        strip_zeros(&format!("{:.*}", (P - 1 - x) as usize, v))
    };
    // C prints "-0" for negative zero and so does Python's %g, so pysam does
    // too -- do not "tidy" it away.
    out
}

fn strip_zeros(s: &str) -> String {
    if !s.contains('.') {
        return s.to_string();
    }
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

//! uLTRA — splice alignment of long transcriptomic reads.
//!
//! PORT STATUS: Stage 1. The CLI contract is implemented and checked against
//! the recorded goldens. Everything past argument validation is not yet
//! written and exits **70** (EX_SOFTWARE), which is a code the reference never
//! produces — so a case that reaches unimplemented territory can never be
//! mistaken for a pass.

mod cli;
mod dump;
mod fasta;
mod gtf;
mod index;
mod text;

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();

    // `dump-index` is a PORT-ONLY subcommand used by the equivalence harness:
    // it builds the index from the inputs and writes the canonical text
    // rendering that bench/dump_reference.py produces from the reference's
    // pickles. It is handled before cli::parse so it cannot perturb the
    // argparse emulation -- the reference has no such subcommand, and the
    // help text and the invalid-choice message must stay exactly as recorded.
    //
    //   uLTRA dump-index <ref.fa> <annot.gtf> <outdir> [--flank_size N]
    //                    [--small_exon_threshold N] [--min_segm N]
    if argv.first().map(|s| s.as_str()) == Some("dump-index") {
        return dump_index(&argv[1..]);
    }

    match cli::parse(&argv) {
        cli::Outcome::Stdout(s) => {
            print!("{s}");
            let _ = std::io::stdout().flush();
            ExitCode::from(0)
        }
        cli::Outcome::UsageError { usage, msg } => {
            eprint!("{usage}");
            eprintln!("{msg}");
            ExitCode::from(2)
        }
        cli::Outcome::Run(args) => {
            // The reference creates the output folder here, AFTER parsing and
            // before any validation of the inputs. Finding 17: argparse
            // failures leave nothing on disk, everything later leaves an empty
            // directory behind. Reproduced deliberately.
            //
            // help_functions.mkdir_p prints "creating <path>" -- but only when
            // it actually creates it, because it swallows EEXIST silently. So
            // the message is conditional on the directory not already existing.
            if !Path::new(&args.outfolder).is_dir() {
                if let Err(e) = std::fs::create_dir_all(&args.outfolder) {
                    eprintln!("uLTRA: could not create output folder {}: {e}", args.outfolder);
                    return ExitCode::from(1);
                }
                println!("creating {}", args.outfolder);
            }

            // Finding 16: the reference prints this and exits 0. Reproduced
            // for now; changing it to exit 1 is a deliberate divergence that
            // needs its own commit and its own golden.
            if args.thinning < 0 || args.thinning > 2 {
                println!("Invalid thinning level. Choose 0, 1 or 2.");
                return ExitCode::from(0);
            }

            // align_reads' first act. Finding 16 again: a detected, clearly
            // reported user error that exits 0. The double space after
            // "specified: " and the "forder" typo are the reference's and are
            // contract until deliberately changed.
            if matches!(args.sub, cli::Sub::Align) && !args.index.is_empty() {
                if !Path::new(&args.index).is_dir() {
                    println!(
                        "The index folder specified for alignment is not found. You specified:  {}",
                        args.index
                    );
                    println!(
                        "Build  the index to this folder, or specify another forder where the index has been built."
                    );
                    return ExitCode::from(0);
                }
            }

            eprintln!(
                "uLTRA: {} is not implemented in the Rust port yet (stage 1: CLI only)",
                args.sub.name()
            );
            ExitCode::from(70)
        }
    }
}

/// Build the index from inputs and render it for the stage oracle.
fn dump_index(args: &[String]) -> ExitCode {
    let mut pos: Vec<&str> = Vec::new();
    let mut flank_size: u64 = 1000;
    let mut small_exon_threshold: u64 = 200;
    let mut min_segm: u64 = 25;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--flank_size" => { i += 1; flank_size = args[i].parse().unwrap_or(flank_size); }
            "--small_exon_threshold" => { i += 1; small_exon_threshold = args[i].parse().unwrap_or(small_exon_threshold); }
            "--min_segm" => { i += 1; min_segm = args[i].parse().unwrap_or(min_segm); }
            "--disable_infer" => {}   // the port never infers; see Finding 23
            other => pos.push(other),
        }
        i += 1;
    }
    if pos.len() != 3 {
        eprintln!("usage: uLTRA dump-index <ref.fa> <annot.gtf> <outdir> [--flank_size N] [--small_exon_threshold N] [--min_segm N]");
        return ExitCode::from(2);
    }
    let (ref_path, gtf_path, out) = (pos[0], pos[1], pos[2]);

    let refs = match fasta::read(ref_path) {
        Ok(r) => r,
        Err(e) => { eprintln!("uLTRA: cannot read reference {ref_path}: {e}"); return ExitCode::from(1); }
    };
    let mut refs_lengths: std::collections::BTreeMap<String, u64> = Default::default();
    for r in &refs {
        refs_lengths.insert(r.name.clone(), r.seq.len() as u64);
    }

    let parsed = match gtf::parse(gtf_path) {
        Ok(g) => g,
        Err(e) => { eprintln!("uLTRA: {e}"); return ExitCode::from(1); }
    };

    let mut ix = index::build(&parsed, &refs_lengths, &index::Params {
        flank_size, small_exon_threshold, min_segm,
    });

    // prep_seqs: extract the sequences, keyed by chr_id for contigs present in
    // BOTH the fasta and the annotation.
    let mut by_id: std::collections::HashMap<u64, &[u8]> = Default::default();
    for r in &refs {
        if let Some(&cid) = ix.chr_to_id.get(&r.name) {
            by_id.insert(cid, &r.seq);
            ix.refs_id_lengths.insert(cid, r.seq.len() as u64);
        }
    }
    let grab = |keys: Vec<index::Key>, by_id: &std::collections::HashMap<u64, &[u8]>|
        -> std::collections::BTreeMap<index::Key, String> {
        let mut out = std::collections::BTreeMap::new();
        for key in keys {
            if let Some(seq) = by_id.get(&key.0) {
                let a = std::cmp::min(key.1 as usize, seq.len());
                let b = std::cmp::min(key.2 as usize, seq.len());
                if a <= b {
                    out.insert(key, String::from_utf8_lossy(&seq[a..b]).into_owned());
                }
            }
        }
        out
    };
    ix.ref_segment_sequences = grab(ix.segment_id_to_choordinates.keys().cloned().collect(), &by_id);
    ix.ref_exon_sequences = grab(ix.exon_choordinates_to_id.iter().cloned().collect(), &by_id);
    ix.ref_flank_sequences = grab(ix.flank_choordinates.iter().cloned().collect(), &by_id);
    // ref_part_sequences is the PART sequences updated with the FLANK ones --
    // uLTRA does `update_nested(ref_part_sequences, ref_flank_sequences)`.
    let mut parts = grab(ix.parts_to_segments.keys().cloned().collect(), &by_id);
    for (k, v) in ix.ref_flank_sequences.iter() {
        parts.insert(*k, v.clone());
    }
    ix.ref_part_sequences = parts;

    match dump::dump(&ix, std::path::Path::new(out)) {
        Ok(w) => {
            for (name, n) in w {
                println!("  {name:34} {n:8} lines");
            }
            ExitCode::from(0)
        }
        Err(e) => { eprintln!("uLTRA: cannot write rendering to {out}: {e}"); ExitCode::from(1) }
    }
}

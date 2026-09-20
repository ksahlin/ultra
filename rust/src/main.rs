//! uLTRA — splice alignment of long transcriptomic reads.
//!
//! PORT STATUS: Stage 1. The CLI contract is implemented and checked against
//! the recorded goldens. Everything past argument validation is not yet
//! written and exits **70** (EX_SOFTWARE), which is a code the reference never
//! produces — so a case that reaches unimplemented territory can never be
//! mistaken for a pass.

mod cli;
mod text;

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();

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

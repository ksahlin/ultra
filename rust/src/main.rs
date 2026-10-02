//! uLTRA — splice alignment of long transcriptomic reads.
//!
//! PORT STATUS: Stage 1. The CLI contract is implemented and checked against
//! the recorded goldens. Everything past argument validation is not yet
//! written and exits **70** (EX_SOFTWARE), which is a code the reference never
//! produces — so a case that reaches unimplemented territory can never be
//! mistaken for a pass.

mod cli;
mod aligndriver;
mod colinear;
mod driver;
mod dump;
mod edlib;
mod parasail;
mod fasta;
mod gtf;
mod index;
mod indexio;
mod mam;
mod namfinder;
mod prefilter;
mod samfmt;
mod mm2;
mod reads;
mod samout;
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

    // Port-only: emit the preprocessed reads uLTRA feeds to namfinder, i.e.
    // the contents of reads_tmp.fa.gz, so the harness can diff them against
    // the reference's.
    //   uLTRA dump-reads <reads.[fa|fq]> [--reduce_read_ployA N] [--strict]
    if argv.first().map(|s| s.as_str()) == Some("dump-reads") {
        return dump_reads(&argv[1..]);
    }

    // Port-only: summarise parsed seed records, to diff against the
    // reference's read_seeds on real output.
    if argv.first().map(|s| s.as_str()) == Some("dump-seeds") {
        let f = std::fs::File::open(&argv[1]).expect("open seeds");
        let gz = flate2::read::GzDecoder::new(f);
        let recs = reads::read_seeds(std::io::BufReader::new(gz)).expect("parse");
        for r in &recs {
            println!("{}\t{}\t{}\t{}", r.acc, r.hits.len(), r.acc_rev, r.hits_rc.len());
        }
        return ExitCode::from(0);
    }

    // Also port-only: a passthrough to the LINKED namfinder, so the harness can
    // diff it against the upstream binary. Not part of the CLI contract.
    if argv.first().map(|s| s.as_str()) == Some("namfinder") {
        let mut a = vec!["namfinder".to_string()];
        a.extend(argv[1..].iter().cloned());
        return ExitCode::from(namfinder::run(&a) as u8);
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

            match args.sub {
                cli::Sub::Index => run_index(&args),
                cli::Sub::Align => run_align(&args),
                cli::Sub::Pipeline => match run_index(&args) {
                    c if c == ExitCode::from(0) => run_align(&args),
                    c => c,
                },
            }
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

/// Reproduce the `reads_tmp.fa.gz` payload: every read, polyA-compressed,
/// written as `>{acc}\n{seq}\n`. uLTRA gzips this at level 1; the harness
/// compares decompressed, so this writes plain text.
fn dump_reads(args: &[String]) -> ExitCode {
    let mut path: Option<&str> = None;
    let mut polya: usize = 8;
    let mut strict = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--reduce_read_ployA" => { i += 1; polya = args[i].parse().unwrap_or(polya); }
            "--strict" => strict = true,
            other => path = Some(other),
        }
        i += 1;
    }
    let path = match path {
        Some(p) => p,
        None => { eprintln!("usage: uLTRA dump-reads <reads.[fa|fq]> [--reduce_read_ployA N] [--strict]"); return ExitCode::from(2); }
    };
    let f = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(e) => { eprintln!("uLTRA: cannot read {path}: {e}"); return ExitCode::from(1); }
    };
    let recs = match reads::readfq(std::io::BufReader::new(f), strict) {
        Ok(r) => r,
        Err(e) => { eprintln!("uLTRA: {e}"); return ExitCode::from(1); }
    };
    let stdout = std::io::stdout();
    let mut w = std::io::BufWriter::new(stdout.lock());
    for r in &recs {
        // remove_read_polyA_ends(seq, qual, args.reduce_read_ployA, 5)
        let s = reads::remove_read_polya_ends(&r.seq, polya, 5);
        if writeln!(w, ">{}", r.name).is_err() || writeln!(w, "{s}").is_err() {
            return ExitCode::from(1);
        }
    }
    let _ = w.flush();
    ExitCode::from(0)
}

/// `uLTRA index`: build the annotation index and write it as one binary file.
fn run_index(args: &cli::Args) -> ExitCode {
    let refs = match fasta::read(&args.reference) {
        Ok(r) => r,
        Err(e) => { eprintln!("uLTRA: cannot read reference {}: {e}", args.reference); return ExitCode::from(1); }
    };
    let mut refs_lengths: std::collections::BTreeMap<String, u64> = Default::default();
    let mut ref_order: Vec<String> = Vec::new();
    for r in &refs {
        refs_lengths.insert(r.name.clone(), r.seq.len() as u64);
        ref_order.push(r.name.clone());
    }

    let parsed = match gtf::parse(&args.gtf) {
        Ok(g) => g,
        Err(e) => { eprintln!("uLTRA: {e}"); return ExitCode::from(1); }
    };
    let mut ix = index::build(&parsed, &refs_lengths, &index::Params {
        flank_size: args.flank_size as u64,
        small_exon_threshold: args.small_exon_threshold as u64,
        min_segm: args.min_segm as u64,
    });
    ix.ref_order = ref_order;
    let mut by_id: std::collections::HashMap<u64, &[u8]> = Default::default();
    for r in &refs {
        if let Some(&cid) = ix.chr_to_id.get(&r.name) {
            by_id.insert(cid, &r.seq);
            ix.refs_id_lengths.insert(cid, r.seq.len() as u64);
        }
    }
    let grab = |keys: Vec<index::Key>| -> std::collections::BTreeMap<index::Key, String> {
        let mut out = std::collections::BTreeMap::new();
        for key in keys {
            if let Some(seq) = by_id.get(&key.0) {
                let a = std::cmp::min(key.1 as usize, seq.len());
                let b = std::cmp::min(key.2 as usize, seq.len());
                if a <= b { out.insert(key, String::from_utf8_lossy(&seq[a..b]).into_owned()); }
            }
        }
        out
    };
    ix.ref_segment_sequences = grab(ix.segment_id_to_choordinates.keys().cloned().collect());
    ix.ref_exon_sequences = grab(ix.exon_choordinates_to_id.iter().cloned().collect());
    ix.ref_flank_sequences = grab(ix.flank_choordinates.iter().cloned().collect());
    let mut parts = grab(ix.parts_to_segments.keys().cloned().collect());
    for (k, v) in ix.ref_flank_sequences.iter() { parts.insert(*k, v.clone()); }
    ix.ref_part_sequences = parts;

    let dest = index_dir(args).join("ultra.idx");
    match indexio::write(&ix, &dest) {
        Ok(()) => { println!("wrote {}", dest.display()); ExitCode::from(0) }
        Err(e) => { eprintln!("uLTRA: cannot write {}: {e}", dest.display()); ExitCode::from(1) }
    }
}

fn index_dir(args: &cli::Args) -> std::path::PathBuf {
    if args.index.is_empty() { std::path::PathBuf::from(&args.outfolder) }
    else { std::path::PathBuf::from(&args.index) }
}

/// `uLTRA align`.
fn run_align(args: &cli::Args) -> ExitCode {
    let idx_path = index_dir(args).join("ultra.idx");
    let ix = match indexio::read(&idx_path) {
        Ok(i) => i,
        Err(e) => { eprintln!("uLTRA: cannot read index {}: {e}", idx_path.display()); return ExitCode::from(1); }
    };
    let out = std::path::PathBuf::from(&args.outfolder);

    // 0. the genomic prefilter. minimap2 is run over the whole genome and the
    //    reads it places outside uLTRA's indexed regions are kept as minimap2
    //    aligned them, never reaching uLTRA's aligner.
    let mut reads_src = args.reads.clone();
    let mut split: Option<(std::path::PathBuf, std::path::PathBuf)> = None;
    if !args.disable_mm2 {
        println!("Filtering reads aligned to unindexed regions with minimap2 ");
        let keys: Vec<index::Key> = ix.ref_part_sequences.keys().cloned().collect();
        let regions = prefilter::IndexedRegions::from_parts(&keys, &ix.id_to_chr);
        let sam = match mm2::run_minimap2(
            &args.reference, &args.reads, &out, args.nr_cores, args.mm2_ksize) {
            Ok(p) => p,
            Err(e) => { eprintln!("uLTRA: {e}"); return ExitCode::from(1); }
        };
        let o = match mm2::filter_reads_to_align(&sam, &regions, &out, args.genomic_frac) {
            Ok(o) => o,
            Err(e) => { eprintln!("uLTRA: {e}"); return ExitCode::from(1); }
        };
        println!("Done filtering. Reads filtered:{}", o.nr_unindexed);
        reads_src = o.reads_path.to_string_lossy().into_owned();
        split = Some((o.indexed_path, o.unindexed_path));
    }

    // 1. refs_sequences.fa -- the part and flank sequences namfinder indexes,
    //    named `chr^start^stop`. The reference writes these in dict order; the
    //    port writes them sorted (PORTING.md Finding 27 established that the
    //    reference's order is hash-dependent and reaches results, and the fix
    //    branch sorts them too).
    let refs_path = out.join("refs_sequences.fa");
    {
        let f = match std::fs::File::create(&refs_path) {
            Ok(f) => f,
            Err(e) => { eprintln!("uLTRA: cannot write {}: {e}", refs_path.display()); return ExitCode::from(1); }
        };
        let mut w = std::io::BufWriter::new(f);
        for (k, seq) in ix.ref_part_sequences.iter() {
            if writeln!(w, ">{}^{}^{}", k.0, k.1, k.2).is_err() || writeln!(w, "{seq}").is_err() {
                eprintln!("uLTRA: cannot write {}", refs_path.display());
                return ExitCode::from(1);
            }
        }
        let _ = w.flush();
    }

    // 2. the reads namfinder sees: polyA-compressed, as uLTRA does before
    //    seeding. Written uncompressed; namfinder reads either.
    let reads_tmp = out.join("reads_tmp.fa");
    let recs = {
        let f = match std::fs::File::open(&reads_src) {
            Ok(f) => f,
            Err(e) => { eprintln!("uLTRA: cannot read {reads_src}: {e}"); return ExitCode::from(1); }
        };
        match reads::readfq(std::io::BufReader::new(f), false) {
            Ok(r) => r,
            Err(e) => { eprintln!("uLTRA: {e}"); return ExitCode::from(1); }
        }
    };
    {
        let f = std::fs::File::create(&reads_tmp).expect("create reads_tmp");
        let mut w = std::io::BufWriter::new(f);
        for r in &recs {
            let s = reads::remove_read_polya_ends(&r.seq, args.reduce_read_polya as usize, 5);
            let _ = writeln!(w, ">{}", r.name);
            let _ = writeln!(w, "{s}");
        }
        let _ = w.flush();
    }

    // 3. seeds, from the LINKED namfinder (Finding 26) -- nothing required on
    //    PATH. It runs in a forked child all the same: namfinder's index is by
    //    far the largest allocation in the run (1 042 MB of a 1 196 MB peak on
    //    sirv-10k) and it does not return that memory to the OS when it
    //    finishes, so calling it in-process would leave the whole alignment
    //    phase running on top of it. Forking gets the subprocess's memory
    //    behaviour back without reintroducing the subprocess's dependency.
    //    See Finding 41.
    let seeds_path = out.join("seeds.txt");
    {
        let argv = namfinder::argv_for(
            refs_path.to_str().unwrap(), reads_tmp.to_str().unwrap(),
            args.nr_cores, args.s, args.thinning);
        let fd = match std::fs::File::create(&seeds_path) {
            Ok(f) => { use std::os::unix::io::IntoRawFd; f.into_raw_fd() }
            Err(e) => { eprintln!("uLTRA: cannot write {}: {e}", seeds_path.display()); return ExitCode::from(1); }
        };
        // The parent is still single-threaded here, and the child only calls
        // into C, so there is nothing for fork to tear in half. Flush first or
        // anything still buffered is written twice, once by each process.
        let _ = std::io::stdout().flush();
        let _ = std::io::stderr().flush();
        let pid = unsafe { libc_fork() };
        if pid == 0 {
            unsafe {
                libc_dup2(fd, 1);
                let rc = namfinder::run(&argv);
                // _exit, not exit: the child must not run atexit handlers or
                // flush streams the parent still owns.
                libc__exit(rc);
            }
            unreachable!();
        }
        unsafe { libc_close(fd) };
        if pid < 0 {
            eprintln!("uLTRA: could not fork for namfinder");
            return ExitCode::from(1);
        }
        let mut status: i32 = 0;
        let rc = unsafe { libc_waitpid(pid, &mut status, 0) };
        // WIFEXITED / WEXITSTATUS, which are macros and so not in the ABI
        let exited = status & 0x7f == 0;
        let code = (status >> 8) & 0xff;
        if rc < 0 || !exited || code != 0 {
            eprintln!("uLTRA: namfinder failed (status {status})");
            return ExitCode::from(1);
        }
    }

    // 4. align, one read at a time over a SHARED index (Finding 4: the
    //    reference loads a copy per worker instead).
    let seed_recs = {
        let f = std::fs::File::open(&seeds_path).expect("open seeds");
        match reads::read_seeds(std::io::BufReader::new(f)) {
            Ok(v) => v,
            Err(e) => { eprintln!("uLTRA: {e}"); return ExitCode::from(1); }
        }
    };
    let by_name: std::collections::HashMap<&str, &reads::SeedRecord> =
        seed_recs.iter().map(|s| (s.acc.as_str(), s)).collect();

    let p = driver::Params {
        max_intron: args.max_intron,
        min_acc: args.min_acc,
        dropoff: args.dropoff,
        max_loc: args.max_loc,
        alignment_threshold: args.alignment_threshold,
        non_covered_cutoff: args.non_covered_cutoff,
        reduce_read_polya: args.reduce_read_polya as usize,
    };

    // @SQ only -- no @PG, no @HD, matching the reference exactly (Part 5)
    let header: Vec<String> = ix.ref_order.iter()
        .map(|name| format!("@SQ\tSN:{}\tLN:{}", name,
                            ix.refs_lengths.get(name).copied().unwrap_or(0)))
        .collect();

    let sam_path = out.join(format!("{}.sam", args.prefix));
    {
        let f = std::fs::File::create(&sam_path).expect("create sam");
        let mut w = std::io::BufWriter::new(f);
        for h in &header {
            let _ = writeln!(w, "{h}");
        }
        align_all(&ix, &recs, &by_name, &p, args.nr_cores, &mut w);
        let _ = w.flush();
    }

    // 5. merge uLTRA's alignments with minimap2's, keeping whichever scored
    //    better per read.
    if let Some((indexed, unindexed)) = split {
        match mm2::output_final_alignments(&sam_path, &indexed, &unindexed, &header) {
            Ok(st) => {
                println!("Total mm2's primary alignments replaced with uLTRA: {}", st.ultra_better);
                println!("{} primary alignments had equal score with alternative aligner.", st.equal_score);
                println!("{} primary alignments had slightly better score with alternative aligner.", st.slightly_worse);
                println!("{} primary alignments had significantly better score with alternative aligner.", st.worse);
                println!("{} reads were unmapped with ultra but not by alternative aligner.", st.ultra_unmapped);
                println!("{} reads were not attempted to be aligned with ultra (unindexed regions), instead alternative aligner was used.", st.not_attempted);
                if st.mm2_absent > 0 {
                    println!("{} reads had no minimap2 record at all; uLTRA's alignments are kept rather than dropped (PORTING.md Finding 39).", st.mm2_absent);
                }
            }
            Err(e) => { eprintln!("uLTRA: merging alignments failed: {e}"); return ExitCode::from(1); }
        }
    }

    if !args.keep_temporary_files {
        let _ = std::fs::remove_file(out.join("refs_sequences.fa"));
        if !args.disable_mm2 {
            for f in ["minimap2.sam", "minimap2_errors.1",
                      "reads_after_genomic_filtering.fastq", "indexed.sam", "unindexed.sam"] {
                let _ = std::fs::remove_file(out.join(f));
            }
        }
        // The reference also lists `reads_tmp.fa`, but it writes
        // `reads_tmp.fa.gz`, so that deletion is a no-op and the file
        // survives. The port writes the uncompressed name and keeps it, to
        // leave the same file present; reconciling the .gz naming is still
        // outstanding (Part 5).
    }
    println!("Done.");
    ExitCode::from(0)
}

extern "C" {
    #[link_name = "dup"]   fn libc_dup(fd: i32) -> i32;
    #[link_name = "dup2"]  fn libc_dup2(a: i32, b: i32) -> i32;
    #[link_name = "close"] fn libc_close(fd: i32) -> i32;
    #[link_name = "fork"]  fn libc_fork() -> i32;
    #[link_name = "_exit"] fn libc__exit(code: i32) -> !;
    #[link_name = "waitpid"] fn libc_waitpid(pid: i32, status: *mut i32, options: i32) -> i32;
}

/// Align every read, on `nr_cores` threads over a SHARED index.
///
/// The reference starts one *process* per core and each re-imports the index
/// from disk (Finding 4). Threads over one `&Index` is the whole reason the
/// port can be faster and smaller at once, so the index is never copied.
///
/// Output order is the input read order, not completion order -- `reads.sam`
/// is a byte-identity target, so a thread count that changed the record order
/// would be a bug, not a detail. Chunks are claimed in order and a chunk that
/// finishes early waits its turn in `pending`.
fn align_all(
    ix: &index::Index,
    recs: &[reads::Record],
    by_name: &std::collections::HashMap<&str, &reads::SeedRecord>,
    p: &driver::Params,
    nr_cores: i64,
    w: &mut (impl Write + Send),
) {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Condvar, Mutex};

    // Small enough that the tail does not idle most threads on a short run,
    // large enough that the handoff is not the cost. Reads vary by an order
    // of magnitude in alignment time, so smaller chunks also even out the
    // load better than they cost.
    const CHUNK: usize = 32;

    let nthreads = nr_cores.max(1) as usize;
    let nchunks = recs.len().div_ceil(CHUNK);
    if nchunks == 0 {
        return;
    }

    let claim = AtomicUsize::new(0);
    // `next` is the chunk the writer wants; `done` holds chunks that finished
    // out of turn. Bounded, or one slow chunk at the front would let every
    // other chunk's output accumulate in memory.
    struct Queue {
        next: usize,
        done: std::collections::BTreeMap<usize, String>,
    }
    let q = Mutex::new(Queue { next: 0, done: std::collections::BTreeMap::new() });
    let room = Condvar::new();
    let backlog = 4 * nthreads;
    let w = Mutex::new(w);

    std::thread::scope(|scope| {
        for _ in 0..nthreads {
            scope.spawn(|| {
                let empty: Vec<String> = Vec::new();
                loop {
                    let i = claim.fetch_add(1, Ordering::Relaxed);
                    if i >= nchunks {
                        break;
                    }
                    let lo = i * CHUNK;
                    let hi = (lo + CHUNK).min(recs.len());
                    let mut buf = String::new();
                    for r in &recs[lo..hi] {
                        let (hits, hits_rc) = match by_name.get(r.name.as_str()) {
                            Some(s) => (&s.hits, &s.hits_rc),
                            None => (&empty, &empty),
                        };
                        for line in driver::align_read(
                            ix, &r.name, &r.seq, r.qual.as_deref(), hits, hits_rc, p) {
                            buf.push_str(&line);
                        }
                    }

                    let mut g = q.lock().unwrap();
                    // Never make the chunk the writer is waiting for wait:
                    // that one is what unblocks everyone else.
                    while g.next != i && g.done.len() >= backlog {
                        g = room.wait(g).unwrap();
                    }
                    g.done.insert(i, buf);
                    let mut out = w.lock().unwrap();
                    while let Some(s) = { let n = g.next; g.done.remove(&n) } {
                        let _ = out.write_all(s.as_bytes());
                        g.next += 1;
                    }
                    drop(out);
                    drop(g);
                    room.notify_all();
                }
            });
        }
    });
}

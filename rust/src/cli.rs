//! An argparse-compatible argument parser for uLTRA.
//!
//! Hand-written, and deliberately not `clap`. The contract is not "accept the
//! same flags" but "produce the same bytes on stdout and stderr and the same
//! exit code", and no general-purpose parser reproduces argparse's wording,
//! wrapping and ordering. The reference's exact messages are in
//! `bench/golden/cli/` and reproduced here.
//!
//! Behaviours pinned by the goldens, each one measured rather than assumed:
//!
//!   * no arguments  -> the top-level help on STDOUT, exit **0**
//!     (uLTRA's own `len(sys.argv)==1` branch, not argparse's)
//!   * `--help`/`-h` -> byte-identical to the no-argument output
//!   * argparse failures -> usage + one error line on STDERR, exit **2**,
//!     and the output folder is NOT created (Finding 17)
//!   * an unrecognised flag reports the **top-level** usage, not the
//!     subcommand's, even when a subcommand was given (Finding: measured in
//!     the `unknown-flag` golden)
//!   * `--alignment_threshold` is `type=int` with a float default, so `0.3`
//!     is rejected (Finding 15)
//!   * the `--ont`/`--isoseq` presets overwrite an explicit `--s` AFTER
//!     parsing (Finding 12)

use crate::text;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sub {
    Pipeline,
    Index,
    Align,
}

impl Sub {
    pub fn name(self) -> &'static str {
        match self {
            Sub::Pipeline => "pipeline",
            Sub::Index => "index",
            Sub::Align => "align",
        }
    }
    fn usage(self) -> &'static str {
        match self {
            Sub::Pipeline => text::USAGE_PIPELINE,
            Sub::Index => text::USAGE_INDEX,
            Sub::Align => text::USAGE_ALIGN,
        }
    }
    fn help(self) -> &'static str {
        match self {
            Sub::Pipeline => text::HELP_PIPELINE,
            Sub::Index => text::HELP_INDEX,
            Sub::Align => text::HELP_ALIGN,
        }
    }
}

/// What `parse` decided. `main` turns this into output and an exit code.
pub enum Outcome {
    /// Print to stdout, exit 0.
    Stdout(String),
    /// Print usage+error to stderr, exit 2. argparse's failure shape.
    UsageError { usage: &'static str, msg: String },
    /// Arguments are valid; run this.
    Run(Box<Args>),
}

#[derive(Debug, Clone)]
pub struct Args {
    pub sub: Sub,
    // positionals
    pub reference: String,
    pub gtf: String,
    pub reads: String,
    pub outfolder: String,
    // shared
    pub nr_cores: i64,
    pub index: String,
    pub prefix: String,
    pub max_intron: i64,
    pub reduce_read_polya: i64,
    pub alignment_threshold: f64,
    pub non_covered_cutoff: i64,
    pub dropoff: f64,
    pub max_loc: f64,
    pub ignore_rc: bool,
    pub min_acc: f64,
    pub disable_mm2: bool,
    pub genomic_frac: f64,
    pub keep_temporary_files: bool,
    pub thinning: i64,
    pub s: i64,
    pub ont: bool,
    pub isoseq: bool,
    // index / pipeline only
    pub min_segm: i64,
    pub flank_size: i64,
    pub small_exon_threshold: i64,
    pub disable_infer: bool,
    // derived after parsing, exactly as the reference does it
    pub mm2_ksize: i64,
}

impl Args {
    fn new(sub: Sub) -> Self {
        // Defaults copied from the reference's add_argument calls.
        Args {
            sub,
            reference: String::new(),
            gtf: String::new(),
            reads: String::new(),
            outfolder: String::new(),
            nr_cores: 3,
            index: String::new(),
            prefix: "reads".into(),
            max_intron: 1_200_000,
            reduce_read_polya: 8,
            // NOTE: the reference declares this type=int with default=0.5.
            // argparse does not coerce a non-string default, so the float
            // default survives while user input must be an int. Finding 15.
            alignment_threshold: 0.5,
            non_covered_cutoff: 15,
            dropoff: 0.95,
            max_loc: 5.0,
            ignore_rc: false,
            min_acc: 0.5,
            disable_mm2: false,
            genomic_frac: 0.1,
            keep_temporary_files: false,
            thinning: 0,
            s: 10,
            ont: false,
            isoseq: false,
            min_segm: 25,
            flank_size: 1000,
            small_exon_threshold: 200,
            disable_infer: false,
            mm2_ksize: 15,
        }
    }

    /// The reference applies the presets AFTER parse_args, overwriting any
    /// explicit --s. Finding 12: `--ont --s 12` silently becomes s=9.
    fn apply_presets(&mut self) {
        if matches!(self.sub, Sub::Align | Sub::Pipeline) {
            self.mm2_ksize = 15;
            if self.ont {
                self.min_acc = 0.6;
                self.mm2_ksize = 14;
                self.s = 9;
            }
            if self.isoseq {
                self.min_acc = 0.8;
                self.s = 10;
            }
        }
    }
}

/// Positional names per subcommand, in order. Drives both consumption and the
/// "the following arguments are required" message.
fn positionals(sub: Sub) -> &'static [&'static str] {
    match sub {
        Sub::Pipeline => &["ref", "gtf", "reads", "outfolder"],
        Sub::Index => &["ref", "gtf", "outfolder"],
        Sub::Align => &["ref", "reads", "outfolder"],
    }
}

fn err(usage: &'static str, prog: &str, msg: String) -> Outcome {
    Outcome::UsageError {
        usage,
        msg: format!("{prog}: error: {msg}"),
    }
}

fn parse_int(v: &str) -> Option<i64> {
    // argparse's int() accepts surrounding whitespace and a leading sign, and
    // rejects anything float-shaped. "0.3" is an error; " 7 " is 7.
    let t = v.trim();
    if t.is_empty() {
        return None;
    }
    t.parse::<i64>().ok()
}

fn parse_float(v: &str) -> Option<f64> {
    let t = v.trim();
    if t.is_empty() {
        return None;
    }
    t.parse::<f64>().ok()
}

pub fn parse(argv: &[String]) -> Outcome {
    // uLTRA's own branch, before argparse gets a say: bare invocation prints
    // the top-level help and exits 0.
    if argv.is_empty() {
        return Outcome::Stdout(text::HELP_TOP.to_string());
    }

    // Top-level options are checked before the subcommand.
    match argv[0].as_str() {
        "-h" | "--help" => return Outcome::Stdout(text::HELP_TOP.to_string()),
        "--version" => return Outcome::Stdout(text::VERSION.to_string()),
        _ => {}
    }

    let sub = match argv[0].as_str() {
        "pipeline" => Sub::Pipeline,
        "index" => Sub::Index,
        "align" => Sub::Align,
        other => {
            return err(
                text::USAGE_TOP,
                "uLTRA",
                format!(
                    "argument {{pipeline,index,align}}: invalid choice: '{other}' \
                     (choose from pipeline, index, align)"
                ),
            )
        }
    };

    let prog = format!("uLTRA {}", sub.name());
    let mut a = Args::new(sub);
    let mut pos: Vec<String> = Vec::new();
    let mut unknown: Vec<String> = Vec::new();
    let mut i = 1usize;

    // Which value-taking flags exist for this subcommand, and how to store them.
    // Splitting int from float matters: the error text differs, and Finding 15
    // lives on exactly this distinction.
    while i < argv.len() {
        let tok = argv[i].clone();

        if tok == "-h" || tok == "--help" {
            return Outcome::Stdout(sub.help().to_string());
        }

        // A flag needing a value; `need` fetches it or reports "expected one argument".
        macro_rules! need {
            () => {{
                i += 1;
                match argv.get(i) {
                    Some(v) => v.clone(),
                    None => {
                        return err(
                            sub.usage(),
                            &prog,
                            format!("argument {tok}: expected one argument"),
                        )
                    }
                }
            }};
        }
        macro_rules! set_int {
            ($field:ident) => {{
                let v = need!();
                match parse_int(&v) {
                    Some(n) => a.$field = n,
                    None => {
                        return err(
                            sub.usage(),
                            &prog,
                            format!("argument {tok}: invalid int value: '{v}'"),
                        )
                    }
                }
            }};
        }
        macro_rules! set_float {
            ($field:ident) => {{
                let v = need!();
                match parse_float(&v) {
                    Some(n) => a.$field = n,
                    None => {
                        return err(
                            sub.usage(),
                            &prog,
                            format!("argument {tok}: invalid float value: '{v}'"),
                        )
                    }
                }
            }};
        }

        let is_align_like = matches!(sub, Sub::Align | Sub::Pipeline);
        let is_index_like = matches!(sub, Sub::Index | Sub::Pipeline);

        match tok.as_str() {
            "--t" if is_align_like => set_int!(nr_cores),
            "--index" if is_align_like => a.index = need!(),
            "--prefix" if is_align_like => a.prefix = need!(),
            "--max_intron" if is_align_like => set_int!(max_intron),
            "--reduce_read_ployA" if is_align_like => set_int!(reduce_read_polya),
            // type=int in the reference, despite the float default. Finding 15.
            "--alignment_threshold" if is_align_like => {
                let v = need!();
                match parse_int(&v) {
                    Some(n) => a.alignment_threshold = n as f64,
                    None => {
                        return err(
                            sub.usage(),
                            &prog,
                            format!("argument {tok}: invalid int value: '{v}'"),
                        )
                    }
                }
            }
            "--non_covered_cutoff" if is_align_like => set_int!(non_covered_cutoff),
            "--dropoff" if is_align_like => set_float!(dropoff),
            "--max_loc" if is_align_like => set_float!(max_loc),
            "--ignore_rc" if is_align_like => a.ignore_rc = true,
            "--min_acc" if is_align_like => set_float!(min_acc),
            "--disable_mm2" if is_align_like => a.disable_mm2 = true,
            "--genomic_frac" if is_align_like => set_float!(genomic_frac),
            "--keep_temporary_files" if is_align_like => a.keep_temporary_files = true,

            "--min_segm" if is_index_like => set_int!(min_segm),
            "--flank_size" if is_index_like => set_int!(flank_size),
            "--small_exon_threshold" if is_index_like => set_int!(small_exon_threshold),
            "--disable_infer" if is_index_like => a.disable_infer = true,

            // Present on all three subcommands.
            "--thinning" => set_int!(thinning),
            "--s" => set_int!(s),

            // The mutually exclusive preset group; align and pipeline only.
            "--ont" if is_align_like => {
                if a.isoseq {
                    return err(
                        sub.usage(),
                        &prog,
                        "argument --ont: not allowed with argument --isoseq".into(),
                    );
                }
                a.ont = true;
            }
            "--isoseq" if is_align_like => {
                if a.ont {
                    return err(
                        sub.usage(),
                        &prog,
                        "argument --isoseq: not allowed with argument --ont".into(),
                    );
                }
                a.isoseq = true;
            }

            // A negative number is a value, not a flag: `--thinning -1` works.
            t if t.starts_with('-') && t.len() > 1 && !t[1..].starts_with(|c: char| c.is_ascii_digit()) => {
                unknown.push(tok.clone());
            }
            _ => pos.push(tok.clone()),
        }
        i += 1;
    }

    // argparse reports unrecognised arguments against the TOP-LEVEL usage,
    // even though a subcommand was parsed. Measured; see the unknown-flag
    // golden.
    if !unknown.is_empty() {
        return err(
            text::USAGE_TOP,
            "uLTRA",
            format!("unrecognized arguments: {}", unknown.join(" ")),
        );
    }

    let names = positionals(sub);
    if pos.len() < names.len() {
        let missing = names[pos.len()..].join(", ");
        return err(
            sub.usage(),
            &prog,
            format!("the following arguments are required: {missing}"),
        );
    }
    if pos.len() > names.len() {
        return err(
            text::USAGE_TOP,
            "uLTRA",
            format!("unrecognized arguments: {}", pos[names.len()..].join(" ")),
        );
    }

    let mut it = pos.into_iter();
    match sub {
        Sub::Pipeline => {
            a.reference = it.next().unwrap();
            a.gtf = it.next().unwrap();
            a.reads = it.next().unwrap();
            a.outfolder = it.next().unwrap();
        }
        Sub::Index => {
            a.reference = it.next().unwrap();
            a.gtf = it.next().unwrap();
            a.outfolder = it.next().unwrap();
        }
        Sub::Align => {
            a.reference = it.next().unwrap();
            a.reads = it.next().unwrap();
            a.outfolder = it.next().unwrap();
        }
    }

    a.apply_presets();
    Outcome::Run(Box::new(a))
}

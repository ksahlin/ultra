//! Replays every case in `bench/cli_cases.tsv` against the recorded goldens.
//!
//! This duplicates what `bench/equivalence.sh cli verify` does, on purpose:
//! `cargo test` is the gate you run every thirty seconds while editing
//! `cli.rs`, and it should not need a conda environment, a corpus or a
//! reference interpreter. The bash harness remains the authority — it is what
//! compares against the *reference*; this compares against what the reference
//! was recorded as doing.

use std::path::{Path, PathBuf};
use std::process::Command;

fn repo_root() -> PathBuf {
    // CARGO_MANIFEST_DIR is <repo>/rust
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().to_path_buf()
}

fn bin() -> PathBuf {
    // Built by cargo before integration tests run.
    let mut p = Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
    p.push(if cfg!(debug_assertions) { "debug" } else { "release" });
    p.push("uLTRA");
    p
}

/// Same scrubbing as bench/equivalence.sh. Keep the two in step.
fn scrub(s: &str, out: &str, root: &str, home: &str) -> String {
    let mut t = s.replace(out, "{OUT}").replace(root, "{ROOT}");
    if !home.is_empty() {
        t = t.replace(home, "{HOME}");
    }
    // "line 123," -> "line N,"
    let mut res = String::with_capacity(t.len());
    let bytes: Vec<char> = t.chars().collect();
    let mut i = 0;
    while i < bytes.len() {
        if t[i..].starts_with("line ") {
            let rest = &t[i + 5..];
            let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
            if !digits.is_empty() && rest[digits.len()..].starts_with(',') {
                res.push_str("line N,");
                i += 5 + digits.len() + 1;
                continue;
            }
        }
        res.push(bytes[i]);
        i += 1;
    }
    res
}

struct Case {
    name: String,
    status: String,
    args: String,
}

fn cases() -> Vec<Case> {
    let tsv = repo_root().join("bench/cli_cases.tsv");
    let body = std::fs::read_to_string(&tsv)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", tsv.display()));
    body.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut f = l.split('\t');
            let name = f.next()?.to_string();
            let status = f.next().unwrap_or("").to_string();
            let args = f.next().unwrap_or("").to_string();
            if name == "name" {
                return None;
            }
            Some(Case { name, status, args })
        })
        .collect()
}

#[test]
fn cli_matches_recorded_goldens() {
    let root = repo_root();
    let root_s = root.to_string_lossy().to_string();
    let home = std::env::var("HOME").unwrap_or_default();
    let binary = bin();
    assert!(
        binary.exists(),
        "binary not built at {} -- run `cargo build` first",
        binary.display()
    );

    let mut failures: Vec<String> = Vec::new();
    let mut checked = 0usize;
    let mut pending: Vec<String> = Vec::new();

    for c in cases() {
        // Every case must be classified; an unrecognised status is a failure,
        // not a skip. Mirrors `bench/equivalence.sh cli_audit`.
        if c.status != "contract" && !c.status.starts_with("pending:") {
            failures.push(format!("{}: unclassified status {:?}", c.name, c.status));
            continue;
        }
        if c.status.starts_with("pending:") {
            pending.push(format!("{} ({})", c.name, c.status));
            continue;
        }
        let gdir = root.join("bench/golden/cli").join(&c.name);
        if !gdir.exists() {
            failures.push(format!("{}: no golden recorded", c.name));
            continue;
        }

        let outdir = std::env::temp_dir().join(format!("ultra-clitest-{}-{}", std::process::id(), c.name));
        let _ = std::fs::remove_dir_all(&outdir);
        let out_s = outdir.to_string_lossy().to_string();

        let expanded = c
            .args
            .replace("{OUT}", &out_s)
            .replace("{REF}", &root.join("test/SIRV_genes.fasta").to_string_lossy())
            .replace("{GTF}", &root.join("test/SIRV_genes_C_170612a.gtf").to_string_lossy())
            .replace("{READS}", &root.join("test/reads.fa").to_string_lossy());
        let argv: Vec<&str> = expanded.split_whitespace().collect();

        let o = Command::new(&binary)
            .args(&argv)
            .current_dir(&root)
            .output()
            .unwrap_or_else(|e| panic!("could not run {}: {e}", binary.display()));

        let got_exit = o.status.code().unwrap_or(-1).to_string();
        let got_out = scrub(&String::from_utf8_lossy(&o.stdout), &out_s, &root_s, &home);
        let got_err = scrub(&String::from_utf8_lossy(&o.stderr), &out_s, &root_s, &home);
        let got_folder = if outdir.is_dir() {
            let n = std::fs::read_dir(&outdir).map(|d| d.count()).unwrap_or(0);
            format!("created=yes entries={n}\n")
        } else {
            "created=no\n".to_string()
        };
        let _ = std::fs::remove_dir_all(&outdir);

        let want = |f: &str| std::fs::read_to_string(gdir.join(f)).unwrap_or_default();
        for (label, got, want) in [
            ("exit", format!("{got_exit}\n"), want("exit")),
            ("stdout", got_out, want("stdout")),
            ("stderr", got_err, want("stderr")),
            ("outfolder", got_folder, want("outfolder")),
        ] {
            if got != want {
                failures.push(format!(
                    "{} [{}]\n     want: {:?}\n     got:  {:?}",
                    c.name,
                    label,
                    want.chars().take(200).collect::<String>(),
                    got.chars().take(200).collect::<String>()
                ));
            }
        }
        checked += 1;
    }

    assert!(checked > 0, "no CLI cases found -- is bench/cli_cases.tsv present?");
    // Pending cases are printed, never hidden. `cargo test -- --nocapture`
    // shows the list; it must shrink to empty as the stages land.
    if !pending.is_empty() {
        println!("{} pending case(s) not yet expected to match:", pending.len());
        for p in &pending {
            println!("    {p}");
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} contract CLI cases differ from the goldens:\n  {}",
        failures.len(),
        checked,
        failures.join("\n  ")
    );
    println!("{checked} contract CLI cases matched the goldens");
}

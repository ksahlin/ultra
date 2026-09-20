//! Fixed output strings.
//!
//! These are **extracted from the recorded goldens** by
//! `bench/extract_cli_text.py`, never typed by hand. argparse's help output is
//! wrapped, indented and ordered in ways that are tedious to reproduce and easy
//! to get subtly wrong; and a 6 KB help block retyped from a screenshot is a
//! byte-identity failure waiting to happen.
//!
//! To refresh after the reference's help changes:
//!     bench/equivalence.sh cli record
//!     python3 bench/extract_cli_text.py

pub const HELP_TOP: &str = include_str!("text/help_top.txt");
pub const HELP_PIPELINE: &str = include_str!("text/help_pipeline.txt");
pub const HELP_INDEX: &str = include_str!("text/help_index.txt");
pub const HELP_ALIGN: &str = include_str!("text/help_align.txt");

pub const USAGE_TOP: &str = include_str!("text/usage_top.txt");
pub const USAGE_PIPELINE: &str = include_str!("text/usage_pipeline.txt");
pub const USAGE_INDEX: &str = include_str!("text/usage_index.txt");
pub const USAGE_ALIGN: &str = include_str!("text/usage_align.txt");

pub const VERSION: &str = include_str!("text/version.txt");

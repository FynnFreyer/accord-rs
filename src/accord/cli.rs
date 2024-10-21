//! This module is responsible for parsing CLI arguments.

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// Path to a reference sequence.
    pub ref_path: String,

    /// Path to a sorted SAM/BAM file containing reads that were aligned against the reference.
    pub aln_path: String,

    /// Optionally, an out path. Defaults to stdout.
    #[arg(short, long, default_value_t = String::from("-"))]
    pub out_path: String,

    /// Alignment quality settings for consensus generation.
    #[arg(short, long, default_value_t = String::new())]
    pub aln_reqs: String,
}

impl Args {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}
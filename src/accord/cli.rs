//! This module is responsible for parsing CLI arguments.

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    pub ref_path: String,

    pub aln_path: String,

    #[arg(short, long, default_value_t = String::from("-"))]
    pub out_path: String,

    #[arg(short, long, default_value_t = String::new())]
    pub aln_reqs: String,
}

impl Args {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}
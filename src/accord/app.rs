//! This module contains the `App` struct, which serves as entry point for the `accord` binary.

use super::utils::{get_fasta_seq, write_file};
use super::calculator::Calculator;
use super::data::settings::AlnQualityReqs;
use super::cli::Args;

pub struct App;

impl App {
    pub fn main() {
        const DEFAULT_REQS: AlnQualityReqs = AlnQualityReqs {
            min_mapq: 10,
            mandatory_flags: 0,
            prohibited_flags: 0,
            indel_cutoff: 0.2,
            save_ends: 24,
            min_observations: 50,
        };

        let args = Args::parse_args();
        let ref_seq = get_fasta_seq(&args.ref_path);
        let aln_path = args.aln_path;

        let mut calculator = Calculator::new(ref_seq, aln_path, DEFAULT_REQS);
        let consensus = calculator.compute_consensus();
        let fasta = consensus.to_fasta();

        let stats = calculator.compute_aln_stats();

        if args.out_path != "-" {
            write_file(&fasta, args.out_path.as_str());
        } else {
            print!("{fasta}");
        }

        println!();
        println!("{stats:?}");
    }
}
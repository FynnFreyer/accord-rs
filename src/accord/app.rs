//! This module contains the `App` struct, which serves as entry point for the `accord` binary.

use super::calculator::Calculator;
use super::cli::Args;
use super::settings::AlnQualityReqs;
use super::utils::write_file;
use crate::accord::data::seq::Seq;

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
        let ref_seqs = Seq::from_file(&args.ref_path);
        let aln_path = args.aln_path;

        let calculator = Calculator::new(DEFAULT_REQS);
        let consensuses = calculator.calculate(ref_seqs, aln_path);

        let mut fastas = Vec::new();
        let mut aln_stats = Vec::new();
        for consensus in consensuses {
            let fasta = consensus.get_consensus_seq().to_fasta();
            fastas.push(fasta);
            let stats = consensus.get_aln_stats().clone();
            aln_stats.push(stats);
        }

        for i in 0..fastas.len() {
            let fasta = fastas.get(i).unwrap();
            let stats = aln_stats.get(i).unwrap();

            if args.out_path != "-" {
                write_file(&fasta, args.out_path.as_str());
            } else {
                println!("{fasta}");
            }

            println!();
            println!("{stats:?}");
        }
    }
}

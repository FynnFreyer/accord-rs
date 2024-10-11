pub mod calculator;
pub mod data;

use clap::Parser;

use crate::utils::{get_fasta_seq, write_file};
use calculator::Calculator;
use data::settings::AlnQualityReqs;


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    ref_path: String,

    aln_path: String,

    #[arg(short, long, default_value_t = String::from("-"))]
    out_path: String,

    #[arg(short, long, default_value_t = String::new())]
    aln_reqs: String,
}


const DEFAULT_REQS: AlnQualityReqs = AlnQualityReqs {
    min_mapq: 10,
    mandatory_flags: 0,
    prohibited_flags: 0,
    indel_cutoff: 0.2,
    save_ends: 24,
    min_observations: 50,
};

pub struct App;

impl App {
    pub fn main() {
        let args = Args::parse();
        let ref_seq = get_fasta_seq(&args.ref_path);
        let aln_path = args.aln_path;

        let mut calculator = Calculator::new(ref_seq, aln_path, DEFAULT_REQS);
        let consensus = calculator.compute_consensus();
        let fasta = consensus.to_fasta();

        let (stats_total, stats_considered) = calculator.compute_aln_stats();

        if args.out_path != "-" {
            write_file(&fasta, args.out_path.as_str());
        } else {
            print!("{fasta}");
        }

        println!();
        println!("{stats_considered:?}");
    }
}

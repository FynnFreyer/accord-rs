//! This module contains settings for consensus generation.

use clap::Args;
use pyo3::{pyclass, pymethods};
use rust_htslib::bam::Record;

/// Requirements for alignment quality.
#[derive(Debug, Clone, Args)]
#[pyclass]
pub struct AlnQualityReqs {
    /// Minimal mapping quality for a read to be considered.
    #[arg(short = 'q', long, alias = "quality", default_value_t = 0)]
    #[pyo3(get)]
    pub min_mapq: u8,

    /// SAM-flags that must be present for a read to be considered.
    #[arg(long, default_value_t = 0)]
    #[pyo3(get)]
    pub mandatory_flags: u16,

    /// SAM-flags that may not be present for a read to be considered.
    #[arg(long, default_value_t = 1540)]
    #[pyo3(get)]
    pub prohibited_flags: u16,

    /// Percentage of coverage that displays an indel for it to be added to the consensus.
    /// E.g. `0.2` means InDels have to appear in 20% of reads that cover that region to be considered.
    #[arg(short = 'c', long)]
    #[pyo3(get)]
    pub indel_cutoff: f64,

    /// Has no purpose at this point.
    ///
    /// Probably useless.
    /// Supposed to avoid artificial lengthening of fragments over multiple iterations.
    /// Included for parity with legacy project.
    #[arg(short, long, default_value_t = 0)]
    #[pyo3(get)]
    pub save_ends: usize,

    /// Minimum coverage needed for considering a position in the consensus calculation.
    #[arg(long, default_value_t = 50)]
    #[pyo3(get)]
    pub min_observations: usize,
}

impl AlnQualityReqs {
    pub fn is_suitable(&self, record: &Record) -> bool {
        //! Calculate whether a given SAM/BAM record is suitable for inclusion while counting bases.
        //!
        //! A read is suitable iff it has the minimum quality, does not have prohibited flags and
        //! has all mandatory flags.

        let qual_ok = record.mapq() >= self.min_mapq;

        let flags = record.flags();
        let mandatory_masked = self.mandatory_flags & flags;
        let has_mandatory = mandatory_masked == self.mandatory_flags;

        let prohibited_masked = self.prohibited_flags & flags;
        let no_prohibited = prohibited_masked == 0;
        let flags_ok = has_mandatory & no_prohibited;

        qual_ok & flags_ok
    }
}

#[pymethods]
impl AlnQualityReqs {
    #[new]
    pub fn new(min_mapq: u8, mandatory_flags: u16, prohibited_flags: u16, indel_cutoff: f64, save_ends: usize, min_observations: usize) -> Self {
        Self { min_mapq, mandatory_flags, prohibited_flags, indel_cutoff, save_ends, min_observations }
    }

    fn __repr__(&self) -> String {
        format!(
            "AlnQualityReqs(min_mapq={}, mandatory_flags={}, prohibited_flags={}, indel_cutoff={}, save_ends={}, min_observations={})",
            self.min_mapq, self.mandatory_flags, self.prohibited_flags, self.indel_cutoff, self.save_ends, self.min_observations
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init_record(mapq: u8, flags: u16) -> Record {
        let mut rec = Record::new();
        rec.set_mapq(mapq);
        rec.set_flags(flags);

        rec
    }

    fn init_with_flags(mandatory: u16, prohibited: u16) -> AlnQualityReqs {
        AlnQualityReqs {
            min_mapq: 10,
            mandatory_flags: mandatory,
            prohibited_flags: prohibited,
            indel_cutoff: 0.2,
            save_ends: 24,
            min_observations: 50,
        }
    }

    #[test]
    fn aln_reqs_mandate_flags() {
        let flags_prohibited: u16 = 0b0101010101010101;
        let flags_neutral: u16 = 0b1010101010101010;

        let rec = init_record(50, flags_neutral);
        let reqs = init_with_flags(0, flags_prohibited);

        assert!(reqs.is_suitable(&rec))
    }

    #[test]
    fn aln_reqs_prohibit_disallows_bad_flags() {
        let flags_prohibited: u16 = 0b0101010101010101;
        let reqs = init_with_flags(0, flags_prohibited);

        let rec_flags: [u16; 8] = [
            0b0000000000000001,
            0b0000000000000100,
            0b0000000000010000,
            0b0000000001000000,
            0b0000000100000000,
            0b0000010000000000,
            0b0001000000000000,
            0b0100000000000000,
        ];

        for flags in rec_flags {
            let rec = init_record(50, flags);
            assert!(!reqs.is_suitable(&rec))
        }
    }

    #[test]
    fn aln_reqs_prohibit_allows_ok_flags() {
        let flags_prohibited: u16 = 0b0101010101010101;
        let reqs = init_with_flags(0, flags_prohibited);

        let rec_flags: [u16; 8] = [
            0b0000000000000010,
            0b0000000000001000,
            0b0000000000100000,
            0b0000000010000000,
            0b0000001000000000,
            0b0000100000000000,
            0b0010000000000000,
            0b1000000000000000,
        ];

        for flags in rec_flags {
            let rec = init_record(50, flags);
            assert!(reqs.is_suitable(&rec))
        }
    }
}

//! This module provides the `Consensus` struct, summarizing the result of a consensus calculation.

use std::collections::HashMap;

use super::indel::InDel;
use super::seq::Seq;
use super::stats::AlnStats;
use pyo3::{pyclass, pymethods};

/// Summarizes the result of calculating a consensus.
#[derive(Debug, Clone)]
#[pyclass]
pub struct Consensus {
    /// The reference sequence.
    #[pyo3(get)]
    ref_seq: Seq,

    /// The generated consensus sequence.
    #[pyo3(get)]
    consensus_seq: Seq,

    /// Statistics of reads, that were considered in the consensus generation.
    #[pyo3(get)]
    aln_stats: AlnStats,

    /// Path to the aligned reads.
    #[pyo3(get)]
    aln_path: String,

    /// Base coverage, relative to the reference sequence.
    #[pyo3(get)]
    coverage: Vec<usize>,

    /// A mapping from base characters to coverage for the respective base, relative to the reference sequence.
    #[pyo3(get)]
    base_counts: ExpandedBaseCounts,

    /// Vector containing the applied indels.
    #[pyo3(get)]
    indels: Vec<InDel>,

    /// Total number of seen reads, including those that were not considered for consensus generation.
    #[pyo3(get)]
    total_reads: usize,
}

#[pymethods]
impl Consensus {
    /// Number of reads considered in the consensus generation.
    #[getter]
    pub fn valid_reads(&self) -> usize {
        self.aln_stats.sample_size()
    }

    /// Number of reads seen, that **were not** considered in the consensus generation.
    #[getter]
    pub fn invalid_reads(&self) -> usize {
        self.mapped_reads() - self.total_reads
    }

    fn __repr__(&self) -> String {
        format!(
            "Consensus(ref_seq={}, aln_path={}, ...)",
            self.ref_seq.get_label(), self.aln_path,
        )
    }
}

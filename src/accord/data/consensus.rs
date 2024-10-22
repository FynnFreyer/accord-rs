//! This module provides the `Consensus` struct, summarizing the result of a consensus calculation.

use std::collections::{HashMap, HashSet};

use pyo3::{pyclass, pymethods};

use super::super::types::{BaseCounts, InDelCounts, Coverage, ExpandedBaseCounts};
use super::seq::Seq;
use super::stats::{AlnData, AlnStats};

/// Intermediary results of analysing aligned reads.
#[derive(Debug, Clone)]
#[pyclass]
pub struct AnalysisResult {
    /// Vector containing valid coverage of the reference genome per base position.
    /// Valid means coverage through aligned reads that suffice the quality criteria.
    #[pyo3(get)]
    pub coverage: Vec<usize>,

    /// Vector with base counts relative to position in reference genome.
    pub base_counts: BaseCounts,

    /// Map with indel counts.
    pub indel_counts: InDelCounts,

    /// Vector containing data for alignments that were considered in consensus generation.
    #[pyo3(get)]
    pub valid_alns: Vec<AlnData>,

    /// Set of IDs of seen reads, regardless of quality.
    /// Used for calculating total number of seen reads.
    #[pyo3(get)]
    pub reads_seen: HashSet<i32>,
}


impl AnalysisResult {
    pub fn new(coverage: Coverage, base_counts: BaseCounts, indel_counts: InDelCounts, valid_alns: Vec<AlnData>, reads_seen: HashSet<i32>) -> Self {
        Self { coverage, base_counts, indel_counts, valid_alns, reads_seen }
    }
}

/// Summarizes the result of calculating a consensus.
#[derive(Debug, Clone)]
#[pyclass]
pub struct Consensus {
    /// The reference sequence.
    #[pyo3(get)]
    ref_seq: Seq,

    /// Path to the aligned reads.
    #[pyo3(get)]
    aln_path: String,

    /// The generated consensus sequence.
    #[pyo3(get)]
    consensus_seq: Seq,

    /// Statistics of reads, that were considered in the consensus generation.
    #[pyo3(get)]
    aln_stats: AlnStats,

    /// Base coverage, relative to the reference sequence.
    #[pyo3(get)]
    coverage: Vec<usize>,

    /// A mapping from base characters to coverage for the respective base, relative to the reference sequence.
    #[pyo3(get)]
    base_counts: ExpandedBaseCounts,

    // /// Vector containing the applied indels.
    // #[pyo3(get)]
    // indels: Vec<InDel>,

    /// Total number of seen reads, including those that were not considered for consensus generation.
    #[pyo3(get)]
    total_reads: usize,
}

#[pymethods]
impl Consensus {
    #[new]
    pub fn new(ref_seq: Seq,
               aln_path: String,
               consensus_seq: Seq,
               aln_stats: AlnStats,
               // indels: Vec<InDel>,
               analysis_result: AnalysisResult) -> Self {
        let coverage = analysis_result.coverage;
        let base_counts = Self::expand_base_counts(&analysis_result.base_counts);
        let total_reads = analysis_result.reads_seen.len();

        Self {
            ref_seq,
            aln_path,
            consensus_seq,
            aln_stats,
            coverage,
            base_counts,
            // indels,
            total_reads,
        }
    }

    /// Number of reads considered in the consensus generation.
    #[getter]
    pub fn valid_reads(&self) -> usize {
        self.aln_stats.sample_size()
    }

    /// Number of reads seen, that **were not** considered in the consensus generation.
    #[getter]
    pub fn invalid_reads(&self) -> usize {
        self.valid_reads() - self.total_reads
    }

    fn __repr__(&self) -> String {
        format!(
            "Consensus(ref_seq={}, aln_path={}, ...)",
            self.ref_seq.get_label(), self.aln_path,
        )
    }
}

impl Consensus {
    fn expand_base_counts(base_counts: &BaseCounts) -> ExpandedBaseCounts {
        let mut counts = HashMap::new();

        for (ref_pos, counter) in base_counts.iter().enumerate() {
            for (base_byte, count) in counter.iter() {
                let base = *base_byte as char;
                let base_coverage = counts.entry(base)
                    .or_insert_with(|| vec![0; base_counts.len()]);

                base_coverage[ref_pos] = *count;
            }
        }

        counts
    }

    pub fn get_ref_seq(&self) -> &Seq { &self.ref_seq }
    pub fn get_aln_path(&self) -> &String { &self.aln_path }
    pub fn get_consensus_seq(&self) -> &Seq { &self.consensus_seq }
    pub fn get_aln_stats(&self) -> &AlnStats{ &self.aln_stats }
    pub fn get_coverage(&self) -> &Coverage { &self.coverage }
    pub fn get_base_counts(&self) -> &ExpandedBaseCounts { &self.base_counts }
    pub fn get_total_reads(&self) -> usize { self.total_reads }
}

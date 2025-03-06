//! This module implements the consensus calculation.

use std::cmp::Ordering;
use std::collections::{HashSet, VecDeque};
use std::iter::Iterator;

use bam::pileup::Indel;
use bam::{Read, Reader};
use counter::Counter;
use itertools::Itertools;
use log::debug;
use pyo3::{pyclass, pymethods};
use rust_htslib::bam;
use rust_htslib::bam::pileup::Alignment;

use super::data;
use super::settings::AlnQualityReqs;
use super::types::{BaseCounts, Coverage, InDelCounts};
use data::consensus::{AnalysisResult, Consensus};
use data::indel::{Deletion, InDel, Insertion};
use data::seq::Seq;
use data::stats::{AlnData, AlnStats};

/// A consensus calculator.
#[derive(Debug)]
#[pyclass]
pub struct Calculator {
    /// Settings for alignment quality.
    /// These determine which reads are considered in the consensus calculation.
    #[pyo3(get)]
    aln_quality_reqs: AlnQualityReqs,
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(aln_quality_reqs: AlnQualityReqs) -> Self {
        Self { aln_quality_reqs }
    }

    #[pyo3(name = "calculate")]
    pub fn py_calculate(&self, ref_path: String, aln_path: String) -> Consensus {
        //! Calculate a consensus for the passed reference and aligned reads.
        //!
        //! - `ref_path: String`: Path to the reference against which the reads were aligned.
        //! - `aln_path: String`: Path to a sorted BAM-file with aligned reads.
        //!
        //! Returns a `Consensus` struct.
        let mut seqs = Seq::from_file(ref_path.clone());

        let ref_seq = match seqs.pop() {
            None => panic!("No sequences in file: {ref_path}"),
            Some(seq) => seq,
        };

        self.calculate(ref_seq, aln_path)
    }
}

impl Calculator {
    pub fn calculate(&self, ref_seq: Seq, aln_path: String) -> Consensus {
        //! Calculate a consensus for the passed reference and aligned reads.
        //!
        //! - `ref_seq: Seq`: The reference against which the reads were aligned.
        //! - `aln_path: String`: Path to a sorted BAM-file with aligned reads.
        //!
        //! Returns a `Consensus` struct.

        let results = self.analyse_alignments(&ref_seq, &aln_path);

        // calculations
        let consensus_seq = self.compute_consensus(&ref_seq, &results);
        let aln_stats = self.compute_aln_stats(&results);

        Consensus::new(ref_seq, aln_path, consensus_seq, aln_stats, results)
    }

    /// Compute the consensus sequence for the seen reads that satisfied the quality criteria.
    fn compute_consensus(&self, ref_seq: &Seq, analysis_result: &AnalysisResult) -> Seq {
        let mut label = String::from(ref_seq.get_label().trim());
        label.push_str(".consensus");

        let base_calling_consensus = self.use_majority_bases(ref_seq, &analysis_result.base_counts);
        let indel_consensus = self.apply_indels(&ref_seq, base_calling_consensus, &analysis_result);

        Seq::new(label, indel_consensus)
    }

    /// Compute alignment statistics for reads considered in the consensus calculation.
    fn compute_aln_stats(&self, analysis_result: &AnalysisResult) -> AlnStats {
        let quantile_factors = vec![
            0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.0,
        ];
        let stats = AlnStats::from_data(&analysis_result.valid_alns, &quantile_factors);

        stats
    }

    fn use_majority_bases(&self, ref_seq: &Seq, base_counts: &BaseCounts) -> Vec<u8> {
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(ref_seq.len());

        for (ref_pos, base_counter) in base_counts.iter().enumerate() {
            // get original base at `ref_pos`
            let reference_base = ref_seq[ref_pos];

            // determine consensus by simple majority
            let consensus_base;
            if base_counter.is_empty() {
                // no coverage -> use reference base
                consensus_base = reference_base;
            } else {
                // has coverage
                let (most_common, observations) = *base_counter.most_common().first().unwrap();
                let sufficient_observations =
                    observations >= self.aln_quality_reqs.min_observations;

                consensus_base = if sufficient_observations {
                    most_common
                } else {
                    reference_base
                };
            }

            consensus_seq.push(consensus_base);
        }

        consensus_seq
    }

    fn analyse_alignments(&self, ref_seq: &Seq, aln_path: &String) -> AnalysisResult {
        let mut alns = match Reader::from_path(aln_path.as_str()) {
            Ok(reader) => reader,
            Err(_e) => panic!("Unable to open BAM file: {_e}"),
        };

        let mut coverage = vec![0; ref_seq.len()];
        let mut base_counts = vec![Counter::new(); ref_seq.len()];
        let mut reads_seen = HashSet::new();
        let mut indel_counts = Counter::new();
        let mut valid_alns = Vec::new();

        for p in alns.pileup() {
            // `pileup` holds references to all reads that were aligned to a specific position
            let pileup = match p {
                Ok(p) => p,
                Err(_e) => panic!("Unable to generate pileup: {_e}"),
            };

            let ref_pos = pileup.pos() as usize;
            debug!("Analysing pileup in position {ref_pos}.");

            for alignment in pileup.alignments() {
                // the SAM record of the aligned read
                let record = alignment.record();

                // register read as seen
                let read_id = record.tid();
                reads_seen.insert(read_id);

                // discard read alignments with insufficient quality, flags, etc.
                if !self.aln_quality_reqs.is_suitable(&record) {
                    let read_name = String::from_utf8_lossy(record.qname());
                    debug!("Skipped low quality alignment for read: {}", read_name);
                    continue;
                }

                // register valid alignment
                let aln_data = AlnData::from_record(&record);
                valid_alns.push(aln_data);

                self.register_position(&alignment, &ref_pos, &mut base_counts, &mut coverage);
                self.register_indels(&alignment, &ref_pos, &mut indel_counts);
            }
        }

        AnalysisResult::new(coverage, base_counts, indel_counts, valid_alns, reads_seen)
    }

    fn register_position(
        &self,
        alignment: &Alignment,
        ref_pos: &usize,
        base_counts: &mut BaseCounts,
        coverage: &mut Coverage,
    ) {
        //! Register alignment for position relative to the reference sequence, and update coverage and base counts.

        let record = alignment.record();
        let seq = record.seq();

        let has_read_pos = !alignment.is_refskip() && !alignment.is_del();
        if has_read_pos {
            // find position in read
            let read_pos = alignment.qpos().unwrap();

            // register the base of this read in this position
            let bases = &mut base_counts[*ref_pos];
            let base = seq[read_pos];
            bases[&base] += 1;

            // increment coverage
            coverage[*ref_pos] += 1;
        }
    }

    fn register_indels(
        &self,
        alignment: &Alignment,
        ref_pos: &usize,
        indel_counts: &mut InDelCounts,
    ) {
        let record = alignment.record();
        let read_name = String::from_utf8_lossy(record.qname());
        let indel = match alignment.indel() {
            Indel::Ins(len) => {
                let ins = Self::compute_insertion(len, *ref_pos, &alignment);
                let start = ins.get_start();
                debug!("{read_name} contains insertion of length {len} after {start}.");
                ins
            }
            Indel::Del(len) => {
                let del = Self::compute_deletion(len, *ref_pos);
                let (start, stop) = (del.get_start(), del.get_stop());
                debug!("{read_name} contains deletion between positions {start} and {stop}.");
                del
            }
            Indel::None => return,
        };
        indel_counts.update([indel]);
    }

    fn compute_insertion(len: u32, ref_pos: usize, alignment: &Alignment) -> InDel {
        // let read_name = String::from_utf8_lossy(record.qname());  // used for logging
        // println!("{}: Insertion of length {} between this and next position.", read_name, len);

        let len = len as usize;
        let record = &alignment.record();
        let seq = record.seq();

        let ins_start = alignment.qpos().unwrap() + 1;
        let mut ins_seq = Vec::with_capacity(len);
        for i in ins_start..ins_start + len {
            let base = seq[i];
            ins_seq.push(base);
        }

        let ins = Insertion::new(ref_pos, ins_seq);
        InDel::Ins(ins)
    }

    fn compute_deletion(len: u32, ref_pos: usize) -> InDel {
        let len = len as usize;

        let del_start = ref_pos + 1;
        let del_stop = del_start + len;

        let del = Deletion::new(del_start, del_stop);
        InDel::Del(del)
    }

    fn apply_indels(
        &self,
        ref_seq: &Seq,
        seq_bytes: Vec<u8>,
        analysis_result: &AnalysisResult,
    ) -> Vec<u8> {
        let applicable_indels =
            self.get_applicable_indels(&analysis_result.indel_counts, &analysis_result.coverage);
        let ref_len = ref_seq.len();

        // we prepend string slices to this vector from which we later construct the consensus
        let mut vd: VecDeque<&[u8]> = VecDeque::new();

        // we get slices from the event stop to the start of the previous event
        // "previous" in the sense of previous iteration, but positionally next
        let mut prev_event_start = ref_len;
        for indel in applicable_indels {
            let event_stop = indel.get_stop();

            // skip if this indel interferes with the last applied indel
            let interferes = prev_event_start < event_stop  // events overlap
                || prev_event_start.abs_diff(event_stop) <= 1; // events are adjacent
            let is_first = prev_event_start == ref_len;
            let skip = interferes && !is_first;
            if skip {
                continue;
            }

            // add unaffected sequence part in between events
            let between_range = event_stop + 1..prev_event_start;
            let between = &seq_bytes[between_range];
            vd.push_front(between);
            // add event sequence
            vd.push_front(indel.get_seq());
            // amend positional cutoff for next iteration
            prev_event_start = indel.get_start();
        }

        // push sequence from absolute start to start of first event
        let rest = &seq_bytes[0..prev_event_start];
        vd.push_front(rest);

        // construct indel consensus by copying the slice bytes into the vector
        let mut consensus = Vec::with_capacity(ref_len);
        for slice in vd {
            for byte in slice {
                consensus.push(*byte);
            }
        }

        consensus
    }

    fn get_applicable_indels<'a>(
        &self,
        indel_counts: &'a InDelCounts,
        coverage: &Coverage,
    ) -> VecDeque<&'a InDel> {
        //! Get a vector of indel references, where indels are filtered by whether they're
        //! applicable, and ordered from back to front, for easy insertion.

        // filter indels by whether they have sufficient observations and
        // by whether they make the percentage cutoff for this positions coverage
        let filtered_by_coverage = indel_counts.iter().filter(|(indel, count)| {
            let count = **count;

            let has_min_obs = count > self.aln_quality_reqs.min_observations;

            let indel_cov = &coverage[indel.range()];
            let total_cov = indel_cov.iter().sum::<usize>() as f64;
            let avg_cov = total_cov / indel_cov.len() as f64;

            let required_cov = avg_cov * self.aln_quality_reqs.indel_cutoff;
            let has_required_cov = required_cov <= count as f64;

            has_min_obs && has_required_cov
        });

        // resolve order preferentially, where importance looks like so:
        // position > count > orf breakage > type
        let ordered_by_preference =
            filtered_by_coverage.sorted_by(|(indel_a, count_a), (indel_b, count_b)| {
                let pos_cmp = indel_a.get_start().cmp(&indel_b.get_start());
                if !matches!(pos_cmp, Ordering::Equal) {
                    return pos_cmp;
                }

                let count_cmp = count_a.cmp(count_b);
                if !matches!(count_cmp, Ordering::Equal) {
                    return count_cmp;
                }

                let pref_a = indel_a.preserves_reading_frame() && indel_b.breaks_reading_frame();
                let pref_b = indel_b.preserves_reading_frame() && indel_a.breaks_reading_frame();
                let orf_breakage;
                if pref_a {
                    orf_breakage = Ordering::Greater;
                } else if pref_b {
                    orf_breakage = Ordering::Less;
                } else {
                    orf_breakage = Ordering::Equal;
                };
                if !matches!(orf_breakage, Ordering::Equal) {
                    return orf_breakage;
                }

                // TODO: ask Britta for proper statement as to why
                // we prefer insertions over deletions (because they "add" information as opposed to dels?)
                let type_preference = match indel_a {
                    InDel::Ins(_) => match indel_b {
                        InDel::Ins(_) => Ordering::Equal,
                        InDel::Del(_) => Ordering::Greater,
                    },
                    InDel::Del(_) => match indel_b {
                        InDel::Ins(_) => Ordering::Less,
                        InDel::Del(_) => Ordering::Equal,
                    },
                };
                type_preference
            });

        // reverse order front to back
        let reversed = ordered_by_preference.rev();

        // remove counts (irrelevant after resolving preference)
        let indels = reversed.map(|(indel, _count)| indel);

        indels.collect::<VecDeque<&InDel>>()
    }
}

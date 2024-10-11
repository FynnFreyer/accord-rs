use crate::accord::data;

use std::cmp::Ordering;
use std::collections::VecDeque;
use std::iter::Iterator;

use bam::pileup::Indel;
use bam::{Read, Reader};
use counter::Counter;
use itertools::Itertools;
use pyo3::pyclass;
use rust_htslib::bam;
use rust_htslib::bam::pileup::Alignment;

use data::indel::{Deletion, InDel, Insertion};
use data::seq::Seq;
use data::settings::AlnQualityReqs;
use data::stats::{AlnData, AlnStats};

/// A list of base counts for every position in the reference sequence.
type BaseCounts = Vec<Counter<u8>>;

/// A map in which encountered insertions point to their respective number of occurrences.
type InDelCounts = Counter<InDel, usize>;


// /// A map in which encountered insertions point to their respective number of occurrences.
// type InsertionCounts = Counter<Insertion, usize>;
//
// /// A map in which encountered deletions point to their respective number of occurrences.
// type DeletionCounts = Counter<Deletion, usize>;

#[pyclass]
pub struct Calculator {
    ref_seq: Seq,
    aln_path: String,
    aln_quality_reqs: AlnQualityReqs,
    coverage: Vec<usize>,
    base_counts: BaseCounts,
    indel_counts: InDelCounts,
    aln_data: Vec<AlnData>,
    // insertion_counts: InsertionCounts,
    // deletion_counts: DeletionCounts,
}

impl Calculator {
    pub fn new(ref_seq: Seq, aln_path: String, aln_quality_reqs: AlnQualityReqs) -> Self {
        let coverage = vec![0; ref_seq.len()];
        let base_counts = vec![Counter::new(); ref_seq.len()];
        let indel_counts = Counter::new();
        let aln_data = Vec::new();

        Self { ref_seq, aln_path, aln_quality_reqs, coverage, base_counts, indel_counts, aln_data }
    }

    // This consumes `self`.
    pub fn compute_aln_stats(self) -> AlnStats {
        AlnStats::from_data(self.aln_data)
    }

    pub fn compute_consensus(&mut self) -> Seq {
        self.analyse_alignments();
        let base_calling_consensus = self.use_majority_bases();
        let indel_consensus = self.apply_indels(base_calling_consensus);

        let mut label = String::from(self.ref_seq.label_string());
        label.push_str(".consensus");

        Seq::new(label, indel_consensus)
    }

    fn use_majority_bases(&self) -> Vec<u8> {
        let mut consensus: Vec<u8> = Vec::with_capacity(self.ref_seq.len());

        for (ref_pos, base_counter) in self.base_counts.iter().enumerate() {
            // get original base at `ref_pos`
            let reference_base = self.ref_seq[ref_pos];

            // determine consensus by simple majority
            let consensus_base;
            if base_counter.is_empty() {  // no coverage -> use reference base
                consensus_base = reference_base;
            } else {  // has coverage
                let (most_common, observations) = *base_counter.most_common().first().unwrap();
                let sufficient_observations = observations >= self.aln_quality_reqs.min_observations;

                consensus_base = if sufficient_observations { most_common } else { reference_base };
            }

            consensus.push(consensus_base);
        }

        consensus
    }

    fn analyse_alignments(&mut self) {
        let mut good_reads: usize = 0;
        let mut bad_reads: usize = 0;

        let mut alns = Reader::from_path(self.aln_path.as_str()).unwrap();
        for p in alns.pileup() {
            // `pileup` holds references to all reads that were aligned to a specific position
            let pileup = match p {
                Ok(p) => p,
                Err(e) => panic!("Unable to generate pileup: {e}"),
            };

            let ref_pos = pileup.pos() as usize;

            for alignment in pileup.alignments() {
                // the SAM record of the aligned read
                let record = alignment.record();
                // TODO only consider good data?
                let aln_data = AlnData::from_record(&record);
                self.aln_data.push(aln_data);

                // discard read alignments with insufficient quality, flags, etc.
                if !self.aln_quality_reqs.is_suitable(&record) {
                    bad_reads += 1;

                    // TODO: use proper logging
                    // let read_name = String::from_utf8_lossy(record.qname());  // used for logging
                    // println!("Skipped low quality alignment for read: {}", read_name);
                    continue;
                }

                good_reads += 1;

                self.update_base_counts(&alignment, &ref_pos);
                self.update_indel_counts(&alignment, &ref_pos);
            }

            // println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
        }
    }

    fn update_base_counts(&mut self, alignment: &Alignment, ref_pos: &usize) {
        let record = alignment.record();
        let seq = record.seq();

        let has_read_pos = !alignment.is_refskip() && !alignment.is_del();
        if has_read_pos {
            // find position in read
            let read_pos = alignment.qpos().unwrap();

            // register the base of this read in this position
            let bases = &mut self.base_counts[*ref_pos];
            let base = seq[read_pos];
            bases[&base] += 1;

            // increment coverage
            self.coverage[*ref_pos] += 1;

            // TODO: use proper logging
            // println!("Base {}@{}", base as char, ref_pos);
        }
    }

    fn update_indel_counts(&mut self, alignment: &Alignment, ref_pos: &usize) {
        let indel = match alignment.indel() {
            Indel::Ins(len) => Self::compute_insertion(len, *ref_pos, &alignment),
            Indel::Del(len) => Self::compute_deletion(len, *ref_pos),
            Indel::None => return
        };
        self.indel_counts.update([indel]);
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
        // let read_name = String::from_utf8_lossy(record.qname());  // used for logging
        // println!("{}: Deletion of length {} between this and next position.", read_name, len);

        let len = len as usize;

        let del_start = ref_pos + 1;
        let del_stop = del_start + len;

        let del = Deletion::new(del_start, del_stop);
        InDel::Del(del)
    }

    //         # twe sort by position, so we can find conflicts
    //         sorted_indel_stats = sorted(filtered_stats, key=attrgetter("indel"))
    //
    //         # now we check for conflicts, by grouping the indels,
    //         # so that conflicting events go into the same sub-list
    //         # e.g. we want [[a], [b, c], [d]], if b and c conflict with each other
    //
    //         if not sorted_indel_stats:  # we're initializing the first group, so we have to bail if there's no data
    //             return apply, discard
    //
    //         groups: list[set[InDelStats]] = []
    //         group: set[InDelStats] = {sorted_indel_stats[0]}
    //         for indel_stat, next_indel_stat in pairwise(sorted_indel_stats):
    //             conflict = indel_stat.indel.interferes_with(next_indel_stat.indel)
    //             if conflict:
    //                 group.add(next_indel_stat)
    //             else:
    //                 groups.append(group)
    //                 group = {next_indel_stat}
    //         groups.append(group)  # add last group
    //
    //         for group in groups:
    //             flip = max(group, key=attrgetter("frequency"))
    //             if flip.frequency >= self.quality_requirements.indel_cutoff:
    //                 apply.add(flip)
    //                 flops = group - {flip}
    //             else:
    //                 flops = group
    //             discard.update(flops)
    //
    //         return apply, discard


    fn apply_indels(&self, seq_bytes: Vec<u8>) -> Vec<u8> {
        let applicable_indels = self.get_applicable_indels();
        let ref_len = self.ref_seq.len();

        // we prepend string slices to this vector from which we later construct the consensus
        let mut vd: VecDeque<&[u8]> = VecDeque::new();

        // we get slices from the event stop to the start of the previous event
        // "previous" in the sense of previous iteration, but positionally next
        let mut prev_event_start = ref_len;
        for indel in applicable_indels {
            let event_stop = indel.get_stop();

            // skip if this indel interferes with the last applied indel
            let interferes
                = prev_event_start < event_stop  // events overlap
                || prev_event_start.abs_diff(event_stop) <= 1;  // events are adjacent
            let is_first = prev_event_start == ref_len;
            let skip = interferes && !is_first;
            if skip { continue; }

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

        // TODO: ask Britta for proper statement as to why
        // we prefer insertions over deletions (because they "add" information as opposed to dels?)
        // while !applicable_indels.is_empty() {
        //     let (indel, indel_count) = match applicable_indels.pop_back() {
        //         None => panic!("Popped from empty deque in `apply_indels`"),
        //         Some(indel) => indel,
        //     };
        //
        //     if indel.type_id() == TypeId::of::<Insertion>() {} else {}
        //
        //     let has_next = !applicable_indels.is_empty();
        // }
        //
        // for (insertion, count) in valid_insertions {
        //     println!("{} @ {}, was seen {} times", insertion.seq_string(), insertion.position, count)
        // }

    }

    fn get_applicable_indels(&self) -> VecDeque<&InDel> {
        //! Get a vector of indel references, where indels are filtered by whether they're
        //! applicable, and ordered from back to front, for easy insertion.

        let iter = self.indel_counts.iter();

        // filter indels by whether they have sufficient observations and
        // by whether they make the percentage cutoff for this positions coverage
        let filtered_by_coverage = iter.filter(
            |(indel, count)| {
                let count = **count;

                let has_min_obs = count > self.aln_quality_reqs.min_observations;

                let indel_cov = &self.coverage[indel.range()];
                let total_cov = indel_cov.iter().sum::<usize>() as f64;
                let avg_cov = total_cov / indel_cov.len() as f64;

                let required_cov = avg_cov * self.aln_quality_reqs.indel_cutoff;
                let has_required_cov = required_cov <= count as f64;

                has_min_obs && has_required_cov
            });

        // resolve order preferentially, where importance looks like so:
        // position > count > orf breakage > type
        let ordered_by_preference = filtered_by_coverage.sorted_by(
            |(indel_a, count_a), (indel_b, count_b)| {
                let pos_cmp = indel_a.get_start().cmp(&indel_b.get_start());
                if !matches!(pos_cmp, Ordering::Equal) { return pos_cmp; }

                let count_cmp = count_a.cmp(count_b);
                if !matches!(count_cmp, Ordering::Equal) { return count_cmp; }

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
                if !matches!(orf_breakage, Ordering::Equal) { return orf_breakage; }

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

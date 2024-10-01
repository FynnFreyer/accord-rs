use std::collections::HashMap;
use rust_htslib::bam::Reader;
use crate::accord::counters::{BaseCounter, InDelCounter};
use crate::accord::data::alphabet::Alphabet;
use crate::accord::data::seq_io::Seq;
use crate::accord::data::settings::AlnQualityReqs;

pub struct AlnCounter<'a> {
    /// The reference sequence against which the reads where aligned.
    ref_seq: Seq,

    /// Quality requirements for aligned reads.
    settings: AlnQualityReqs,

    /// A struct for counting bases.
    base_counts: BaseCounter<'a>,

    /// A struct for considering indels.
    indel_counts: InDelCounter<'a>,

    /// The alphabet.
    alphabet: Alphabet,
}

impl<'a> AlnCounter<'a> {
    pub fn new(ref_seq: Seq, settings: AlnQualityReqs, alphabet: Alphabet) -> Self {
        let obj = Self { ref_seq, base_counts, indel_counts, settings, alphabet };

        let ref_seq_reference: &'a Seq = &obj.ref_seq;
        let settings_reference: &'a AlnQualityReqs = &obj.settings;
        let alphabet_reference: &'a Alphabet = &obj.alphabet;

        let base_counts: BaseCounter<'a>
            = BaseCounter::new(ref_seq_reference,
                               settings_reference,
                               alphabet_reference);

        let indel_counts = InDelCounter::new();

        obj
    }

    pub fn count(&mut self, mut hts_reader: Reader) -> () {
        for result in hts_reader.records() {
            let record = result.expect("Invalid record in HTS file");
            // update base and indel counts
            self.base_counts.update(&record);
            self.indel_counts.update(&record);
        }
    }

    pub fn get_base_counts(&self) -> &Vec<Vec<usize>> {
        self.base_counts.get_counts()
    }
}


pub fn foo(ref_seq: &str, sam_path: &str) {
    let mut bam = Reader::from_path(sam_path).unwrap();

    for result in bam.rc_records() {
        let record = result.unwrap();
        let mapq = record.mapq();
        let flags = record.flags();

        println!("{mapq}");
    }

    // let header = bam::Header::from_template(bam.header());
    //
    // // print header records to the terminal, akin to samtool
    // for (key, records) in header.to_hashmap() {
    //     for record in records {
    //         //println!("@{}\tSN:{}\tLN:{}", key, record["SN"], record["LN"]);
    //         println!(key)
    //     }
    // }
}


struct BaseCounter {
    alphabet: Alphabet,
    counts: BaseCounts,
}

impl BaseCounter {
    pub fn new(ref_seq: &Seq, alphabet: Alphabet) -> Self {
        let counts = vec![vec![0; ref_seq.seq_len()]; alphabet.size()];

        Self { alphabet, counts }
    }

    pub fn get_counts(&self) -> &Vec<Vec<usize>> { &self.counts }

    pub fn update(&mut self, mut record: &Record) {
        let offset = record.pos() as usize;
        for (i, byte) in record.seq().as_bytes().iter().enumerate() {
            let symbol = *byte as char;
            let sym_pos = self.alphabet[symbol];
            let seq_pos = i + offset;
            self.counts[sym_pos][seq_pos] += 1;
        }
    }
}

struct InDelCounter {
    insertions: InsertionCounts,
    deletions: DeletionCounts,
}

impl InDelCounter {
    pub fn new() -> Self {
        let insertions = HashMap::new();
        let deletions = HashMap::new();

        Self { insertions, deletions }
    }

    fn register_insertions(&self, record: &Record) -> Vec<Insertion> {
        let insertions = Vec::new();
        let cigar = record.cigar();
        insertions
    }
    fn register_deletions(&self, record: &Record) -> Vec<Deletion> {
        let deletions = Vec::new();
        deletions
    }

    pub fn update(&mut self, record: &Record) {
        let insertions = self.register_insertions(record);
        for insertion in insertions {
            let current = match self.insertions.get(&insertion) {
                None => 0,
                Some(v) => *v
            };
            self.insertions.insert(insertion, current + 1);
        }
        let deletions = self.register_deletions(record);
        for deletion in deletions {
            let current = match self.deletions.get(&deletion) {
                None => 0,
                Some(v) => *v
            };
            self.deletions.insert(deletion, current + 1);
        }
    }
}


// Old Count method on Calculator

// pub fn count(mut self) -> (BaseCounts, InsertionCounts, DeletionCounts) {
//     self.pile();
//
//     for result in self.aln_reader.records() {
//         // get the SAM data record for a specific read
//         let record = match result {
//             Ok(record) => record,
//             Err(e) => panic!("Encountered invalid record in SAM-file: {e}")
//         };
//
//         // skip aligned reads with insufficient quality etc.
//         if !self.aln_quality_reqs.is_suitable(&record) { continue; }
//
//         // update base and indel counts
//         self.base_counter.update(&record);
//         self.indel_counter.update(&record);
//     }
//
//     // TODO calculate base counts
//     let counts = vec![vec![0; self.ref_seq.seq_len()]; 4];
//     (counts, self.indel_counter.insertions, self.indel_counter.deletions)
// }


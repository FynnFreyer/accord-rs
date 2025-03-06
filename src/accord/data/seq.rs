//! This module provides a `Seq` struct for working with sequence data.

use itertools::Itertools;
use pyo3::types::PyType;
use pyo3::{pyclass, pymethods, Bound};
use std::cmp::min;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::fs;
use std::ops::Index;
use std::slice::SliceIndex;

#[derive(Debug, Clone)]
#[pyclass]
pub struct Seq {
    /// The sequence name.
    #[pyo3(get)]
    label: String,

    /// The sequence bytes as a vector name.
    sequence: Vec<u8>, // `Vec<u8>` allows for O(1) indexing, as opposed to `String`.
}

impl Seq {
    pub fn new(label: String, sequence: Vec<u8>) -> Self {
        Self { label, sequence }
    }

    pub fn from_string(label: String, seq_string: String) -> Self {
        let sequence = Vec::from(seq_string.as_bytes());
        Self::new(label, sequence)
    }

    pub fn label_str(&self) -> &str {
        &self.label.as_str()
    }

    /// Get the sequence length.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Parse `Seq`s from a FASTA string.
    pub fn from_fasta(fasta: String) -> Vec<Self> {
        // vector to collect all sequences in this fasta
        let mut seqs = Vec::new();

        // defined outside the loop, so we can reference it between iterations
        let mut label_string = String::new();

        // we group lines by wether they are labels, i.e., start with `>`, or not
        let groups = &fasta.lines().into_iter().chunk_by(|l| l.starts_with(">"));
        for (is_label, lines_group) in groups {
            // collect the actual lines into a vector
            let lines = lines_group.collect::<Vec<&str>>();
            if is_label {
                // this is a label group
                if lines.len() > 1 {
                    // throw if the fasta has more than one label line
                    panic!("More than one label line:\n{:?}", lines)
                }

                // remove the leading `>` and trim start and end
                let clean_label = lines[0][1..].trim();
                // set the label string
                label_string = String::from(clean_label);
            } else {
                // this group contains sequence data
                // clone the label, so we can move it
                let label = label_string.clone();
                // get the sequence data as bytes
                let sequence = Vec::from(lines.join("").as_bytes());

                // construct a Seq and push it into the sequence vector
                seqs.push(Self { label, sequence });
            }
        }

        seqs
    }

    pub fn from_file(file: String) -> Vec<Self> {
        let content = match fs::read_to_string(file) {
            Ok(content) => content,
            Err(e) => panic!("{e}"),
        };
        Self::from_fasta(content)
    }

    /// Get the label of the `Seq`.
    pub fn get_label(&self) -> &String {
        &self.label
    }

    pub fn get_sequence(&self) -> &Vec<u8> {
        &self.sequence
    }
    pub fn get_sequence_as_string(&self) -> String {
        String::from_utf8(self.sequence.clone()).unwrap()
    }
}

#[pymethods]
impl Seq {
    #[new]
    fn py_new(label: String, sequence: String) -> Self {
        Self::from_string(label, sequence)
    }

    #[classmethod]
    #[pyo3(name = "from_fasta")]
    fn py_from_fasta(_cls: &Bound<'_, PyType>, fasta: String) -> Vec<Self> {
        Self::from_fasta(fasta)
    }

    /// Get the sequence data as a string, not a vector.
    #[getter]
    #[pyo3(name = "sequence")]
    fn py_sequence(&self) -> String {
        let seq_bytes = self.sequence.clone();
        String::from_utf8(seq_bytes).unwrap()
    }

    #[classmethod]
    #[pyo3(name = "from_file")]
    pub fn py_from_file(_cls: &Bound<'_, PyType>, file: String) -> Vec<Self> {
        Self::from_file(file)
    }

    /// Convert a `Seq` into a FASTA string.
    pub fn to_fasta(&self) -> String {
        let mut fasta = String::new();

        // write label
        fasta.push('>');
        fasta.push_str(self.label.as_str());
        fasta.push('\n');

        // write sequence
        let width = 80; // lines wrap after 80 chars
        let lines = self.sequence.len().div_ceil(width);
        let len = self.sequence.len();
        for line in 0..lines {
            let start = line * width;
            let full_end = start + width;
            let end = min(full_end, len);

            let line_slice = &self.sequence[start..end];
            let line_content = match std::str::from_utf8(line_slice) {
                Ok(s) => s,
                Err(e) => panic!("Unable to convert Seq into string: {}", e),
            };

            fasta.push_str(line_content);
            fasta.push('\n');
        }

        fasta
    }

    /// Python dunder method for `len`.
    /// Returns sequence length, ignoring the label.
    fn __len__(&self) -> usize {
        self.len()
    }

    /// Python dunder method for `str`.
    /// Returns sequence as string, ignoring the label.
    fn __str__(&self) -> String {
        self.py_sequence()
    }

    fn __repr__(&self) -> String {
        let seq_string = self.py_sequence();
        let mut seq = String::new();
        if seq_string.len() > 20 {
            seq_string.as_str()[..20].clone_into(&mut seq);
            seq.push_str("...");
        } else {
            seq_string.clone_into(&mut seq)
        }

        format!("Seq(label='{}', sequence='{seq}')", self.label)
    }
}

impl Display for Seq {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let seq = self.to_fasta();
        f.write_str(seq.as_str())
    }
}

impl<Idx: SliceIndex<[u8]>> Index<Idx> for Seq {
    type Output = Idx::Output;
    fn index(&self, index: Idx) -> &Self::Output {
        self.sequence.index(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init_seq() -> Seq {
        let label = String::from("Nucleotides");
        let sequence = String::from("GATTACA");
        Seq::from_string(label, sequence)
    }

    #[test]
    fn test_seq() {
        use std::fs;
        let fasta = match fs::read_to_string("assets/K03455.fasta") {
            Ok(content) => content,
            Err(e) => panic!("{e}"),
        };

        let mut seqs = Seq::from_fasta(fasta);

        // correct number of sequences
        assert_eq!(seqs.len(), 1);

        let reference = match seqs.pop() {
            None => panic!("No ref (wat?!)"), // shouldn't happen because we check len before
            Some(seq) => seq,
        };

        // label was read correctly
        assert_eq!(
            reference.label,
            "HXB2_K03455 (Human_immunodeficiency_virus_type_1-complete_genome)"
        );

        // sequence was read correctly
        assert_eq!(reference.sequence.len(), 9719);
        let ref_seq_string = reference.get_sequence_as_string();
        assert!(ref_seq_string.starts_with("TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCC"));
        assert!(ref_seq_string.ends_with("CCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA"));
    }
}

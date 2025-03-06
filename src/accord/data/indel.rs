//! Module that provides structs for working with indels.

use pyo3::{pyclass, pymethods};
use std::hash::Hash;
use std::ops::RangeInclusive;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[pyclass]
pub enum InDel {
    Ins(Insertion),
    Del(Deletion),
}

#[pymethods]
impl InDel {
    // TODO: rephrase this in terms of 5' and 3' ends
    /// Get the starting position of this indel event, independent of read direction.
    /// That means the left side start of the event when considering forward direction.
    pub fn get_start(&self) -> usize {
        match self {
            InDel::Ins(ins) => ins.position,
            InDel::Del(del) => del.start,
        }
    }

    /// Get the ending position of this indel, in forward reading direction.
    pub fn get_stop(&self) -> usize {
        match self {
            InDel::Ins(ins) => ins.position + 1,
            InDel::Del(del) => del.stop,
        }
    }

    fn __len__(&self) -> usize {
        self.len()
    }

    /// Get a byte slice corresponding to this events sequence.
    /// I.e. what should be inserted between event start and stop.
    pub fn get_seq(&self) -> &[u8] {
        match self {
            InDel::Ins(ins) => ins.sequence.as_slice(),
            InDel::Del(_) => &[],
        }
    }

    /// Whether this indel preserves the reading frame by only shifting it by a multiple of three.
    pub fn preserves_reading_frame(&self) -> bool {
        self.len() % 3 == 0
    }

    /// Whether this indel breaks the reading frame,  by shifting it by a non-multiple of three.
    pub fn breaks_reading_frame(&self) -> bool {
        !self.preserves_reading_frame()
    }
}

impl InDel {
    /// Base positions spanning the event site as an inclusive range `start..=stop`.
    /// Start and stop are independent of read direction, and you may assume order `start <= stop`.
    pub fn range(&self) -> RangeInclusive<usize> {
        self.get_start()..=self.get_stop()
    }

    /// The length of this indel event. For Insertions, how long the inserted sequence is,
    /// and for deletions, how many bases are spanned by the deletion.
    pub fn len(&self) -> usize {
        match self {
            InDel::Ins(ins) => ins.sequence.len(),
            InDel::Del(del) => del.start.abs_diff(del.stop),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[pyclass]
pub struct Insertion {
    /// Base position directly to the left of the insertion in the forward sequence.
    #[pyo3(get)]
    position: usize,

    /// The sequence bytes that have been inserted to the right of the position.
    sequence: Vec<u8>,
}

#[pymethods]
impl Insertion {
    #[new]
    pub fn new(position: usize, sequence: Vec<u8>) -> Self {
        Self { position, sequence }
    }

    fn __repr__(&self) -> String {
        format!(
            "Insertion(position={}, sequence='{}')",
            self.position,
            self.py_sequence()
        )
    }

    #[getter]
    #[pyo3(name = "sequence")]
    fn py_sequence(&self) -> String {
        let seq = String::from_utf8_lossy(self.sequence.as_slice());
        format!("{seq}")
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[pyclass]
pub struct Deletion {
    /// Position of the first base that was affected by this deletion.
    #[pyo3(get)]
    start: usize,

    /// Position of the last base that was affected by this deletion.
    #[pyo3(get)]
    stop: usize,
}

#[pymethods]
impl Deletion {
    #[new]
    pub fn new(start: usize, stop: usize) -> Self {
        Self { start, stop }
    }

    fn __repr__(&self) -> String {
        format!("Deletion(start={}, stop={})", self.start, self.stop)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn indel_test() {
        todo!("Write Insertion and Deletion tests!")
    }

    #[test]
    fn indel_interference() {
        todo!("Write Insertion and Deletion tests!")
    }

    #[test]
    fn indel_range() {
        todo!("Write Insertion and Deletion tests!")
    }
}

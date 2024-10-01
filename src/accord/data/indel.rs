use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::ops::RangeInclusive;

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Hash)]
pub enum InDel {
    Ins(Insertion),
    Del(Deletion),
}

impl InDel {
    /// Get the starting position of this indel event, independent of read direction.
    /// That means the left side start of the event when considering forward direction.
    pub fn get_start(&self) -> usize {
        match self {
            InDel::Ins(ins) => { ins.position }
            InDel::Del(del) => { del.start }
        }
    }

    pub fn get_stop(&self) -> usize {
        match self {
            InDel::Ins(ins) => { ins.position + 1 }
            InDel::Del(del) => { del.stop }
        }
    }

    /// The length of this indel event. For Insertions, how long the inserted sequence is,
    /// and for deletions, how many bases are spanned by the deletion.
    pub fn len(&self) -> usize {
        match self {
            InDel::Ins(ins) => { ins.sequence.len() }
            InDel::Del(del) => { del.start.abs_diff(del.stop) }
        }
    }

    /// Get a byte slice corresponding to this events sequence.
    /// I.e. what should be inserted between event start and stop.
    pub fn get_seq(&self) -> &[u8] {
        match self {
            InDel::Ins(ins) => { ins.sequence.as_slice() }
            InDel::Del(_) => { &[] }
        }
    }

    /// Base positions spanning the event site as an inclusive range `start..=stop`.
    /// Start and stop are independent of read direction, and you may assume order `start <= stop`.
    pub fn range(&self) -> RangeInclusive<usize> { self.get_start()..=self.get_stop() }

    // /// InDels interfere with one another, if they cover the same area.
    // /// Or in mixed cases (i.e. insertion and deletion), even if they're just directly adjacent.
    // fn interferes_with(&self, other: &impl InDel) -> bool;

    // /// Apply the InDel to a sequence. This may or may not change the sequence in-place.
    // fn apply(&self, seq: &mut String) -> &String;

    /// Whether this indel preserves the reading frame by only shifting it by a multiple of three.
    pub fn preserves_reading_frame(&self) -> bool { self.len() % 3 == 0 }

    /// Whether this indel breaks the reading frame,  by shifting it by a non-multiple of three.
    pub fn breaks_reading_frame(&self) -> bool { !self.preserves_reading_frame() }
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Hash)]
pub struct Insertion {
    /// Base position directly to the left of the insertion in the forward sequence.
    pub position: usize,

    /// The sequence bytes that have been inserted to the right of the position.
    sequence: Vec<u8>,
}

impl Insertion {
    pub fn new(position: usize, sequence: Vec<u8>) -> Self { Self { position, sequence } }
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Hash)]
pub struct Deletion {
    /// Position of the first base that was affected by this deletion.
    start: usize,

    /// Position of the last base that was affected by this deletion.
    stop: usize,
}

impl Deletion {
    pub fn new(start: usize, stop: usize) -> Self { Self { start, stop } }
}


#[cfg(test)]
mod tests {
    #[test]
    fn indel_test() {
        todo!("Write Insertion and Deletion tests!")
    }
    fn indel_interference() {
        todo!("Write Insertion and Deletion tests!")
    }
    fn indel_range() {
        todo!("Write Insertion and Deletion tests!")
    }
}

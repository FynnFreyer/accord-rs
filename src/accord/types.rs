use crate::accord::data::indel::InDel;
use counter::Counter;
use std::collections::HashMap;

/// A list of base counts for every position in the reference sequence.
pub type BaseCounts = Vec<Counter<u8>>;

/// A mapping from base characters to coverage for the respective base, relative to the reference sequence.
pub type ExpandedBaseCounts = HashMap<char, Coverage>;

/// A map in which encountered insertions point to their respective number of occurrences.
pub type InDelCounts = Counter<InDel, usize>;

/// Vector containing coverage of a reference genome per base position.
pub type Coverage = Vec<usize>;

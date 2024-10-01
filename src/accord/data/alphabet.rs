use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::ops::Index;

#[derive(Serialize, Deserialize, Debug)]
pub struct Alphabet {
    set: HashSet<char>,
    gap: char,
    idx: HashMap<char, usize>,
}

impl Alphabet {
    pub fn new(symbols: &str, gap: char) -> Self {
        let set = symbols
            .chars()
            .collect::<HashSet<char>>();
        let idx = Self::index_alphabet(&set, &gap);

        Self { set, gap, idx }
    }

    pub fn size(&self) -> usize {
        self.idx.len()
    }

    fn index_alphabet<'a>(alphabet: &HashSet<char>, gap: &char) -> HashMap<char, usize> {
        // create lookup table for symbol positions
        let mut sym_index = HashMap::new();
        for (symbol_pos, symbol) in alphabet.iter().enumerate() {
            sym_index.insert(symbol.clone(), symbol_pos);
        }
        sym_index.insert(gap.clone(), sym_index.len());
        sym_index
    }
}

impl Index<char> for Alphabet {
    type Output = usize;

    fn index(&self, index: char) -> &Self::Output {
        if !self.set.contains(&index) && index != self.gap {
            panic!("Char '{index}' not contained in alphabet")
        }
        self.idx.index(&index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alphabet_from_string_fills_set() {
        let nucleotides = "ACGT";
        let alphabet = Alphabet::new(&nucleotides, '-');
        assert_eq!(alphabet.set.len(), nucleotides.len())
    }

    #[test]
    fn alphabet_from_string_is_indexed() {
        let nucleotides = "ACGT";
        let alphabet = Alphabet::new(&nucleotides, '-');
        assert_eq!(alphabet.idx.len(), nucleotides.len() + 1)
    }

    #[test]
    fn alphabet_produces_template() {
        todo!("Move this to BaseCount tests");
        let nucleotides = "ACGT";
        let alphabet = Alphabet::new(&nucleotides, '-');
        let a_pos = alphabet.idx[&'A'];

        let reference = String::from("GATTACA");
        let base_counts = alphabet.get_base_count_template(&reference);
        assert_eq!(base_counts.len(), nucleotides.len() + 1);
        assert_eq!(base_counts[a_pos].len(), reference.len());
    }

    #[test]
    fn alphabet_is_indexable() {
        let nucleotides = "ACGT";
        let alphabet = Alphabet::new(&nucleotides, '-');
        for sym in nucleotides {
            let pos = alphabet.idx[&'A'];
            assert!(pos >= 0)
        }
    }
}


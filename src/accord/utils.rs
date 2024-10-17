//! Module containing utility functions for file IO.

use super::data::seq::Seq;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

pub fn write_file(content: &String, path: &str) -> () {
    //! Helper function to write files.

    let path = Path::new(path);
    let loc = path.display();  // used in error messages

    // open file for writing
    let mut file = match File::create_new(&path) {
        Err(e) => panic!("Failed to create {loc}: {e}"),
        Ok(f) => f,
    };

    match file.write_all(content.as_bytes()) {
        Err(e) => panic!("Failed to read {loc}: {e}"),
        Ok(_) => {}
    }
}

pub fn read_file(path: &str) -> String {
    //! Helper function to read files.

    let path = Path::new(path);
    let loc = path.display();  // used in error messages

    // open file for reading
    let mut file = match File::open(&path) {
        Err(e) => panic!("Failed to open {loc}: {e}"),
        Ok(f) => f,
    };

    // read file
    let mut raw = String::new();
    match file.read_to_string(&mut raw) {
        Err(e) => panic!("Failed to read {loc}: {e}"),
        Ok(_) => {}
    }

    raw
}

pub fn get_fasta_seq(fasta_path: &String) -> Seq {
    let ref_string = read_file(fasta_path.as_str());

    let seqs = Seq::from_fasta(ref_string);
    let ref_seq = match seqs.first() {
        None => panic!("No sequence found in file: {}", fasta_path),
        Some(seq) => seq.clone()
    };

    ref_seq
}


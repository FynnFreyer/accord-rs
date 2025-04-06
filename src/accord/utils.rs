//! Module containing utility functions for file IO.

use std::fs::File;
use std::io::Write;
use std::path::Path;

pub fn write_file(content: &String, path: &str) -> () {
    //! Helper function to write files.

    let path = Path::new(path);
    let loc = path.display(); // used in error messages

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

pub fn change_suffix(path: &String, new_suffix: &str) -> String {
    let path = Path::new(path);
    let new_path = path.with_extension(new_suffix);
    new_path.to_str().unwrap().to_string()
}

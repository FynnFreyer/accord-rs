mod accord;
mod utils;

use pyo3::prelude::*;

use accord::calculator::Calculator;

/// The internals of `accord-rs` that are implemented in Rust.
#[pymodule(name = "_internal")]
mod py_accord {
    use super::*;

    #[pymodule_export]
    use Calculator;
}

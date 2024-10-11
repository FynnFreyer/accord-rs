mod accord;
mod utils;

use pyo3::prelude::*;

use accord::calculator::Calculator;

/// The internals of `accord-rs` that are implemented in Rust.
#[pymodule(name = "_internal")]
fn py_accord(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Calculator>()?;
    Ok(())
}

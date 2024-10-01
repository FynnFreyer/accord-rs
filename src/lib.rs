mod accord;
mod utils;

use accord::App;
use accord::calculator::Calculator;
use utils::get_fasta_seq;

use pyo3::prelude::*;

#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyfunction]
fn build_consensus(ref_path: String, aln_path: String) -> PyResult<String> {
    let ref_seq = get_fasta_seq(&ref_path);
    let consensus = App::calculate_consensus(ref_seq, aln_path);
    Ok(consensus.sequence_string())
}

/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name="_accord")]
fn py_accord(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(build_consensus, m)?)?;
    m.add_class::<Calculator>()?;
    Ok(())
}

mod accord;
mod utils;

use pyo3::prelude::*;

use accord::calculator;
use accord::data;

/// The internals of `accord-rs`. These are implemented in Rust.
#[pymodule(name = "_internal")]
mod py_accord {
    use super::*;

    #[pymodule_export]
    use calculator::Calculator;

    #[pymodule(name = "data")]
    mod py_data {
        use super::*;

        #[pymodule_export]
        use data::seq::Seq;

        #[pymodule_export]
        use data::settings::AlnQualityReqs;

        #[pymodule(name = "indel")]
        mod py_indel {
            use super::*;

            #[pymodule_export]
            use data::indel::Insertion;
            #[pymodule_export]
            use data::indel::Deletion;
            #[pymodule_export]
            use data::indel::InDel;
        }
    }
}

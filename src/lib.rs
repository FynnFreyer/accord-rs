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

    /// Holds classes for handling sequence data, settings, results etc.
    #[pymodule(name = "data")]
    mod py_data {
        use super::*;
        use data::seq;
        use data::settings;
        use data::stats;

        #[pymodule_export]
        use seq::Seq;

        #[pymodule_export]
        use settings::AlnQualityReqs;

        /// Classes for working with InDels.
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

        /// Classes for working with statistical data.
        #[pymodule(name = "stats")]
        mod py_stats {
            use super::*;

            #[pymodule_export]
            use stats::Quantile;
            #[pymodule_export]
            use stats::DistributionStats;
            #[pymodule_export]
            use stats::AlnStats;
        }
    }
}

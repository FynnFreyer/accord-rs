//! This module provides structs for working with, and summarizing statistical data concerning aligned reads.

use itertools::Itertools;
use log::warn;
use pyo3::types::PyType;
use pyo3::{pyclass, pymethods, Bound};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;


/// Relevant data for an aligned read.
#[derive(Debug, Clone)]
#[pyclass]
pub struct AlnData {
    /// Length of the aligned sequence.
    #[pyo3(get)]
    length: usize,

    /// Quality of the alignment.
    #[pyo3(get)]
    mapq: u8,

    /// SAM flags of the alignment.
    #[pyo3(get)]
    flags: u16,

    /// Alignment score from 'AS' field.
    #[pyo3(get)]
    score: usize,

    /// Editing distance from 'NM' field.
    #[pyo3(get)]
    distance: usize,
}

impl AlnData {
    /// Create an `AlnData` object by extracting relevant data from an htslib BAM record.
    pub fn from_record(record: &Record) -> Self {
        let length = record.seq_len();
        let mapq = record.mapq();
        let flags = record.flags();
        let score = Self::extract_unisgned(record, b"AS");
        let distance = Self::extract_unisgned(record, b"NM");

        Self { length, mapq, flags, score, distance }
    }

    fn extract_unisgned(record: &Record, tag: &[u8]) -> usize {
        let _tag_name = String::from_utf8_lossy(tag);
        match record.aux(tag) {
            Ok(value) => {
                if let Aux::U8(v) = value {
                    v as usize
                } else if let Aux::U16(v) = value {
                    v as usize
                } else if let Aux::U32(v) = value {
                    v as usize
                } else {
                    panic!("Value in field '{_tag_name}' is not of an unsigned type: {value:?}")
                }
            }
            Err(_e) => panic!("Error extracting value from field '{_tag_name}': {_e}")
        }
    }
}

#[pymethods]
impl AlnData {
    fn __repr__(&self) -> String {
        format!("AlnData(length={}, mapq={}, flags={}, score={}, distance={})",
                self.length, self.mapq, self.flags, self.score, self.distance)
    }
}

/// Struct for quantile data of unsigned integral values.
///
/// Combines a quantile factor and the quantile value.
/// E.g., with `{ factor: 0.2, value: 3 }` 20 % of values are lower or equal to 3.
#[derive(Debug, Clone)]
#[pyclass]
pub struct Quantile {
    /// Floating point number that describes the percentage cutoff for this quantile.
    /// Must be in interval $[0, 1]$.
    #[pyo3(get)]
    factor: f64,

    /// The quantile value.
    #[pyo3(get)]
    value: usize,
}

#[pymethods]
impl Quantile {
    fn __repr__(&self) -> String {
        format!("Quantile(factor='{}', value='{}')", self.factor, self.value)
    }
}

/// Metrics describing the distribution of *unsigned, integral* data.
#[derive(Debug, Clone)]
#[pyclass]
pub struct DistStats {
    /// List of quantile values for the distribution.
    #[pyo3(get)]
    quantiles: Vec<Quantile>,

    /// The sample size $n$.
    #[pyo3(get)]
    sample_size: usize,

    /// The sample mean $\bar x$.
    #[pyo3(get)]
    mean: f64,

    /// The sum of squares $\sum_{i = 1}^{n}(x_i - \bar x)^2$.
    #[pyo3(get)]
    sum_of_squares: f64,
}

impl DistStats {
    /// Determine the distribution of some numbers.
    pub fn from_numbers(numbers: Vec<usize>, quantile_factors: &Vec<f64>) -> Self {
        let quantiles = Self::calculate_quants(&numbers, quantile_factors);
        let sample_size = numbers.len();

        let total = numbers.iter().sum::<usize>();
        let mean = total as f64 / sample_size as f64;

        let sum_of_squares = numbers.iter().map(|num| {
            (*num as f64 - mean).powi(2)
        }).sum();

        Self { quantiles, sample_size, mean, sum_of_squares }
    }

    fn calculate_quants(numbers: &Vec<usize>, factors: &Vec<f64>) -> Vec<Quantile> {
        let n = numbers.len();
        let sorted_nums = numbers.iter().sorted().collect_vec();

        if n < factors.len() {
            warn!("Trying to determine more quantiles than numbers in sequence.");
        }

        let mut quantiles = Vec::new();
        for factor in factors {
            // we determine the rank n, that is the n_th element of the sorted list, we're interested in
            let rank = (n as f64 * factor).round() as usize;

            // we subtract one to translate rank to index (first element -> index zero, etc.)
            // also, we make sure to not undershoot by subtracting, by enforcing a minimum of zero
            let index = rank.saturating_sub(1);

            // construct the quantile
            let factor = factor.clone();
            let value = sorted_nums[index].clone();
            let quantile = Quantile { factor, value };

            quantiles.push(quantile);
        }

        quantiles
    }
}

#[pymethods]
impl DistStats {
    #[classmethod]
    #[pyo3(name = "from_numbers")]
    fn py_from_numbers(_cls: &Bound<'_, PyType>, numbers: Vec<usize>, factors: Vec<f64>) -> Self {
        Self::from_numbers(numbers, &factors)
    }

    #[getter]
    pub fn std_deviation(&self) -> f64 {
        self.variance().sqrt()
    }

    #[getter]
    pub fn variance(&self) -> f64 {
        self.sum_of_squares / self.sample_size as f64
    }

    fn __repr__(&self) -> String {
        let mut quant_reprs = Vec::new();
        for quant in &self.quantiles {
            quant_reprs.push(quant.__repr__())
        }
        let mut quant_repr = String::from("[");
        quant_repr.push_str(quant_reprs.join(", ").as_str());
        quant_repr.push(']');

        format!(
            "DistStats(quantiles={quant_repr}, sample_size={}, mean={}, sum_of_squares={})",
            self.sample_size, self.mean, self.sum_of_squares
        )
    }
}

/// Statistical data of the seen alignments.
#[derive(Debug, Clone)]
#[pyclass]
pub struct AlnStats {
    /// Statistics describing the distribution of alignment lengths.
    #[pyo3(get)]
    length_distribution: DistStats,

    /// Statistics describing the distribution of alignment quality.
    #[pyo3(get)]
    quality_distribution: DistStats,

    /// Statistics describing the distribution of alignment scores.
    #[pyo3(get)]
    score_distribution: DistStats,

    /// Statistics describing the distribution of alignment editing distances.
    #[pyo3(get)]
    editing_distance_distribution: DistStats,
}

impl AlnStats {
    pub fn from_data(aln_data: &Vec<AlnData>, quantile_factors: &Vec<f64>) -> Self {
        let (
            length_distribution,
            quality_distribution,
            score_distribution,
            editing_distance_distribution
        ) = Self::calculate_distributions(aln_data, quantile_factors);

        Self {
            length_distribution,
            quality_distribution,
            score_distribution,
            editing_distance_distribution,
        }
    }

    fn calculate_distributions(aln_data: &Vec<AlnData>, quantile_factors: &Vec<f64>)
                               -> (DistStats, DistStats, DistStats, DistStats) {
        // split aln data into seperate vectors
        let mut lengths = Vec::new();
        let mut qualities = Vec::new();
        let mut scores = Vec::new();
        let mut distances = Vec::new();

        for data in aln_data {
            lengths.push(data.length);
            qualities.push(data.mapq as usize);
            scores.push(data.score);
            distances.push(data.distance);
        }

        // calculate dist stats
        let length_distribution = DistStats::from_numbers(lengths, quantile_factors);
        let quality_distribution = DistStats::from_numbers(qualities, quantile_factors);
        let score_distribution = DistStats::from_numbers(scores, quantile_factors);
        let editing_distance_distribution = DistStats::from_numbers(distances, quantile_factors);

        (length_distribution, quality_distribution, score_distribution, editing_distance_distribution)
    }
}

#[pymethods]
impl AlnStats {
    #[classmethod]
    #[pyo3(name = "from_data")]
    fn py_from_data(_cls: &Bound<'_, PyType>, data: Vec<AlnData>, factors: Vec<f64>) -> Self {
        Self::from_data(&data, &factors)
    }

    /// Number of reads that these statistics were generated from.
    #[getter]
    pub fn sample_size(&self) -> usize {
        self.length_distribution.sample_size
    }

    fn __repr__(&self) -> String {
        format!(
            "AlnStats(length_distribution={}, quality_distribution={}, \
            score_distribution={}, editing_distance_distribution={})",
            self.length_distribution.__repr__(), self.quality_distribution.__repr__(),
            self.score_distribution.__repr__(), self.editing_distance_distribution.__repr__(),
        )
    }
}

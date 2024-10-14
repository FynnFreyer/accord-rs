use itertools::Itertools;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use serde::{Deserialize, Serialize};


/// Relevant data for an aligned read.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct AlnData {
    length: usize,
    mapq: u8,
    flags: u16,
    score: usize,
    distance: usize,
}

impl AlnData {
    pub fn from_record(record: &Record) -> Self {
        let length = record.seq_len();
        let mapq = record.mapq();
        let flags = record.flags();
        let score = Self::extract_unisgned(record, b"AS");
        let distance = Self::extract_unisgned(record, b"NM");

        Self { length, mapq, flags, score, distance }
    }

    fn extract_unisgned(record: &Record, tag: &[u8]) -> usize {
        let tag_name = String::from_utf8_lossy(tag);
        match record.aux(tag) {
            Ok(value) => {
                if let Aux::U8(v) = value {
                    v as usize
                } else if let Aux::U16(v) = value {
                    v as usize
                } else if let Aux::U32(v) = value {
                    v as usize
                } else {
                    panic!("Value in field '{tag_name}' is not of an unsigned type: {value:?}")
                }
            }
            Err(e) => panic!("Error extracting value from field '{tag_name}': {e}")
        }
    }
}

/// Struct for quantile data of unsigned integral values.
///
/// Combines a quantile factor and the quantile value.
/// E.g., with `{ factor: 0.2, value: 3 }` 20 % of values are lower or equal to 3.
#[derive(Serialize, Deserialize, Debug, Clone)]
struct Quantile {
    factor: f64,
    value: usize,
}

/// Metrics describing the distribution of *unsigned, integral* data.
#[derive(Serialize, Deserialize, Debug, Clone)]
struct DistributionStats {
    quantiles: Vec<Quantile>,
    sample_size: usize,
    mean: f64,
    sum_of_squares: f64,
}

impl DistributionStats {
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

    pub fn std_deviation(self) -> f64 {
        self.variance().sqrt()
    }

    pub fn variance(self) -> f64 {
        self.sum_of_squares / self.sample_size as f64
    }

    fn calculate_quants(numbers: &Vec<usize>, factors: &Vec<f64>) -> Vec<Quantile> {
        let n = numbers.len();
        let sorted_nums = numbers.iter().sorted().collect_vec();

        if n < factors.len() {
            panic!("Trying to determine more quantiles than numbers in sequence.")
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

/// Statistical data of the seen alignments.
/// The distribution stats
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct AlnStats {
    length_distribution: DistributionStats,
    quality_distribution: DistributionStats,
    score_distribution: DistributionStats,
    editing_distance_distribution: DistributionStats,
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

    pub fn get_sample_size(&self) -> usize {
        self.length_distribution.sample_size
    }

    fn calculate_distributions(aln_data: &Vec<AlnData>, quantile_factors: &Vec<f64>) -> (
        DistributionStats,
        DistributionStats,
        DistributionStats,
        DistributionStats,
    ) {
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

        let length_distribution = DistributionStats::from_numbers(lengths, quantile_factors);
        let quality_distribution = DistributionStats::from_numbers(qualities, quantile_factors);
        let score_distribution = DistributionStats::from_numbers(scores, quantile_factors);
        let editing_distance_distribution = DistributionStats::from_numbers(distances, quantile_factors);

        (length_distribution, quality_distribution, score_distribution, editing_distance_distribution)
    }
}

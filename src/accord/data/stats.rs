use std::iter::Sum;

use itertools::Itertools;
use num::traits::real::Real;
use num::traits::AsPrimitive;
use num::{Num, Saturating};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use serde::{Deserialize, Serialize};


#[derive(Serialize, Deserialize, Debug)]
pub struct AlnData {
    length: usize,
    mapq: u8,
    flags: u16,
    score: u8,
    distance: u8,
}

impl AlnData {
    pub fn from_record(record: &Record) -> Self {
        let length = record.seq_len();
        let mapq = record.mapq();
        let flags = record.flags();
        let score = Self::extract_u8(record, b"AS");
        let distance = Self::extract_u8(record, b"NM");

        Self { length, mapq, flags, score, distance }
    }

    fn extract_u8(record: &Record, tag: &[u8]) -> u8 {
        let tag_name = String::from_utf8_lossy(tag);
        match record.aux(tag) {
            Ok(value) => {
                if let Aux::U8(v) = value {
                    v
                } else {
                    panic!("Value in field '{tag_name}' is not of type u8: {value:?}")
                }
            }
            Err(e) => panic!("Error extracting field '{tag_name}': {e}")
        }
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct Quantile<T> {
    factor: f64,
    value: T,
}

/// Distribution of *integral* data.
#[derive(Serialize, Deserialize, Debug)]
struct DistributionStats<T> {
    quantiles: Vec<Quantile<T>>,
    total: usize,
    sample_size: usize,
    sum_of_squares: f64,
    mean: f64,
}

impl<T: Num + AsPrimitive<f64> + AsPrimitive<usize> + Sum<T> + Ord + Clone> DistributionStats<T> {
    pub fn from_numbers(numbers: Vec<T>) -> Self {
        let quantile_factors = vec![0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.0];
        let quantiles = Self::calculate_quants(&numbers, quantile_factors);
        let sample_size = numbers.len();
        let total = numbers.iter().fold(0, |accum, num| {
            let num_u: usize = num.as_();
            accum + num_u
        });
        let total_f: f64 = total.as_();
        let mean = total_f / sample_size as f64;
        let sum_of_squares = numbers.iter().map(|num| {
            let num_f: f64 = num.as_();
            let diff: f64 = num_f - mean;
            diff.powi(2)
        }).sum();

        Self { quantiles, total, sample_size, sum_of_squares, mean }
    }

    pub fn std_deviation(self) -> f64 {
        self.variance().sqrt()
    }

    pub fn variance(self) -> f64 {
        self.sum_of_squares / self.sample_size as f64
    }

    fn calculate_quants(numbers: &Vec<T>, factors: Vec<f64>) -> Vec<Quantile<T>> {
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
            let value = sorted_nums[index].clone();
            let quantile = Quantile { factor, value };
            quantiles.push(quantile);
        }

        quantiles
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct AlnStats {
    length_distribution: DistributionStats<usize>,
    quality_distribution: DistributionStats<u8>,
    score_distribution: DistributionStats<u8>,
    editing_distance_distribution: DistributionStats<u8>,
    total_reads: usize,
    mapped_reads: usize,
}

impl AlnStats {
    pub fn from_data(aln_data: Vec<AlnData>) -> Self {
        let (
            length_distribution,
            quality_distribution,
            score_distribution,
            editing_distance_distribution
        ) = Self::calculate_distributions(aln_data);

        // TODO consider mapped and total reads
        Self { length_distribution, quality_distribution, score_distribution, editing_distance_distribution, mapped_reads: 0, total_reads: 0 }
    }

    fn calculate_distributions(aln_data: Vec<AlnData>) -> (
        DistributionStats<usize>,
        DistributionStats<u8>,
        DistributionStats<u8>,
        DistributionStats<u8>
    ) {
        let mut lengths = Vec::new();
        let mut qualities = Vec::new();
        let mut scores = Vec::new();
        let mut distances = Vec::new();

        for data in aln_data {
            lengths.push(data.length);
            qualities.push(data.mapq);
            scores.push(data.score);
            distances.push(data.distance);
        }

        let length_distribution = DistributionStats::from_numbers(lengths);
        let quality_distribution = DistributionStats::from_numbers(qualities);
        let score_distribution = DistributionStats::from_numbers(scores);
        let editing_distance_distribution = DistributionStats::from_numbers(distances);

        (length_distribution, quality_distribution, score_distribution, editing_distance_distribution)
    }
}

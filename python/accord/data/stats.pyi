from typing import Sequence


class AlnData:
    length: int
    mapq: int
    flags: int
    score: int
    distance: int


class Quantile:
    factor: float
    value: int


class DistStats:
    quantiles: Sequence[Quantile]
    sample_size: int
    mean: float
    sum_of_squares: float
    std_deviation: float
    variance: float


class AlnStats:
    length_distribution: DistStats
    quality_distribution: DistStats
    score_distribution: DistStats
    editing_distance_distribution: DistStats

    def from_data(self, data: Sequence[AlnData], factors: Sequence[float], total_reads: int) -> "AlnStats": ...

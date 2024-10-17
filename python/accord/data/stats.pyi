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
    quantiles: list[Quantile]
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
    total_reads: int
    mapped_reads: int
    unmapped_reads: int

    def from_data(self, data: list[AlnData], factors: list[float], total_reads: int) -> "AlnStats": ...

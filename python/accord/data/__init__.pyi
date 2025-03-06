from typing import Sequence, Set, Mapping

from .stats import AlnData, AlnStats


class AlnQualityReqs:
    min_mapq: int
    mandatory_flags: int
    prohibited_flags: int
    indel_cutoff: float
    save_ends: int
    min_observations: int

    def __init__(self, min_mapq: int, mandatory_flags: int, prohibited_flags: int,
                 indel_cutoff: float, save_ends: int, min_observations: int): ...


class Seq:
    label: str
    sequence: str

    def __init__(self, label: str, sequence: str): ...

    @classmethod
    def from_fasta(cls, fasta: str) -> Sequence["Seq"]: ...

    @classmethod
    def from_file(cls, file: str) -> Sequence["Seq"]: ...

    def to_fasta(self) -> str: ...


class AnalysisResult:
    coverage: Sequence[int]
    valid_alns: Sequence[AlnData]
    reads_seen: Set[int]


class Consensus:
    ref_seq: Seq
    aln_path: str
    consensus_seq: Seq
    aln_stats: AlnStats
    coverage: Sequence[int]
    base_counts: Mapping[str, Sequence[int]]
    total_reads: int
    valid_reads: int
    invalid_reads: int

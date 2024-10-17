from .data import Seq, AlnQualityReqs
from .data.stats import AlnData, AlnStats


class Calculator:
    ref_seq: Seq
    aln_path: str
    aln_quality_reqs: AlnQualityReqs
    coverage: list[int]
    aln_data: list[AlnData]
    reads_seen: set[int]

    def __init__(self, ref_path: str, aln_path: str, reqs: AlnQualityReqs): ...

    def calculate(self) -> tuple[Seq, AlnStats]: ...

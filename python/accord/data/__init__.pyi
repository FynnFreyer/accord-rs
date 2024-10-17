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
    def from_fasta(cls, fasta: str) -> list["Seq"]: ...

    def to_fasta(self) -> str: ...

# accord-rs

To bring something into accord (*/əˈkɔrd/*) is to make it agree or correspond.
The `accord-rs` library calculates a consensus from reads that where aligned against a reference sequence.

This library is in the early alpha stage, and its API may be subject to change without notice.

## Python

Accord provides bindings for consensus calculation from Python code.
Install with `pip install accord-rs`, and then use like so:

```python
from accord import Calculator
from accord.data import AlnQualityReqs, Seq

# settings for consensus calculation
reqs = AlnQualityReqs(
    min_mapq=10,            # minimum mapping quality
    mandatory_flags=0,      # required SAM flags, see: https://broadinstitute.github.io/picard/explain-flags.html
    prohibited_flags=1540,  # no unmapped, failing quality or optical duplicate reads
    indel_cutoff=0.1,       # only indels contained in at least 10 % of reads covering that position are considered
    save_ends=0,            # this has no function yet - only relevant when using PCR fragments
    min_observations=50,    # base coverage needs to be at least 50
)

# generate the consensus and alignment statistics
calc = Calculator(reqs)
consensus = calc.calculate("/path/to/reference.fasta", "/path/to/aln.bam")

aln_path = "/home/fynn/Desktop/HIV/HIV-Pipeline/out/02_mapping/alignments/200626_20-07149_HIV20-00278_F483_9_MiS85_S1_L000_001.cyc01.bam"
ref_path = "/home/fynn/Desktop/HIV/HIV-Pipeline/tst/assets/references/HXB2_K03455.fasta"
# print the consensus as a FASTA record
print(consensus.consensus_seq.to_fasta())
```

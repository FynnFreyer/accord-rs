# Data

Test data.

## Generating Alignments

Index the reference, then run `bwa mem` and `samtools sort`.
E.g. for the full-length data:

```
bwa index K03455.fasta
bwa mem K03455.fasta reads/FL* | samtools sort > alns/FL.bam
```

# Quality control of the sequencing data.

Several tools available to do so. For this workshop, we will use fastqc.

```
fastqc -o QC/ -f WT1_R1.fastq.gz

```

## How to act on fastq after QC.

We can do several trimming:

  * on quality using Phred score: we want an accuracy of 99%. What will be the Phred score?
  * on the sequences, if they contain adaptor sequences.

To do so, we can use on tools: cutadapt.

```
cutadapt -q 20 -a file:QC/adaptors.fa -A file:QC/adaptors.fa -o trimmed/WT1_R1_trimmed.fastq -p trimmed/WT1_R2_trimmed.fastq sequencing/WT1_R1.fastq.gz sequencing/WT1_R2.fastq.gz

```

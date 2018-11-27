## Getting a count matrix:

To get a matrix of count we will use [featureCounts ](https://www.ncbi.nlm.nih.gov/pubmed/23558742)

```
mkdir count

cd mapping/

featureCounts -a /mnt/RNAseq_Workshop_Data/genome/Ensembl_NCBIM37.67_chr5.gtf -o ../count/CountMat_chr5.tsv  -t exon \
              -g gene_id -Q 10 -s 1 -p -T 8 WT1_sorted.bam WT2_sorted.bam WT3_sorted.bam KO1_sorted.bam KO2_sorted.bam KO3_sorted.bam

```

Important options:

  * -a annotation file
  * -t What is the feature used to count the overlap
  * -g feature to summarise against.
  * -s strandness of the hit: 0: unstrand, 1: forward stranded, 2: reverse stranded.
  * -p count paired reads as __one fragment__

## featureCounts output:

First, we need to see how many output does featureCounts produce:

```
ls ../count/
CountMat_chr5.tsv 
CountMat_chr5.tsv.summary
```

So 2 files:

  * Count matrix: CountMat_chr5.tsv
  * Summary of the counting: CountMat_chr5.tsv.summary

```
head ../count/CountMat_chr5.tsv

#or to have Excel look and feel:

gnumeric ../count/CountMat_chr5.tsv 
```


\# Program:featureCounts v1.6.3; Command:"featureCounts" "-a" "/mnt/RNAseq_Workshop_Data/genome/Ensembl_NCBIM37.67_chr5.gtf" "-o" "../count/CountMat_chr5.tsv" "-t" "exon" "-g" "gene_id" "-Q" "10" "-s" "1" "-p" "-T" "8" "WT1_sorted.bam" "WT2_sorted.bam" "WT3_sorted.bam" "KO1_sorted.bam" "KO2_sorted.bam" "KO3_sorted.bam"
Geneid | Chr | Start | End | Strand | Length | WT1_sorted.bam | WT2_sorted.bam | WT3_sorted.bam | KO1_sorted.bam | KO2_sorted.bam | KO3_sorted.bam
------ | --- | ----- | --- | ------ | ------ | -------------- | -------------- | -------------- | -------------- | -------------- | --------------
ENSMUSG00000090577 | chr5;chr5;chr5 | 3030238;3031011;3031960 | 3030381;3031194;3032082 | +;+;+ | 451 | 0 | 0 | 0 | 0 | 0 | 0
------------------ | -------------- | ----------------------- | ----------------------- | ----- | --- | - | - | - | - | - | -

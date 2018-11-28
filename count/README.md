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

File header:

\# Program:featureCounts v1.6.3; Command:"featureCounts" "-a" "/mnt/RNAseq_Workshop_Data/genome/Ensembl_NCBIM37.67_chr5.gtf" "-o" "../count/CountMat_chr5.tsv" "-t" "exon" "-g" "gene_id" "-Q" "10" "-s" "1" "-p" "-T" "8" "WT1_sorted.bam" "WT2_sorted.bam" "WT3_sorted.bam" "KO1_sorted.bam" "KO2_sorted.bam" "KO3_sorted.bam"



Geneid | Chr | Start | End | Strand | Length | WT1_sorted.bam | WT2_sorted.bam | WT3_sorted.bam | KO1_sorted.bam | KO2_sorted.bam | KO3_sorted.bam
------ | --- | ----- | --- | ------ | ------ | -------------- | -------------- | -------------- | -------------- | -------------- | --------------
ENSMUSG00000090577 | chr5;chr5;chr5 | 3030238;3031011;3031960 | 3030381;3031194;3032082 | +;+;+ | 451 | 0 | 0 | 0 | 0 | 0 | 0

Now, let's have a look at the summary file:

```
cat CountMat_chr5.tsv.summary

Status	WT1_sorted.bam	WT2_sorted.bam	WT3_sorted.bam	KO1_sorted.bam	KO2_sorted.bamKO3_sorted.bam
Assigned	11793	9508	10617	9639	10278	9653
Unassigned_Unmapped	80	46	82	93	58	64
Unassigned_MappingQuality	89	31	48	60	46	43
Unassigned_Chimera	0	0	0	0	0	0
Unassigned_FragmentLength	0	0	0	0	0	0
Unassigned_Duplicate	0	0	0	0	0	0
Unassigned_MultiMapping	23248	15861	22103	13815	29572	20306
Unassigned_Secondary	0	0	0	0	0	0
Unassigned_Nonjunction	0	0	0	0	0	0
Unassigned_NoFeatures	96383	88984	106545	88874	56047	69737
Unassigned_Overlapping_Length	0	0	0	0	0	0
Unassigned_Ambiguity	95	63	71	78	81	79
```

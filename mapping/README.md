#How to align sequencing reads to genome

## Preparation of the genome

To be able to map (align) sequencing reads on the genome, the genome needs to be indexed first. In this workshop we will use [HISAT2](https://www.nature.com/articles/nmeth.3317).
Note for speed reason, the reads will be aligned on the chr5 of the mouse genome.

```
# need to go where the sequence is stored:
cd genome
#to list what is in your directory:
ls
chr5.fa
# only chromosome 5 available in this directory

# index file:
hisat2-build -p 2 -f chr5.fa chr5

# take around 3 mins
#list what is in the directory:
ls
chr5.1.ht2  chr5.3.ht2  chr5.5.ht2  chr5.7.ht2  chr5.fa
chr5.2.ht2  chr5.4.ht2  chr5.6.ht2  chr5.8.ht2

```
How many files were created during the indexing process?

## Alignment on the genome (chromosome5)

Now that the genome is prepared. Sequencing reads can be aligned.

Information required:

  * Where the sequence information is stored (e.g. fastq files ...) ?
  * What kind of sequencing: Single End or Paired end ?
  * Where are stored the indexes and the genome? 
  * Where will the mapping files be stored?

```
# syntax:
hisat2 -x genome/chr5 -1 sequencing/WT1_R1.fastq.gz -2 sequencing/WT1_R1.fastq.gz -S mapping/mapping_WT1.sam

```
Now we need to align all the rest of the samples.

Hint: looping over the files?


```
# let's use a counter to process our samples:
# for each iteration $i will change:  first iteration i=1, then 2 then 3.
for i in `seq 1 3`; 
do
# processing WT first:
hisat2 -x genome/chr5 -1 sequencing/WT$i\_R1.fastq.gz -2 sequencing/WT$\i_R2.fastq.gz -S mapping/WT$i.sam
# processing KO
hisat2 -x genome/chr5 -1 sequencing/KO$i\_R1.fastq.gz -2 sequencing/KO$\i_R2.fastq.gz -S mapping/KO$i.sam

done

```

Now we can explore our SAM files.

## Converting SAM files to BAM files

This SAM to BAM files can be achieved by several tools. Today we will focus on samtools which the most commonly used.

To convert sam to bam 

```
samtools view -b mapping/WT1.sam -o mapping/WT1.bam

```

We need to perform this step for all the samples. 

Hint: loop???

```
# let's re-use a counter to process our samples:
# for each iteration $i will change:  first iteration i=1, then 2 then 3.
for i in `seq 1 3`;
do
# processing WT first:
samtools view -b mapping/WT$i.sam -o  mapping/WT$i.bam
# processing KO
samtools view -b mapping/KO$i.sam -o  mapping/KO$i.bam
done

```

## Getting a count matrix:

To get a matrix of count we will use [featureCounts ](https://www.ncbi.nlm.nih.gov/pubmed/23558742)

```
cd mapping/

featureCounts -a genome/Ensembl_mm9_chr5.gtf -o ../count/CountMat_chr5.tsv  -t exon \
              -g gene_id -Q 10 -s 2 -p -T 8 WT1.bam WT2.bam WT3.bam KO1.bam KO2.bam KO3.bam

```

Important options:

  * -a annotation file 
  * -t What is the feature used to count the overlap
  * -g feature to summarise against.
  * -s strandness of the hit: 0: unstrand, 1: forward stranded, 2: reverse stranded.
  * -p count paired reads as __one fragment__

Other useful options (not used in this workshop)

  * -J Count number of reads that span exon-exon junctions
  * -L long reads



# Getting the data

## Starting material

The workshop will be based on the data published by Corley et al 2016.

Basically the experimental design has a total of 6 samples. 

  * 3 biological replicates of Gtf2ird1 knockout out mouse (KO) located on Chr5.
  * 3 biological replicates of Wild Type (WT).

## Accessing the data:

  1. Pubmed ID: [27295951](https://www.ncbi.nlm.nih.gov/pubmed/27295951)

  2. Accession ID: [GSE81082](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81082)

  3. SRA: [SRP074330](https://www.ncbi.nlm.nih.gov/sra?term=SRP074330)

  4. Click on [Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_1485027_130.14.18.97_5555_1542666014_3715348358_0MetA0_S_HStore&query_key=5)

  5. Save RunInfo table: SRA_list.tsv.

## RunInfo table

|                                                                                                                                                         | 
|---------------------------------------------------------------------------------------------------------------------------------------------------------| 
| SRR3473989      SAMN04939190    GSM2142347      192     SRX1741375      2016-05-03      5,064   2,367   Gtf2ird1 KO     Lip tissue, Gtf2ird1-KO         | 
|         SRR3473988      SAMN04939189    GSM2142346      193     SRX1741374      2016-05-03      6,692   3,139   Gtf2ird1 KO     Lip tissue, Gtf2ird1-KO | 
|         SRR3473987      SAMN04939188    GSM2142345      192     SRX1741373      2016-11-30      5,022   2,344   Gtf2ird1 KO     Lip tissue, Gtf2ird1-KO | 
|         SRR3473986      SAMN04939187    GSM2142344      192     SRX1741372      2016-05-03      5,751   2,650   Wild type       Lip tissue, WT          | 
|         SRR3473985      SAMN04939186    GSM2142343      193     SRX1741371      2016-05-03      5,212   2,397   Wild type       Lip tissue, WT          | 
|         SRR3473984      SAMN04939185    GSM2142342      192     SRX1741370      2016-05-03      6,438   2,970   Wild type       Lip tissue, WT          | 


RunInfo contains the meta-data required for the RNAseq analysis such as sequence name (SRA) and conditions.

## Downloading the fastq file

It needs to have SRA-Toolkit installed. [How to install SRA-Toolkit](https://ncbi.github.io/sra-tools/install_config.html).

Here is an example of bash script that download Paired end files using SRA ID stored in the 2nd column.

```

#!/bin/sh

# load appropriate modules
module purge
module load SRA-Toolkit

SRAfile=$PWD/SRA_list.tsv

outdir=$PWD/sra_paired_fastq
if [ -e $SRAfile ];then
   echo "getting data from $SRAfile"
   SRRs=`cat $SRAfile | cut -f2`
   for SRR in $SRRs; do
      echo "processing: $SRR"
      fastq-dump -A $SRR --split-3 --gzip -O $outdir
      echo fastq-dump done.....
   done
else
 echo missing $SRAfile
fi

```
Note: Output directory: sra_paired_fastq will contain 3 files per SRAID. Here is the nomenclature:
 
  * SRAID_1.fastq.gz (Read1)
  * SRAID_2.fastq.gz (Read2)
  * SRAID.fastq.gz (singleton)

## Subsetting the data

In order to make the hands-on quicker we will work on a subset of data, only reads that are mapping to chr5.

### Aligning to mouse genome.

#### Getting the mouse genome version mm9:

Several different places where mouse genome can be found. We will use the [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz)

```

# better create a directory to put the files to

mkdir mm9

cd mm9

# get it on the platform using wget:

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz

# or using curl if available on the system:

curl -O http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz

# This is tarball containing one fastA file per chromomosome.
# extracting the files:
gunzip -dc chromFa.tar.gz | tar -xof -

# It needs to be concatenated in one file:

cat chr*.fa > genome.fa

```

#### Indexing the genome 

Once the genome file is prepared then it needs to be indexed to use HISAT2. 

```
#!/bin/bash

module purge
module load HISAT2

hisat2-build -p 8 -f genome.fa genome

# note1 that -p is the number of threads that will be used. 
# note2 it should create 8 files with the suffix .ht2
```

#### Aligning on mm9

Finally aligning on mm9 :

```
#!/bin/bash

# load appropriate modules
module purge
module load HISAT2 SAMtools

Hisat2Idx=genome
SRAfile=$PWD/SRA_list.tsv

seqdir=$PWD/sra_paired_fastq
mappingdir=$PWD/mapping

if [ -d $mappingdir ] ; then
   echo "found mapping dir: $mappingdir"
else
   echo "creating mapping dir: $mappingdir"
   mkdir $mappingdir
fi

if [ -e $SRAfile ];then
   echo "getting data from $SRAfile"
   SRRs=`cat $SRAfile | cut -f2`
   for SRR in $SRRs; do
      echo "processing: $SRR"
      # only considering paired end reads:
      R1=$seqdir/$SRR\_1.fastq.gz
      R2=$seqdir/$SRR\_2.fastq.gz
      if [ -f $R1 ] && [ -f $R2 ]; then
         bam=$mappingdir/$SRR.bam
         hisat2 -x $Hisat2Idx  -p 8 -1 $R1 -2 $R2 | samtools view -Sb - | samtools sort -@ 8 - -o $bam 
         # align to genome convert to bam and sort by coordinates
      else
        echo Missing either $R1 or $R2 or both
      fi
   done
else
 echo missing $SRAfile
fi

```

### extracting reads mapped on chr5:

```
#!/bin/bash

# load appropriate modules
module purge
module load SAMtools

SRAfile=$PWD/SRA_list.tsv

seqdir=$PWD/sra_paired_fastq
mappingdir=$PWD/mapping
chr5dir=$PWD/chr5fastq
region="chr5:1-152537259"
if [ -d $chr5dir ] ; then
   echo "found mapping dir: $chr5dir"
else
   echo "creating mapping dir: $mappingdir"
   mkdir $chr5dir
fi

if [ -e $SRAfile ];then
   echo "getting data from $SRAfile"
   SRRs=`cat $SRAfile | cut -f2`
   for SRR in $SRRs; do
      echo "processing: $SRR"
      # only considering paired end reads:
      bam=$mappingdir/$SRR.bam
      if [ -f $bam ]; then
         # extracting only region and all paired end reads (note that we will take also singleton).
         chr5_R1=$chr5dir/$SRR\_chr5_R1.fastq.gz
	 chr5_R2=$chr5dir/$SRR\_chr5_R2.fastq.gz
         samtools view -b $bam $region | samtools fastq -c 6 -1 $chr5_R1 -2 $chr5_R2  -
      else
        echo Missing BAM file: $bam
      fi
   done
else
 echo missing $SRAfile
fi

```


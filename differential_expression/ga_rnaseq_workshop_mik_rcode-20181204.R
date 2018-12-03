## ------------------------------------------------------------------------
readcounts <- read.table("Data/CountMat_NCBIM37.67.dat", header = TRUE)

## ------------------------------------------------------------------------
dim(readcounts)
names(readcounts)

## ------------------------------------------------------------------------
head(readcounts)[,-c(2:6)]

## ------------------------------------------------------------------------
View( head(readcounts) )

## ------------------------------------------------------------------------
## Use the transcript IDs as the row names
row.names(readcounts) <- readcounts$Geneid
## Remove the non-count columns (and convert to a matrix)
counts <- as.matrix(readcounts[ , -c(1:6)])
## Remove the ".bam" suffix from teh name of each column
colnames(counts) <- gsub(".bam", "", colnames(counts))
## Show the first 6 rows of the counts object
head(counts)

## ------------------------------------------------------------------------
## Take logs of data (add 0.5 to avoid log(0) issues)
logCounts <- log2(counts + 0.5)

## ------------------------------------------------------------------------
## ## Boxplots of original and logged data
boxplot(counts ~ col(counts), names=colnames(counts))
boxplot(logCounts ~ col(logCounts), names=colnames(counts))

## Plot the density (distribution) for the first sample
plot( density(logCounts[,1]) )
## Use a loop to add the densities for the other samples
for(i in 2:6) lines( density(logCounts[,i]) )

## ------------------------------------------------------------------------
library(edgeR)
library(limma)

## ------------------------------------------------------------------------
## Create DGEList object - use counts and gene legnth information
dge <- DGEList(counts=counts, genes=data.frame( length=readcounts$Length) )
names(dge)
dge$counts[1:2,]
dge$genes[1:2,]

## ------------------------------------------------------------------------
dge$samples

## ------------------------------------------------------------------------
groups <- rep(c("WT", "KO"), c(3,3))
groups
design <- model.matrix(~groups)
design

## ------------------------------------------------------------------------
## Figure out what to keep: output is TRUE/FALSE for each gene
keep <- filterByExpr(dge, design)
## Make table of TRUE and FALSE
table(keep)
## Apply filtering and recalculate library sizes
dge <- dge[keep, keep.lib.sizes=FALSE]
## Calculate normalisation factor (i.e., account for total reads per sample)
dge <- calcNormFactors(dge)

## ------------------------------------------------------------------------
dge

## ------------------------------------------------------------------------
## Make log count data
logCounts <- log2(dge$counts + 0.5)
## Boxplots
boxplot(logCounts ~ col(logCounts), names=colnames(logCounts))
## Density plots
plot(density(logCounts[,1]))
for(i in 2:ncol(logCounts)) lines(density(logCounts[,i]), col=i)

## ------------------------------------------------------------------------
## Generate "counts per million", based on gene lengths
logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
head(logCPM)

## ------------------------------------------------------------------------
## Use genelength data from dge to calculate rpkm values (not on log scale)
rpkmData <- rpkm(dge)
head(rpkmData)

## ------------------------------------------------------------------------
## The function "t.test"" performs a two sample t-test in R.
## Perform t-test on first gene in data set (first row of matrix:
## first three values are WT, second three are KO):
t.test(logCPM[1,] ~ groups)

## ------------------------------------------------------------------------
## Very first gene does not appear to be 
## differentially expressed between the two groups:
t.test(logCPM[1,]~groups)$p.value

## How big is the change?
2^(2.886-2.544) 

## ------------------------------------------------------------------------
## Examine expression graphically using a dotplot:
lattice::dotplot(logCPM[1,]~groups, main=rownames(logCPM)[1],
           ylab='Log2 Expression')

## ------------------------------------------------------------------------
pvalues = c()
for(i in 1:10) pvalues[i] = t.test( logCPM[i,] ~ groups )$p.value
sort( round(pvalues,4) )

## ------------------------------------------------------------------------
round(p.adjust(sort(pvalues)),4)
round(p.adjust(sort(pvalues),"BH"),4)  # BH is the FDR method

## ------------------------------------------------------------------------
design

## ------------------------------------------------------------------------
## Load the limma package
library(limma)
## Fit linear model
fit = lmFit(logCPM, design)
fit = eBayes(fit)
tt = topTable(fit, coef=2, adjust="BH", n=nrow(logCPM))
options(digits=4)
tt[1:5,]

## ------------------------------------------------------------------------
topGene = match(rownames(tt)[1], rownames(logCPM))
t.test(logCPM[topGene,] ~ groups)
# "Raw" data agrees with logFC (WT - KO): -4.405 - 2.613 = -7.018

## ------------------------------------------------------------------------
lattice::dotplot(logCPM[topGene,] ~ groups, main=rownames(logCPM)[topGene],
           ylab='Log2 Expression')

## ------------------------------------------------------------------------
sum( tt$adj.P.Val < 0.05 )

## ------------------------------------------------------------------------
limmaPadj <- tt[tt$adj.P.Val <= 0.05, ] 

## ------------------------------------------------------------------------
## Plot log fold-change _versus_ -log(P-value)
## (i.e., higher number = lower p-value):
volcanoplot(fit, coef=2)

## ------------------------------------------------------------------------
library(org.Mm.eg.db)
sigGenes <- rownames(limmaPadj)
select(org.Mm.eg.db, keys = head(sigGenes), column = c("SYMBOL","GENENAME"), 
    keytype="ENSEMBL")

## ------------------------------------------------------------------------
allSigGenes <- select(org.Mm.eg.db, keys = sigGenes, column = "SYMBOL", 
                   keytype="ENSEMBL")[,2]
na.omit(allSigGenes)

## ------------------------------------------------------------------------
options(width=100)
ls("package:org.Mm.eg.db")

## ------------------------------------------------------------------------
na.omit( select(org.Mm.eg.db, keys = head(sigGenes), 
                column = c("SYMBOL", "ENSEMBL", "GO"), keytype="ENSEMBL") )

## ------------------------------------------------------------------------
library(DESeq2)
# Specify "conditions" (groups: WT and KO)
# Create object of class CountDataSet derived from eSet class
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = data.frame(groups),
                              design = ~groups)
counts(dds)[1:4, ]

## ------------------------------------------------------------------------
## Fit DESeq model to identify DE transcripts
dds <- DESeq(dds)
res <- DESeq2::results(dds)
## Remove rows with NAs
res = na.omit(res)
options(digits=2, width=100)
head(res, 4)

## ------------------------------------------------------------------------
sum(res$padj <= 0.05)

## Get the rows of "res" with significant adjusted p-values
resPadj<-res[res$padj <= 0.05 , ]

## ------------------------------------------------------------------------
library(edgeR)
# Construct DGEList object
y <- DGEList(counts=counts, group=groups) 

# Calculate library size (counts per sample)
y <- calcNormFactors(y)

# Estimate common dispersion (overall variability)
y <- estimateCommonDisp(y) 

# Estimate tagwise dispersion (per gene variability)
y <- estimateTagwiseDisp(y) 

# Compute exact test for the negative binomial distribution.
et <- exactTest(y) 

## ------------------------------------------------------------------------
topTags(et, n=4)$table

edge <- as.data.frame(topTags(et, n=nrow(counts))) 
sum(p.adjust(edge$FDR <= 0.05))

## Get the rows of "edge" with significant adjusted p-values
edgePadj <- edge[edge$FDR <= 0.05, ]

## ------------------------------------------------------------------------
library(gplots)
venn(list(edgeR=rownames(edgePadj), DESeq2=rownames(resPadj), limma=rownames(limmaPadj)))

## ------------------------------------------------------------------------
## Create DGEList object from count data
dge <- DGEList(counts=counts)

## Figure out which genes to keep
keep <- filterByExpr(dge, design)

## Apply filtering and recalculate library sizes
dge <- dge[keep, keep.lib.sizes=FALSE]

## Calculate normalisation factor (i.e., account for total reads per sample)
dge <- calcNormFactors(dge)

## ------------------------------------------------------------------------
v <- voom(dge, design, plot=TRUE)

## ------------------------------------------------------------------------
par(mfrow=c(1,3))
for(i in 1:3){
  plot(logCPM[,i], v$E[,i], xlab="LogCPM", ylab="Voom", main=colnames(logCPM)[i])
  abline(0,1)
}

## ------------------------------------------------------------------------
options(digits=4)
fit <- lmFit(v, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef=ncol(design), n=nrow(v))
head(tt)
## Get the rows of top table with significant adjusted p-values
limmaVoomPadj <- tt[tt$adj.P.Val <= 0.05, ]


## ------------------------------------------------------------------------
venn( list(limmaVoom = rownames(limmaVoomPadj), limma = rownames(limmaPadj)) )

## ------------------------------------------------------------------------
venn(list(edgeR = rownames(edgePadj), DESeq2 = rownames(resPadj),
          limmaVoom = rownames(limmaVoomPadj)))

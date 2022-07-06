## Batch of 192- Untrimmed Hippocampus STAR Aligned Data run trough featureCounts with mouse nlcRNAs GTF (Gencode) File
## Starting this on 7-1-2022
## Will import using DESEq2 from matrix
## Count matrix is from featureCounts
## Using https://lashlock.github.io/compbio/R_presentation.html as a guide



library(BiocManager)
library(tximport)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)
library(ggplot2)
library('DESeq2')
library(pheatmap)
library(dendextend)
library(readr)
library(AnnotationDbi)

## Import count matrix from Excel
Hippo_nlcRNAs <- read_excel("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/B192_Hippocampus_nlcRNAs_only_6_29_2022.xlsx")

## Delete columns with strand and chromosome information
Hippo_nlcRNAs <- Hippo_nlcRNAs[-c(2:6)]

## Create a matrix and get rid of first column (geneids)
Hippo_nlcRNAs_matrix <- as.matrix(Hippo_nlcRNAs[,-1])

## Use the GeneId olumn from Hippo_nlcRNAs as rownames for the matrix
rownames(Hippo_nlcRNAs_matrix) <- Hippo_nlcRNAs$Geneid

## Import the metadata
metadata <- read_excel("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Batch_of_192_metadata.xlsx")
## Turn the metadata into a matrix so can be imported with DESeq2
metadata <- as.matrix(metadata)


## Create a DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = Hippo_nlcRNAs_matrix,
                              colData = metadata,
                              design = ~APP)
## Look at the DESeqDataSet
dds


## Run DGE Analysis on the data
dds <- DESeq(dds)


## Create a results of the results of the DGE analysis
results <- results(dds)
## See the Top of the results
head(results)

## Get a summary of the results
summary(results)
## Sadly there are no differentially expressed nlcRNAs


## Sort the results based on the adjust p-value
results_padj_sorted <- results[order(results$padj),]
head(results_padj_sorted)

## Sort the results based on the p-value
results_pval_sorted <- results[order(results$pvalue),]
head(results_pval_sorted)


## Graph the 6 genes with the lowest p-values
par(mfrow=c(2,3))
plotCounts(dds, gene = "ENSMUSG00000107585", intgroup = "APP")
plotCounts(dds, gene = "ENSMUSG00000086395", intgroup = "APP")
plotCounts(dds, gene = "ENSMUSG00000086733", intgroup = "APP")
plotCounts(dds, gene = "ENSMUSG00000107472", intgroup = "APP")
plotCounts(dds, gene = "ENSMUSG00000105251", intgroup = "APP")
plotCounts(dds, gene = "ENSMUSG00000086331", intgroup = "APP")
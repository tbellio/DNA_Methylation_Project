## Batch of 192-Cortex Data
## This is the data that was run on Salmon 6-6-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 6/7/2022
## Trying to Use DESeq2 to get DEGs for WT vs APP
## Will import data using tximeta to create a summarized experiment First
## Using https://www.reneshbedre.com/blog/deseq2.html and https://lashlock.github.io/compbio/R_presentation.html as a guide


library(BiocManager)
library(tximport)
library(tximportData)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192")
sample_metadata <- read_excel("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Batch_of_192_metadata.xlsx")  #reads in excel file with data about each sample
sample_metadata
sample_metadata <- as.data.frame(sample_metadata)  ##turns the samples into a data frame
sample_metadata
class(sample_metadata$Sample)  #3checking class of the samples column in the samples data frame
class(sample_metadata$Sex)
sample_metadata <- transform(sample_metadata, Age = as.character(Age))  ##changes the class of the Age column in samples to character
class(sample_metadata$Age)
rownames(sample_metadata) <- sample_metadata$Sample   ##turns the row names of samples data frame into what is said in the sample column
sample_metadata$names <- sample_metadata$Sample
sample_metadata

dir <- setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Trimmed_FASTQ_Aligned_and_Quant_Salmon_6_6_2022")
files <- file.path(dir, sample_metadata$names,  "quant.sf" )
file.exists(files)

cortex_coldata <- data.frame(files, names= sample_metadata$names, APP= sample_metadata$APP, Age= sample_metadata$Age, Diet= sample_metadata$Diet, Sex= sample_metadata$Sex, stringsAsFactors = F)
cortex_coldata

## Now we have the information for tximeta, lets load it in to create a summarized experiment with all the sample metadata
## This should also find all the matching transcriptome information-- Ensembl Mus musculus release 104
cortex_se <- tximeta(cortex_coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe


suppressPackageStartupMessages(library(SummarizedExperiment))
## check the column data to see if everything loaded properly; should have all the names, APP, Age, Diet, and Sex information
colData(cortex_se)
## Check the assay names; should be 3- counts, abundance, length
assayNames(cortex_se)
## Check the row ranges
rowRanges(cortex_se)
## Check the seqinfo
seqinfo(hippo_se)

## get the database that was used to match the Salmon input quants
edb <- retrieveDb(cortex_se)
class(edb)

##adding exons to the summarized experiment (hippo_se)
cortex_se.exons <- addExons(cortex_se)

rowRanges(cortex_se.exons)

## Summarize read counts to the gene level
cortex_gene_se <- summarizeToGene(cortex_se)
rowRanges(cortex_gene_se)

## Now use DESeq2 to look for DEGs
library('DESeq2')
dds_cortex_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ APP)
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_gene_se_10filtered <- dds_cortex_gene_se[rowSums(counts(dds_cortex_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 21234

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_cortex_gene_se_10filtered$APP <- relevel(dds_cortex_gene_se_10filtered$APP, ref = "WT")

## Run DGE analysis
dge_dds_cortex_gene_se_10filtered <- DESeq(dds_cortex_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_cortex_gene_se_10filtered)

## Get gene expression table as DESeqResults
gene_expression_from_dge_dds_cortex_gene_se_10filtered <- results(dge_dds_cortex_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
gene_expression_from_dge_dds_cortex_gene_se_10filtered

## Sort the gene expression data by FDR(adjusted p-value)
gene_expression_padj_ordered_from_dge_dds_cortex_gene_se_10filtered <- gene_expression_from_dge_dds_cortex_gene_se_10filtered[order(gene_expression_from_dge_dds_cortex_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_cortex_gene_se_10filtered

##export DGE Analysis to csv File
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/Bioinformatic_Analysis_Git/DNA_Methylation_Project/Cortex-R-Analysis")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_cortex_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_padj_sorted_6_7_2022.csv")


## Get a summary of DEGs with FDR (p adj) of <0.05
summary(results(dge_dds_cortex_gene_se_10filtered, alpha = 0.05))
## Output is that of the 21234 genes, 3337 genes were upregulated in APP and 3471 were downregulated 

## Get normalized counts; DESeq2 uses median of ratios method
## Takes into account sequencing depth and RNA composition but not gene length (ok if just comparing between samples, not within)
## What it does is 1) it gets the geometric mean for each gene and makes this the 'reference' 2)calculates the ratio of each sample to the reference for each gene 3) calculates the normalization factor for each sample (so adds up all things from 2) 4)calculates the normalized count values using the median from 3
normalized_counts_dge_dds_cortex_gene_se_10filtered <- counts(dge_dds_cortex_gene_se_10filtered, normalized=TRUE)
normalized_counts_dge_dds_cortex_gene_se_10filtered

## Could export the normalized counts to an excel file for Jmp analysis if wanted but won't do that today (6-7-2022)


## use plotCounts to compare normalized counts between APP and WT for top 5 genes
par(mfrow=c(2,3))
plotCounts(dge_dds_cortex_gene_se_10filtered, gene = "ENSMUSG00000068129", intgroup = "APP")
plotCounts(dge_dds_cortex_gene_se_10filtered, gene = "ENSMUSG00000079293", intgroup = "APP")
plotCounts(dge_dds_cortex_gene_se_10filtered, gene = "ENSMUSG00000030789", intgroup = "APP")
plotCounts(dge_dds_cortex_gene_se_10filtered, gene = "ENSMUSG00000018927", intgroup = "APP")
plotCounts(dge_dds_cortex_gene_se_10filtered, gene = "ENSMUSG00000000982", intgroup = "APP")
plotCounts(dge_dds_cortex_gene_se_10filtered, gene = "ENSMUSG00000118672", intgroup = "APP")

## visualize log2 fold changes
par(mfrow=c(1,1))
plotMA(dge_dds_cortex_gene_se_10filtered, ylim=c(-2,2))

## Shrinkage of effect size (LFC)- useful for visualization and ranking
LFC_dge_dds_cortex_gene_se_10filtered <- lfcShrink(dge_dds_cortex_gene_se_10filtered, coef = "APP_APP_vs_WT", type = "apeglm")
plotMA(LFC_dge_dds_cortex_gene_se_10filtered, ylim=c(-2,2))

par(mfrow=c(1,2))
plotMA(dge_dds_cortex_gene_se_10filtered, ylim=c(-2,2))
plotMA(LFC_dge_dds_cortex_gene_se_10filtered, ylim=c(-2,2))
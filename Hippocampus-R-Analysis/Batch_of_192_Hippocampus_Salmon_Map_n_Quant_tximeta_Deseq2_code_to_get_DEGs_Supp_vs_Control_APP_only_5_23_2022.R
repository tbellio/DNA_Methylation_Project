## Batch of 192-Hippocampus Data 
## Starting this on 5/23/2022
## Trying to Use DESeq2 to get DEGs for WT vs APP at all time points
## Will import data using tximeta to create a summarized experiment First
## Using https://www.reneshbedre.com/blog/deseq2.html and https://lashlock.github.io/compbio/R_presentation.html as a guide

##trying to use tximeta to create abundances table with gene symbol

library(BiocManager)
library(tximport)
library(tximportData)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)

getwd()
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

dir <- setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_aligned_Against_Mouse_cdna_Output_5_2022")
files <- file.path(dir, sample_metadata$names,  "quant.sf" )
file.exists(files)

## Create dataframe that has the location of the quant files and all sample metadata
hippo_coldata <- data.frame(files, names= sample_metadata$names, APP= sample_metadata$APP, Age= sample_metadata$Age, Diet= sample_metadata$Diet, Sex= sample_metadata$Sex, stringsAsFactors = F)
hippo_coldata

## Subset data to only include APP animals
hippo_app_only_coldata <- hippo_coldata[hippo_coldata$APP== c("APP"), ]
## Subset data to only include animals 6-12 months old
hippo_older_than_3_app_only_coldata <- hippo_app_only_coldata[hippo_app_only_coldata$Age %in% c(6, 9, 12), ]


## Now we have the information for tximeta, lets load it in to create a summarized experiment with all the sample metadata
## This should also find all the matching transcriptome information-- Ensembl Mus musculus release 104
hippo_se <- tximeta(hippo_older_than_3_app_only_coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe


suppressPackageStartupMessages(library(SummarizedExperiment))
## check the column data to see if everything loaded properly; should have all the names, APP, Age, Diet, and Sex information
colData(hippo_se)
## Check the assay names; should be 3- counts, abundance, length
assayNames(hippo_se)
## Check the row ranges
rowRanges(hippo_se)
## Check the seqinfo
seqinfo(hippo_se)

## get the database that was used to match the Salmon input quants
edb <- retrieveDb(hippo_se)
class(edb)

##adding exons to the summarized experiment (hippo_se)
hippo_se.exons <- addExons(hippo_se)

rowRanges(hippo_se.exons)

## Summarize read counts to the gene level
hippo_gene_se <- summarizeToGene(hippo_se)
rowRanges(hippo_gene_se)


## Now use DESeq2 to look for DEGs
library('DESeq2')
dds_hippo_gene_se <- DESeqDataSet(hippo_gene_se, design = ~ Diet)
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor

## Filter so that only genes that have more than 10 reads are inculuded
dds_hippo_gene_se_10filtered <- dds_hippo_gene_se[rowSums(counts(dds_hippo_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 20764

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_hippo_gene_se_10filtered$Diet <- relevel(dds_hippo_gene_se_10filtered$Diet, ref = "Control")

## Run DGE analysis
dge_dds_hippo_gene_se_10filtered <- DESeq(dds_hippo_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just supp vs control diet
resultsNames(dge_dds_hippo_gene_se_10filtered)

## Get gene expression table as DESeqResults
gene_expression_from_dge_dds_hippo_gene_se_10filtered <- results(dge_dds_hippo_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so supp vs control for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
gene_expression_from_dge_dds_hippo_gene_se_10filtered

## Sort the gene expression data by FDR(adjusted p-value)
gene_expression_padj_ordered_from_dge_dds_hippo_gene_se_10filtered <- gene_expression_from_dge_dds_hippo_gene_se_10filtered[order(gene_expression_from_dge_dds_hippo_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_hippo_gene_se_10filtered

##export DGE Analysis to csv File
getwd()
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_hippo_gene_se_10filtered), file = "Batch_of_192_Hippocampus_DESeq2_SuppvsCont_greater_than_3mo_DGE_Analysis_padj_sorted_5_23_2022.csv")


## Get a summary of DEGs with FDR (p adj) of <0.05
summary(results(dge_dds_hippo_gene_se_10filtered, alpha = 0.05))
## Output is that of the 19408 genes, 1 gene was upregulated in Supplemented and 4 were downregulated 

##Make the results into a dataframe
results <- as.data.frame(results(dge_dds_hippo_gene_se_10filtered))
##order the results dataframe by pvalue instead of FDR
results <- results[order(results$pvalue),]
##export results dataframe to Execel
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(results, file = "Batch_of_192_Hippocampus_DESeq2_SuppvsCont_greater_than_3mo_DGE_Analysis_pvalue_sorted_5_25_2022.csv")

## Get normalized counts; DESeq2 uses median of ratios method
## Takes into account sequencing depth and RNA composition but not gene length (ok if just comparing between samples, not within)
## What it does is 1) it gets the geometric mean for each gene and makes this the 'reference' 2)calculates the ratio of each sample to the reference for each gene 3) calculates the normalization factor for each sample (so adds up all things from 2) 4)calculates the normalized count values using the median from 3
normalized_counts_dge_dds_hippo_gene_se_10filtered <- counts(dge_dds_hippo_gene_se_10filtered, normalized=TRUE)
normalized_counts_dge_dds_hippo_gene_se_10filtered


## Visualize the 5 genes that are DE
par(mfrow=c(2,3))
plotCounts(dge_dds_hippo_gene_se_10filtered, gene = "ENSMUSG00000079224", intgroup = "Diet")
plotCounts(dge_dds_hippo_gene_se_10filtered, gene = "ENSMUSG00000102070", intgroup = "Diet")
plotCounts(dge_dds_hippo_gene_se_10filtered, gene = "ENSMUSG00000040794", intgroup = "Diet")
plotCounts(dge_dds_hippo_gene_se_10filtered, gene = "ENSMUSG00000049382", intgroup = "Diet")
plotCounts(dge_dds_hippo_gene_se_10filtered, gene = "ENSMUSG00000052837", intgroup = "Diet")

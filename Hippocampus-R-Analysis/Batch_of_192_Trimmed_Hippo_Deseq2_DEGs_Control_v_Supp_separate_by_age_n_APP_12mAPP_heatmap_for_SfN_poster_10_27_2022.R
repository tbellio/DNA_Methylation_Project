## Batch of 192- Trimmed Hippocampus Data 
## This is the data that was run on Salmon 6-6-2022 and 6-7-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 10/27/2022
## Will import data using tximeta to create a summarized experiment First
## Using https://www.reneshbedre.com/blog/deseq2.html and https://lashlock.github.io/compbio/R_presentation.html as a guide
## Trying to Use DESeq2 to get DEGs for Control vs Supp diet at all time points but separated into just APP and just WT (so a comparison of Control vs Supplemented at 3, 6, 9, and 12 months old)
## Then taking the DEGs from 12M APP animals to make heatmap to see if we can separate out diet by clustering



library(BiocManager)
library(tximport)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)
library(ggplot2)
library(pheatmap)

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192")
sample_metadata <- read_excel("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Batch_of_192_metadata_with_RIN.xlsx")  #reads in excel file with data about each sample
sample_metadata
sample_metadata <- as.data.frame(sample_metadata)  ##turns the samples into a data frame
sample_metadata
class(sample_metadata$Sample)  #checking class of the samples column in the samples data frame
class(sample_metadata$Sex)
rownames(sample_metadata) <- sample_metadata$Sample   ##turns the row names of samples data frame into what is said in the sample column
sample_metadata$names <- sample_metadata$Sample
sample_metadata <- transform(sample_metadata, Age = as.factor(Age))  ##changes the class of the Age column in samples to character
class(sample_metadata$Age)
sample_metadata <- transform(sample_metadata, APP = as.factor(APP))  ##changes the class of the Age column in samples to character
class(sample_metadata$APP)
sample_metadata <- transform(sample_metadata, Diet = as.factor(Diet))  ##changes the class of the Age column in samples to character
class(sample_metadata$Diet)
sample_metadata <- transform(sample_metadata, Sex = as.factor(Sex))  ##changes the class of the Age column in samples to character
class(sample_metadata$Sex)
class(sample_metadata$Hippocampus.RIN)
sample_metadata <- transform(sample_metadata, Hippocampus.RIN = as.factor(Hippocampus.RIN))
class(sample_metadata$Hippocampus.RIN)
class(sample_metadata$Cortex.RIN)
sample_metadata <- transform(sample_metadata, Cortex.RIN=as.factor(Cortex.RIN))
class(sample_metadata$Cortex.RIN)
sample_metadata

dir <- setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Trimmed_FASTQ_Aligned_and_Quant_Salmon_6_6_2022")
files <- file.path(dir, sample_metadata$names,  "quant.sf" )
file.exists(files)

## Create dataframe that has the location of the quant files and all sample metadata
hippo_coldata <- data.frame(files, names= sample_metadata$names, APP= sample_metadata$APP, Age= sample_metadata$Age, Diet= sample_metadata$Diet, Sex= sample_metadata$Sex, RIN= sample_metadata$Hippocampus.RIN, stringsAsFactors = F)
hippo_coldata

## Subset data to only include 3, 6, 9, and 12 month animals
## Make sure only 24 rows for each
hippo_three_only_coldata <- hippo_coldata[hippo_coldata$Age== c("3"), ]
hippo_six_only_coldata <- hippo_coldata[hippo_coldata$Age== c("6"),]
hippo_nine_only_coldata <- hippo_coldata[hippo_coldata$Age== c("9"),]
hippo_twelve_only_coldata <- hippo_coldata[hippo_coldata$Age== c("12"),]

##Subset for APP Only and WT Only
APP_hippo_coldata <- hippo_coldata[hippo_coldata$APP== c("APP"),]
APP_hippo_three_only_coldata <- hippo_three_only_coldata[hippo_three_only_coldata$APP== c("APP"),]
APP_hippo_six_only_coldata <- hippo_six_only_coldata[hippo_six_only_coldata$APP== c("APP"),]
APP_hippo_nine_only_coldata <- hippo_nine_only_coldata[hippo_nine_only_coldata$APP== c("APP"),]
APP_hippo_twelve_only_coldata <- hippo_twelve_only_coldata[hippo_twelve_only_coldata$APP== c("APP"),]

WT_hippo_coldata <- hippo_coldata[hippo_coldata$APP== c("WT"),]
WT_hippo_three_only_coldata <- hippo_three_only_coldata[hippo_three_only_coldata$APP== c("WT"),]
WT_hippo_six_only_coldata <- hippo_six_only_coldata[hippo_six_only_coldata$APP== c("WT"),]
WT_hippo_nine_only_coldata <- hippo_nine_only_coldata[hippo_nine_only_coldata$APP== c("WT"),]
WT_hippo_twelve_only_coldata <- hippo_twelve_only_coldata[hippo_twelve_only_coldata$APP== c("WT"),]


## Now we have the information for tximeta, lets load it in to create a summarized experiment with all the sample metadata
## This should also find all the matching transcriptome information-- Ensembl Mus musculus release 104
APP_hippo_se <- tximeta(APP_hippo_coldata)
APP_hippo__three_se <- tximeta(APP_hippo_three_only_coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe
APP_hippo_six_se <- tximeta(APP_hippo_six_only_coldata)
APP_hippo_nine_se <- tximeta(APP_hippo_nine_only_coldata)
APP_hippo_twelve_se <- tximeta(APP_hippo_twelve_only_coldata)

WT_hippo_se <- tximeta(WT_hippo_coldata)
WT_hippo__three_se <- tximeta(WT_hippo_three_only_coldata)
WT_hippo_six_se <- tximeta(WT_hippo_six_only_coldata)
WT_hippo_nine_se <- tximeta(WT_hippo_nine_only_coldata)
WT_hippo_twelve_se <- tximeta(WT_hippo_twelve_only_coldata)

hippo_se <- tximeta(hippo_coldata)


suppressPackageStartupMessages(library(SummarizedExperiment))
##adding exons to the summarized experiment (hippo__three_se, etc.)
APP_hippo_se.exons <- addExons(APP_hippo_se)
APP_hippo_three_se.exons <- addExons(APP_hippo__three_se)
APP_hippo_six_se.exons <- addExons(APP_hippo_six_se)
APP_hippo_nine_se.exons <- addExons(APP_hippo_nine_se)
APP_hippo_twelve_se.exons <- addExons(APP_hippo_twelve_se)

WT_hippo_se.exons <- addExons(WT_hippo_se)
WT_hippo_three_se.exons <- addExons(WT_hippo__three_se)
WT_hippo_six_se.exons <- addExons(WT_hippo_six_se)
WT_hippo_nine_se.exons <- addExons(WT_hippo_nine_se)
WT_hippo_twelve_se.exons <- addExons(WT_hippo_twelve_se)

hippo_se.exons <- addExons(hippo_se)

##Check rowranges--shouldnt have to make sure it works for all four as same code
rowRanges(APP_hippo_three_se.exons)

## Summarize read counts to the gene level
APP_hippo_gene_se <- summarizeToGene(APP_hippo_se.exons)
APP_hippo_three_gene_se <- summarizeToGene(APP_hippo_three_se.exons)
APP_hippo_six_gene_se <- summarizeToGene(APP_hippo_six_se.exons)
APP_hippo_nine_gene_se <- summarizeToGene(APP_hippo_nine_se.exons)
APP_hippo_twelve_gene_se <- summarizeToGene(APP_hippo_twelve_se.exons)

WT_hippo_gene_se <- summarizeToGene(WT_hippo_se.exons)
WT_hippo_three_gene_se <- summarizeToGene(WT_hippo_three_se.exons)
WT_hippo_six_gene_se <- summarizeToGene(WT_hippo_six_se.exons)
WT_hippo_nine_gene_se <- summarizeToGene(WT_hippo_nine_se.exons)
WT_hippo_twelve_gene_se <- summarizeToGene(WT_hippo_twelve_se.exons)

hippo_gene_se <- summarizeToGene(hippo_se.exons)

rowRanges(APP_hippo_three_gene_se)


## Now use DESeq2 to look for DEGs at each time point
library('DESeq2')
dds_APP_hippo_gene_se <- DESeqDataSet(APP_hippo_gene_se, design = ~ Age + Sex + RIN + Diet)
dds_APP_hippo_three_gene_se <- DESeqDataSet(APP_hippo_three_gene_se, design = ~ Sex + RIN + Diet)
dds_APP_hippo_six_gene_se <- DESeqDataSet(APP_hippo_six_gene_se, design = ~ Sex + RIN + Diet)
dds_APP_hippo_nine_gene_se <- DESeqDataSet(APP_hippo_nine_gene_se, design = ~ Sex + RIN + Diet)
dds_APP_hippo_twelve_gene_se <- DESeqDataSet(APP_hippo_twelve_gene_se, design = ~  Sex + RIN + Diet)

dds_WT_hippo_gene_se <- DESeqDataSet(WT_hippo_gene_se, design = ~ Age + Sex + RIN + Diet)
dds_WT_hippo_three_gene_se <- DESeqDataSet(WT_hippo_three_gene_se, design = ~ Sex + RIN + Diet)
dds_WT_hippo_six_gene_se <- DESeqDataSet(WT_hippo_six_gene_se, design = ~ Sex + Diet)
dds_WT_hippo_nine_gene_se <- DESeqDataSet(WT_hippo_nine_gene_se, design = ~ Sex + RIN + Diet)
dds_WT_hippo_twelve_gene_se <- DESeqDataSet(WT_hippo_twelve_gene_se, design = ~ Sex + RIN + Diet)


dds_hippo_gene_se <- DESeqDataSet(hippo_gene_se, design = ~ Age + Sex + APP + RIN + Diet)
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; if not it will change it to a factor (this code has it as a character so it knows to change it)

## Filter so that only genes that have more than 10 reads are inculuded
dds_APP_hippo_gene_se_10filtered <-  dds_APP_hippo_gene_se[rowSums(counts(dds_APP_hippo_gene_se)) >= 10,]
dds_APP_hippo_three_gene_se_10filtered <- dds_APP_hippo_three_gene_se[rowSums(counts(dds_APP_hippo_three_gene_se)) >= 10,]
dds_APP_hippo_six_gene_se_10filtered <- dds_APP_hippo_six_gene_se[rowSums(counts(dds_APP_hippo_six_gene_se)) >= 10,]
dds_APP_hippo_nine_gene_se_10filtered <- dds_APP_hippo_nine_gene_se[rowSums(counts(dds_APP_hippo_nine_gene_se)) >= 10,]
dds_APP_hippo_twelve_gene_se_10filtered <- dds_APP_hippo_twelve_gene_se[rowSums(counts(dds_APP_hippo_twelve_gene_se)) >= 10,]

dds_WT_hippo_gene_se_10filtered <- dds_WT_hippo_gene_se[rowSums(counts(dds_WT_hippo_gene_se)) >= 10,]
dds_WT_hippo_three_gene_se_10filtered <- dds_WT_hippo_three_gene_se[rowSums(counts(dds_WT_hippo_three_gene_se)) >= 10,]
dds_WT_hippo_six_gene_se_10filtered <- dds_WT_hippo_six_gene_se[rowSums(counts(dds_WT_hippo_six_gene_se)) >= 10,]
dds_WT_hippo_nine_gene_se_10filtered <- dds_WT_hippo_nine_gene_se[rowSums(counts(dds_WT_hippo_nine_gene_se)) >= 10,]
dds_WT_hippo_twelve_gene_se_10filtered <- dds_WT_hippo_twelve_gene_se[rowSums(counts(dds_WT_hippo_twelve_gene_se)) >= 10,]

dds_hippo_gene_se_10filtered <- dds_hippo_gene_se[rowSums(counts(dds_hippo_gene_se)) >= 10,]
## this filtered out ~18000 genes; went from 35682 to 18000-19000

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the contorl diet mice
dds_APP_hippo_gene_se_10filtered$Diet <- relevel(dds_APP_hippo_gene_se_10filtered$Diet, ref = "Control")
dds_APP_hippo_three_gene_se_10filtered$Diet <- relevel(dds_APP_hippo_three_gene_se_10filtered$Diet, ref = "Control")
dds_APP_hippo_six_gene_se_10filtered$Diet <- relevel(dds_APP_hippo_six_gene_se_10filtered$Diet, ref = "Control")
dds_APP_hippo_nine_gene_se_10filtered$Diet <- relevel(dds_APP_hippo_nine_gene_se_10filtered$Diet, ref = "Control")
dds_APP_hippo_twelve_gene_se_10filtered$Diet <- relevel(dds_APP_hippo_twelve_gene_se_10filtered$Diet, ref = "Control")

dds_WT_hippo_gene_se_10filtered$Diet <- relevel(dds_WT_hippo_gene_se_10filtered$Diet, ref = "Control")
dds_WT_hippo_three_gene_se_10filtered$Diet <- relevel(dds_WT_hippo_three_gene_se_10filtered$Diet, ref = "Control")
dds_WT_hippo_six_gene_se_10filtered$Diet <- relevel(dds_WT_hippo_six_gene_se_10filtered$Diet, ref = "Control")
dds_WT_hippo_nine_gene_se_10filtered$Diet <- relevel(dds_WT_hippo_nine_gene_se_10filtered$Diet, ref = "Control")
dds_WT_hippo_twelve_gene_se_10filtered$Diet <- relevel(dds_WT_hippo_twelve_gene_se_10filtered$Diet, ref = "Control")

dds_hippo_gene_se_10filtered$Diet <- relevel(dds_hippo_gene_se_10filtered$Diet, ref = "Control")

## Run DGE analysis
dge_dds_APP_hippo_gene_se_10filtered <- DESeq(dds_APP_hippo_gene_se_10filtered)
dge_dds_APP_hippo_three_gene_se_10filtered <- DESeq(dds_APP_hippo_three_gene_se_10filtered)
dge_dds_APP_hippo_six_gene_se_10filtered <- DESeq(dds_APP_hippo_six_gene_se_10filtered)
dge_dds_APP_hippo_nine_gene_se_10filtered <- DESeq(dds_APP_hippo_nine_gene_se_10filtered)
dge_dds_APP_hippo_twelve_gene_se_10filtered <- DESeq(dds_APP_hippo_twelve_gene_se_10filtered)

dge_dds_WT_hippo_gene_se_10filtered <- DESeq(dds_WT_hippo_gene_se_10filtered)
dge_dds_WT_hippo_three_gene_se_10filtered <- DESeq(dds_WT_hippo_three_gene_se_10filtered)
dge_dds_WT_hippo_six_gene_se_10filtered <- DESeq(dds_WT_hippo_six_gene_se_10filtered)
dge_dds_WT_hippo_nine_gene_se_10filtered <- DESeq(dds_WT_hippo_nine_gene_se_10filtered)
dge_dds_WT_hippo_twelve_gene_se_10filtered <- DESeq(dds_WT_hippo_twelve_gene_se_10filtered)

dge_dds_hippo_gene_se_10filtered <- DESeq(dds_hippo_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just Control vs Supp
resultsNames(dge_dds_APP_hippo_gene_se_10filtered)
resultsNames(dge_dds_APP_hippo_three_gene_se_10filtered)
resultsNames(dge_dds_APP_hippo_six_gene_se_10filtered)
resultsNames(dge_dds_APP_hippo_nine_gene_se_10filtered)
resultsNames(dge_dds_APP_hippo_twelve_gene_se_10filtered)

resultsNames(dge_dds_WT_hippo_gene_se_10filtered)
resultsNames(dge_dds_WT_hippo_three_gene_se_10filtered)
resultsNames(dge_dds_WT_hippo_six_gene_se_10filtered)
resultsNames(dge_dds_WT_hippo_nine_gene_se_10filtered)
resultsNames(dge_dds_WT_hippo_twelve_gene_se_10filtered)

resultsNames(dge_dds_hippo_gene_se_10filtered)

## Get gene expression table as DESeqResults
gene_expression_from_dge_dds_APP_hippo_gene_se_10filtered <- results(dge_dds_APP_hippo_gene_se_10filtered)
gene_expression_from_dge_dds_APP_hippo_three_gene_se_10filtered <- results(dge_dds_APP_hippo_three_gene_se_10filtered)
gene_expression_from_dge_dds_APP_hippo_six_gene_se_10filtered <- results(dge_dds_APP_hippo_six_gene_se_10filtered)
gene_expression_from_dge_dds_APP_hippo_nine_gene_se_10filtered <- results(dge_dds_APP_hippo_nine_gene_se_10filtered)
gene_expression_from_dge_dds_APP_hippo_twelve_gene_se_10filtered <- results(dge_dds_APP_hippo_twelve_gene_se_10filtered)

gene_expression_from_dge_dds_WT_hippo_gene_se_10filtered <- results(dge_dds_WT_hippo_gene_se_10filtered)
gene_expression_from_dge_dds_WT_hippo_three_gene_se_10filtered <- results(dge_dds_WT_hippo_three_gene_se_10filtered)
gene_expression_from_dge_dds_WT_hippo_six_gene_se_10filtered <- results(dge_dds_WT_hippo_six_gene_se_10filtered)
gene_expression_from_dge_dds_WT_hippo_nine_gene_se_10filtered <- results(dge_dds_WT_hippo_nine_gene_se_10filtered)
gene_expression_from_dge_dds_WT_hippo_twelve_gene_se_10filtered <- results(dge_dds_WT_hippo_twelve_gene_se_10filtered)

gene_expression_from_dge_dds_hippo_gene_se_10filtered <- results(dge_dds_hippo_gene_se_10filtered)

## See what the gene expression table says
## Will tell you what it was comparing (so APP vs WT for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
gene_expression_from_dge_dds_APP_hippo_gene_se_10filtered
gene_expression_from_dge_dds_APP_hippo_three_gene_se_10filtered
gene_expression_from_dge_dds_APP_hippo_six_gene_se_10filtered
gene_expression_from_dge_dds_APP_hippo_nine_gene_se_10filtered
gene_expression_from_dge_dds_APP_hippo_twelve_gene_se_10filtered

## Sort the gene expression data by FDR(adjusted p-value)
gene_expression_padj_ordered_from_dge_dds_APP_hippo_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_APP_hippo_three_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_three_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_three_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_APP_hippo_six_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_six_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_six_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_APP_hippo_nine_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_nine_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_APP_hippo_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_twelve_gene_se_10filtered$padj),]

gene_expression_padj_ordered_from_dge_dds_WT_hippo_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_WT_hippo_three_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_three_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_three_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_WT_hippo_six_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_six_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_six_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_WT_hippo_nine_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_nine_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_WT_hippo_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_twelve_gene_se_10filtered$padj),]

gene_expression_padj_ordered_from_dge_dds_hippo_gene_se_10filtered <- gene_expression_from_dge_dds_hippo_gene_se_10filtered[order(gene_expression_from_dge_dds_hippo_gene_se_10filtered$padj),]

gene_expression_padj_ordered_from_dge_dds_APP_hippo_six_gene_se_10filtered 


## Sort by p-value
gene_expression_pval_ordered_from_dge_dds_APP_hippo_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_APP_hippo_three_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_three_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_three_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_APP_hippo_six_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_six_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_six_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_APP_hippo_nine_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_nine_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_APP_hippo_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_APP_hippo_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_hippo_twelve_gene_se_10filtered$pvalue),]

gene_expression_pval_ordered_from_dge_dds_WT_hippo_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_WT_hippo_three_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_three_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_three_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_WT_hippo_six_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_six_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_six_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_WT_hippo_nine_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_nine_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_WT_hippo_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_WT_hippo_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_hippo_twelve_gene_se_10filtered$pvalue),]

gene_expression_pval_ordered_from_dge_dds_hippo_gene_se_10filtered <- gene_expression_from_dge_dds_hippo_gene_se_10filtered[order(gene_expression_from_dge_dds_hippo_gene_se_10filtered$pvalue),]


## All Together
summary(results(dge_dds_hippo_gene_se_10filtered, alpha = 0.05))
## 5 genes were upregulated and 20 genes were downrefgulated in supplemented mice
summary(results(dge_dds_hippo_gene_se_10filtered))
## 7 genes were upregulated and 22 genes were downregulated in supplemented mice

## Get a summary of DEGs with FDR (padj) of <0.05
summary(results(dge_dds_APP_hippo_gene_se_10filtered, alpha = 0.05))
## 14 genes were upregulated, 22 genes were upregulated in supp APP mice
summary(results(dge_dds_WT_hippo_gene_se_10filtered, alpha = 0.05))
## 9 gene was upregulated and 7 gene downregulated in supp WT mice
## Get a summary of DEGs with FDR (padj) of <0.05
summary(results(dge_dds_APP_hippo_gene_se_10filtered))
## 14 genes were upregulated, 27 genes were upregulated in supp APP mice
summary(results(dge_dds_WT_hippo_gene_se_10filtered))
## 10 gene was upregulated and 10 genes downregulated in supp WT mice




## 3 Months

## Get a summary of DEGs with FDR (p adj) of <0.05
summary(results(dge_dds_APP_hippo_three_gene_se_10filtered, alpha = 0.05))
## Output is that of the 18837 genes,1 genes were upregulated and 12 were downregulated in supplemented 3 motnh old APP mice
summary(results(dge_dds_WT_hippo_three_gene_se_10filtered, alpha = 0.05))
## Output is that of the 18837 genes,296 genes were upregulated and 188 was downregulated in supplemented 3 motnh old WT mice
## Get a summary of DEGs with FDR (p adj) of <0.05
summary(results(dge_dds_APP_hippo_three_gene_se_10filtered))
## Output is that of the 18837 genes,2 genes were upregulated and 14 were downregulated in supplemented 3 motnh old APP mice
summary(results(dge_dds_WT_hippo_three_gene_se_10filtered))
## Output is that of the 18837 genes,397 genes were upregulated and 318 was downregulated in supplemented 3 motnh old WT mice


## 6 Months

summary(results(dge_dds_APP_hippo_six_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19242 genes, 1 genes were upregulated and 2 genes were downregulated in supp mice 6mo APP mice
summary(results(dge_dds_WT_hippo_six_gene_se_10filtered, alpha = 0.05))
## Output is that of the 19242 genes,0 genes were upregulated and 1 were downregulated in supplemented 6 motnh old WT mice
summary(results(dge_dds_APP_hippo_six_gene_se_10filtered))
##Output is that of the 19242 genes, 1 genes were upregulated and 2 genes were downregulated in supp mice 6mo APP mice
summary(results(dge_dds_WT_hippo_six_gene_se_10filtered))
## Output is that of the 19242 genes,0 genes were upregulated and 1 were downregulated in supplemented 6 motnh old WT mice



## 9 Months

summary(results(dge_dds_APP_hippo_nine_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19225 genes, 2 genes were upregulated and 61 genes were downregulated in supplemented 9mo APP mice
summary(results(dge_dds_WT_hippo_nine_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19225 genes, no genes were differentially expressed in supplemented 9mo WT mice
summary(results(dge_dds_APP_hippo_nine_gene_se_10filtered))
##Output is that of the 19225 genes, 2 genes were upregulated and 7 genes were downregulated in supplemented 9mo APP mice
summary(results(dge_dds_WT_hippo_nine_gene_se_10filtered))
##Output is that of the 19225 genes, no genes were differentially expressed in supplemented 9mo WT mice



## 12 Months

summary(results(dge_dds_APP_hippo_twelve_gene_se_10filtered, alpha = 0.05))
##Output is that of the 18741 genes, 0 genes was upregulated and 4 gene was downregulated in supp APP mice at 12 months
summary(results(dge_dds_WT_hippo_twelve_gene_se_10filtered, alpha = 0.05))
##Output is that of the 18741 genes, 1 genes were upregulated and 0 were downregulated in supp WT mice at 12 months
summary(results(dge_dds_APP_hippo_twelve_gene_se_10filtered))
##Output is that of the 18741 genes, 0 genes was upregulated and 5 gene was downregulated in supp APP mice at 12 months
summary(results(dge_dds_WT_hippo_twelve_gene_se_10filtered))
##Output is that of the 18741 genes, 1 genes were upregulated and 3 were downregulated in supp WT mice at 12 months



##export DGE Analysis to csv File
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_gene_se_10filtered), file = 'Trimmed_Hippo_DESeq2_Control_v_Supp_APP_DGE_Analysis_Age_sex_RIN_controlled_padj_sorted_10_27_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_three_gene_se_10filtered), file = 'Trimmed_Hippo_DESeq2_Control_v_Supp_3mo_APP_DGE_Analysis_sex_RIN_controlled_padj_sorted_10_27_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_six_gene_se_10filtered), file = "Trimmed_Hippo_DESeq2_Control_v_Supp_6mo_APP_DGE_Analysis_sex_RIN_controlled_padj_sorted_10_27_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_nine_gene_se_10filtered), file = "Trimmed_Hippo_DESeq2_Control_v_Supp_9mo_APP_DGE_Analysis_sex_RIN_controlled_padj_sorted_10_27_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_twelve_gene_se_10filtered), file = "Trimmed_Hippo_DESeq2_Control_v_Supp_12mo_APP_DGE_Analysis_sex_RIN_controlled_padj_sorted_10_27_2022.csv")

write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_gene_se_10filtered), file = 'Trimmed_Hippo_DESeq2_Control_v_Supp_WT_DGE_Analysis_age_sex_RIN_controlled_padj_sorted_10_27_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_three_gene_se_10filtered), file = 'Trimmed_Hippo_DESeq2_Control_v_Supp_3mo_WT_DGE_Analysis_sex_RIN_controlled_padj_sorted_10_27_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_six_gene_se_10filtered), file = "Trimmed_Hippo_DESeq2_Control_v_Supp_6mo_WT_DGE_Analysis_sex_RIN_controlled_padj_sorted_10_27_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_nine_gene_se_10filtered), file = "Trimmed_Hippo_DESeq2_Control_v_Supp_9mo_WT_DGE_Analysis_sex_RIN_controlled_padj_sorted_10_27_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_twelve_gene_se_10filtered), file = "Trimmed_Hippo_DESeq2_Control_v_Supp_12mo_WT_DGE_Analysis_sex_RIN_controlled_padj_sorted_10_27_2022.csv")

write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_hippo_gene_se_10filtered), file = "Trimmed_Hippo_DESeq2_Control_v_Supp_All_animals_DGE_Analysis_sex_age_app_RIN_controlled_padj_sorted_10_27_2022.csv")

write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_hippo_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_APP_DGE_Analysis_age_sex_controlled_pval_sorted_7_23_2022.csv')
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_hippo_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_3mo_APP_DGE_Analysis_sex_controlled_pval_sorted_7_23_2022.csv')
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_hippo_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_6mo_APP_DGE_Analysis_sex_controlled_pval_sorted_7_23_2022.csv")
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_hippo_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_9mo_APP_DGE_Analysis_sex_controlled_pval_sorted_7_23_2022.csv")
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_hippo_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_12mo_APP_DGE_Analysis_sex_controlled_pval_sorted_7_23_2022.csv")

write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_hippo_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_WT_DGE_Analysis_age_sex_controlled_pval_sorted_7_23_2022.csv')
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_hippo_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_3mo_WT_DGE_Analysis_sex_controlled_pval_sorted_7_23_2022.csv')
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_hippo_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_6mo_WT_DGE_Analysis_sex_controlled_pval_sorted_7_23_2022.csv")
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_hippo_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_9mo_WT_DGE_Analysis_sex_controlled_pval_sorted_7_23_2022.csv")
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_hippo_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_12mo_WT_DGE_Analysis_sex_controlled_pval_sorted_7_23_2022.csv")

write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_hippo_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_All_animals_DGE_Analysis_sex_age_app_controlled_pval_sorted_7_23_2022.csv")




## Trying to make heatmap of DEGs from all APP
gene_expression_padj_ordered_from_dge_dds_APP_hippo_twelve_gene_se_10filtered
ALL_APP_Diet_DEGs <- gene_expression_padj_ordered_from_dge_dds_APP_hippo_gene_se_10filtered
## Only 41 DEGs so subset out first 41 rows
top41 <- row.names(ALL_APP_Diet_DEGs[1:41,])
top41


## Make matrix of all APP animal abundances
filtered_ALL_APP_gene_abundances <- as.data.frame(dge_dds_APP_hippo_gene_se_10filtered@assays@data@listData[["abundance"]])
filtered_ALL_APP_gene_abundances
filtered_ALL_APP_gene_abundances_matrix <- as.matrix(filtered_ALL_APP_gene_abundances)
filtered_ALL_APP_gene_abundances_matrix


## Making the rownames gene symbols instead of ENSEMBL names
gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, row.names(filtered_ALL_APP_gene_abundances_matrix), "SYMBOL", keytype = "ENSEMBL")
##this creates gene_symbols matrix that has 2 columns (ensembl gene name and gene symbol); however there are more genes in this annotation than we have in our filtered_gene_abundances matrix
length(unique(gene_symbols$ENSEMBL))  #19541 unique emsembl gene names
length(unique(gene_symbols$SYMBOL))  ##18210 unique gene symbols
## so we need to adjust gene_symbols so we dont have any overlapping genes
gene_symbols[duplicated(gene_symbols$ENSEMBL),]  #gives you the duplicated ensembl gene names
##use dplyr to select just the genes that are not duplicated and call it nonduplicated_gene_symbols
library(dplyr)
nonduplicated_gene_symbols <- distinct(gene_symbols, ENSEMBL, .keep_all = T)
##now make sure that the nonduplicated_gene_symbols is the same length as your filtered_gene_abundances
nrow(filtered_ALL_APP_gene_abundances_matrix)
nrow(nonduplicated_gene_symbols)
## or 
nrow(filtered_ALL_APP_gene_abundances_matrix) ==  nrow(nonduplicated_gene_symbols) ##make sure this is TRUE

##now you need to make sure that the row names of your abundance data set matches the data in your annotation
test <- row.names(filtered_ALL_APP_gene_abundances_matrix) %in% nonduplicated_gene_symbols$ENSEMBL
sum(test) #make sure that what this returns matches your number of rows in both your abundance data set and your annotation (gene symbols)

##now you can set the row.names of your filtered_gene_abundancs to the gene symbols from nonduplicated_gene_symbols since you know they match
row.names(filtered_ALL_APP_gene_abundances_matrix) <- nonduplicated_gene_symbols$SYMBOL
## Check to see if row names have been renamed
filtered_ALL_APP_gene_abundances_matrix ## They have

## Switch top 41 list from ENSEMBL names to Gene symbols
top41 <- as.data.frame(top41)
top41
top41_gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, top41$top41, "SYMBOL", keytype = "ENSEMBL")
top41_gene_symbols
## Need to make sure all values in symbol column have  a value and if it is NA then switch to the ENSEMBL name
top41_gene_symbols$SYMBOL[is.na(top41_gene_symbols$SYMBOL)] <- top41_gene_symbols$ENSEMBL[is.na(top41_gene_symbols$SYMBOL)]
top41_gene_symbols
row.names(top41_gene_symbols) <- top41_gene_symbols$SYMBOL
top41_gene_symbols
top41 <- row.names(top41_gene_symbols)
top41

## now need to subset the filtered_cortex_gene_abundances_matrix so that we only have the top 50 genes remaining
top41_filtered_hippo <- filtered_ALL_APP_gene_abundances_matrix[rownames(filtered_ALL_APP_gene_abundances_matrix) %in% top41,]
top41_filtered_hippo
## this actually drops it down to 25 DEGs as it doesnt have a gene name but lets move on and not worry about it

## make heatmap
pheatmap(top41_filtered_hippo)
## this looks awful because things aren't normalized so lets normalize
## make function to z-score
cal_z_score <- function(x) {
  (x- mean(x)) / sd(x)
}

## apply z-score function to top 100 gene abundances
top41_filtered_hippo_gene_abundances_matrix_z_score <- t(apply(top41_filtered_hippo,1, cal_z_score))

## Make heatmap of top 100 DEGs APP vs WT (all ages) in z-score
pheatmap(top41_filtered_hippo_gene_abundances_matrix_z_score)

## Make annotation information with metadata
my_sample_col <- data.frame(APP = sample_metadata$APP, Diet= sample_metadata$Diet, Age= sample_metadata$Age)
rownames(my_sample_col) <- sample_metadata$Sample
my_sample_col

## Make top 41 Diet DEGs heatmap
pheatmap(top41_filtered_hippo_gene_abundances_matrix_z_score, annotation_col = my_sample_col)
annot_colors_gradient= list(APP=c(APP="red1", WT= "royalblue1"), Diet= c(Control= "gold", Supplemented= "goldenrod3"), Age= c('3'='gray99', '6'='gray75', '9'= 'gray50', '12'= 'gray25'))
palettelength2 <- 50
myColor_ex2 <- colorRampPalette(c("#762A83", "#9970AB", "#C2A5CF", '#E7D4E8', "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B"))(palettelength2)
myColor_ex2
hippo_top50_heatmap_overflow <- pheatmap(top41_filtered_hippo_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors = annot_colors_gradient,
                                         color = myColor_ex2, border_color = "gray60",
                                         main = "Hippocampus Hierarchical Clustering of Top 50 Diet DEGs in 12M APP NLGF Animals",
                                         fontsize = 18,
                                         fontsize_row =  10,
                                         fontsize_col= 10)  ## well it doesnt look great but honestly wasnt expecting it to look fantastic






## ok lets try 3 Month WT
## Trying to make heatmap of top 50 Diet DEGs in 3month WT
gene_expression_padj_ordered_from_dge_dds_WT_hippo_three_gene_se_10filtered
WT3_Diet_DEGs <- gene_expression_padj_ordered_from_dge_dds_WT_hippo_three_gene_se_10filtered
## Use top 50 DEGs
top50 <- row.names(WT3_Diet_DEGs[1:50,])
top50


## Make matrix of 3M WT animal abundances
filtered_WT3_gene_abundances <- as.data.frame(dge_dds_WT_hippo_three_gene_se_10filtered@assays@data@listData[["abundance"]])
filtered_WT3_gene_abundances
filtered_WT3_gene_abundances_matrix <- as.matrix(filtered_WT3_gene_abundances)
filtered_WT3_gene_abundances_matrix

## Switch top 50 list from ENSEMBL names to Gene symbols
top50 <- as.data.frame(top50)
top50
top50_gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, top50$top50, "SYMBOL", keytype = "ENSEMBL")
top50_gene_symbols
## Need to make sure all values in symbol column have  a value and if it is NA then switch to the ENSEMBL name
top50_gene_symbols$SYMBOL[is.na(top50_gene_symbols$SYMBOL)] <- top50_gene_symbols$ENSEMBL[is.na(top50_gene_symbols$SYMBOL)]
top50_gene_symbols
row.names(top50_gene_symbols) <- top50_gene_symbols$SYMBOL
top50_gene_symbols
top50 <- row.names(top50_gene_symbols)
top50

## Subset 3M WT animal abundances so only have top 50 in matrix 
top50_WT3_gene_abundances_matrix <- filtered_WT3_gene_abundances_matrix[top50_gene_symbols$ENSEMBL,]

## Combine top50_gene_symbols and matrix so can easily rename the rownames of the matrix
top50_WT3_gene_abundances_df <- as.data.frame(top50_WT3_gene_abundances_matrix)
top50_WT3_gene_abundances_df$ENSEMBL <- row.names(top50_WT3_gene_abundances_df)
top50_WT3_gene_abundances_df <- merge(top50_WT3_gene_abundances_df, top50_gene_symbols, by='ENSEMBL' )
row.names(top50_WT3_gene_abundances_df) <- top50_WT3_gene_abundances_df$SYMBOL
top50_WT3_gene_abundances_df
top50_WT3_gene_abundances_df <- top50_WT3_gene_abundances_df[,2:13]
top50_WT3_gene_abundances_df
top50_WT3_gene_abundances_matrix <- as.matrix(top50_WT3_gene_abundances_df)
top50_WT3_gene_abundances_matrix

## Make pheatmap
pheatmap(top50_WT3_gene_abundances_matrix)
## this looks awful because things aren't normalized so lets normalize
## make function to z-score
cal_z_score <- function(x) {
  (x- mean(x)) / sd(x)
}

## apply z-score function to top 100 gene abundances
top50_WT3_gene_abundances_matrix_z_score <- t(apply(top50_WT3_gene_abundances_matrix,1, cal_z_score))

## Make heatmap of top 100 DEGs APP vs WT (all ages) in z-score
pheatmap(top50_WT3_gene_abundances_matrix_z_score)


## Make annotation information with metadata
my_sample_col <- data.frame(APP = sample_metadata$APP, Diet= sample_metadata$Diet, Age= sample_metadata$Age)
rownames(my_sample_col) <- sample_metadata$Sample
my_sample_col

## Make top 50 Diet DEGs heatmap
pheatmap(top50_WT3_gene_abundances_matrix_z_score, annotation_col = my_sample_col)
annot_colors_gradient= list(APP=c(APP="red1", WT= "royalblue1"), Diet= c(Control= "gold", Supplemented= "goldenrod3"), Age= c('3'='gray99', '6'='gray75', '9'= 'gray50', '12'= 'gray25'))
palettelength2 <- 50
myColor_ex2 <- colorRampPalette(c("#762A83", "#9970AB", "#C2A5CF", '#E7D4E8', "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B"))(palettelength2)
myColor_ex2
hippo_top50_heatmap_overflow <- pheatmap(top50_WT3_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors = annot_colors_gradient,
                                         color = myColor_ex2, border_color = "gray60",
                                         main = "Hippocampus Hierarchical Clustering of Top 50 Diet DEGs in 3M WT Animals",
                                         fontsize = 18,
                                         fontsize_row =  10,
                                         fontsize_col= 10)
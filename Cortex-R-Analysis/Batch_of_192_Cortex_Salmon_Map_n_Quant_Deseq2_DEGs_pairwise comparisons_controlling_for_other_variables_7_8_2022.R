## Batch of 192- Trimmed Cortex Data
## This is the data that was run on Salmon 6-6-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 7/8/2022
## Using DESeq2 to get DEGs for WT vs APP, trying to control for age, sex, diet
## Will import data using tximeta to create a summarized experiment First
## Using http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions as a guide


library(BiocManager)
library(tximport)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)
library(ggplot2)
library(DESeq2)

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192")
sample_metadata <- read_excel("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Batch_of_192_metadata.xlsx")  #reads in excel file with data about each sample
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
sample_metadata

dir <- setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Trimmed_FASTQ_Aligned_and_Quant_Salmon_6_6_2022")
files <- file.path(dir, sample_metadata$names,  "quant.sf" )
file.exists(files)

cortex_coldata <- data.frame(files, names= sample_metadata$names, APP= sample_metadata$APP, Age= sample_metadata$Age, Diet= sample_metadata$Diet,
                             Sex= sample_metadata$Sex, stringsAsFactors = F)
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
seqinfo(cortex_se)

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
dds_cortex_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ Age + Diet + Sex + APP)
dds_cortex_gene_se$group <- factor(paste0(dds_cortex_gene_se$APP, dds_cortex_gene_se$Age, dds_cortex_gene_se$Diet))
design(dds_cortex_gene_se) <- ~ group
dds_cortex_gene_se_interactions <- DESeq(dds_cortex_gene_se)
resultsNames(dds_cortex_gene_se_interactions)
results(dds_cortex_gene_se_interactions, contrast=c("group", "APP3Control", "APP12Control"))
APP_Control_3_v_12_results <- results(dds_cortex_gene_se_interactions, contrast=c("group", "APP3Control", "APP12Control"))
summary(APP_Control_3_v_12_results)
APP_Control_3_v_12_results_padj_sorted <- APP_Control_3_v_12_results[order(APP_Control_3_v_12_results$padj),]
APP_Control_3_v_12_results_padj_sorted

APP_Supp_v_Control_12mo_results <- results(dds_cortex_gene_se_interactions, contrast=c("group", "APP12Control", "APP12Supplemented"))
APP_Supp_v_Control_12mo_results
summary(APP_Supp_v_Control_12mo_results)
APP_Supp_v_Control_12mo_results_padj_sorted <- APP_Supp_v_Control_12mo_results[order(APP_Supp_v_Control_12mo_results$padj),]
APP_Supp_v_Control_12mo_results_padj_sorted

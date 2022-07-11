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
dds_cortex_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ Age + Diet + Sex + APP) ## this will create a DESEqDataSet with the design saying look for the effects of APP while controlling for age, diet, sex
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_gene_se_10filtered <- dds_cortex_gene_se[rowSums(counts(dds_cortex_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 21234

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_cortex_gene_se_10filtered$APP <- relevel(dds_cortex_gene_se_10filtered$APP, ref = "WT")

dds_cortex_gene_se_10filtered$Age <- relevel(dds_cortex_gene_se_10filtered$Age, ref = "3")

## Run DGE analysis
dge_dds_cortex_gene_se_10filtered <- DESeq(dds_cortex_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_cortex_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
APP_cortex_results <- results(dge_dds_cortex_gene_se_10filtered, contrast = c("APP", "APP", "WT"))
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_cortex_results

## see how many genes reach threshold
summary(results(dge_dds_cortex_gene_se_10filtered, alpha= 0.05))
## 3367 genes are upregulated and 3530 genes are downegulated in APP mice compared to WT with padj < 0.05; this is controlling for age, sex, diet
summary(APP_cortex_results)
## 3972 genes are upregulated and 4120 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_cortex_results_padj_ordered <- APP_cortex_results[order(APP_cortex_results$padj),]
APP_cortex_results_padj_ordered

APP_cortex_results_pval_ordered <- APP_cortex_results[order(APP_cortex_results$pvalue),]
APP_cortex_results_pval_ordered

APP_cortex_results_log2fc_ordered <- APP_cortex_results[order(APP_cortex_results$log2FoldChange),]
APP_cortex_results_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_cortex_results_padj_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_adjust_factors_padj_sorted_7_8_2022.csv")



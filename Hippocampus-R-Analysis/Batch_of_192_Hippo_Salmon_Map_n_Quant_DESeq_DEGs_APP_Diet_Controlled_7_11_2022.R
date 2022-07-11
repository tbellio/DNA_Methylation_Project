## Batch of 192- Trimmed Hippocampus Data
## This is the data that was run on Salmon 6-6-2022 and 6-7-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 7/6/2022
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

dir <- setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Trimmed_FASTQ_Aligned_and_Quant_Salmon_6_6_2022")
files <- file.path(dir, sample_metadata$names,  "quant.sf" )
file.exists(files)

hippo_coldata <- data.frame(files, names= sample_metadata$names, APP= sample_metadata$APP, Age= sample_metadata$Age, Diet= sample_metadata$Diet, Sex= sample_metadata$Sex, stringsAsFactors = F)
hippo_coldata

## Now we have the information for tximeta, lets load it in to create a summarized experiment with all the sample metadata
## This should also find all the matching transcriptome information-- Ensembl Mus musculus release 104
hippo_se <- tximeta(hippo_coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe


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



## Now use DESeq2 to look for DEGs for APP vs WT while controlling for nothing
dds_hippo_gene_se <- DESeqDataSet(hippo_gene_se, design = ~APP) ## this will create a DESEqDataSet with the design saying look for the effects of APP while controlling for age
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_hippo_gene_se_10filtered <- dds_hippo_gene_se[rowSums(counts(dds_hippo_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 21234

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_hippo_gene_se_10filtered$APP <- relevel(dds_hippo_gene_se_10filtered$APP, ref = "WT")

## Run DGE analysis for APP vs WT while controlling for Age
dge_dds_hippo_gene_se_10filtered <- DESeq(dds_hippo_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_hippo_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
APP_hippo_results <- results(dge_dds_hippo_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_hippo_results

## see how many genes reach threshold
summary(APP_hippo_results, alpha= 0.05)
## 2092 genes are upregulated and 2243 genes are downegulated in APP mice compared to WT with padj < 0.05
summary(APP_hippo_results)
## 2670 genes are upregulated and 2738 genes are downregulated in APP mice compared to WT with padj < 0.1 

## Sort the gene expression data by FDR(adjusted p-value)
APP_hippo_results_age_controlled_padj_ordered <- APP_hippo_age_controlled_results[order(APP_hippo_age_controlled_results$padj),]
APP_hippo_results_age_controlled_padj_ordered

APP_hippo_results_age_controlled_pval_ordered <- APP_hippo_age_controlled_results[order(APP_hippo_age_controlled_results$pvalue),]
APP_hippo_results_age_controlled_pval_ordered

APP_hippo_results_age_controlled_log2fc_ordered <- APP_hippo_age_controlled_results[order(APP_hippo_age_controlled_results$log2FoldChange),]
APP_hippo_results_age_controlled_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_hippo_results_age_controlled_padj_ordered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_APP_v_WT_DGE_Analysis_age_controlled_padj_sorted_7_11_2022.csv")





## Now use DESeq2 to look for DEGs for APP vs WT while controlling for Age
library('DESeq2')
dds_hippo_gene_se <- DESeqDataSet(hippo_gene_se, design = ~ Age + APP) ## this will create a DESEqDataSet with the design saying look for the effects of APP while controlling for age
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_hippo_gene_se_10filtered <- dds_hippo_gene_se[rowSums(counts(dds_hippo_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 21234

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_hippo_gene_se_10filtered$APP <- relevel(dds_hippo_gene_se_10filtered$APP, ref = "WT")

## Run DGE analysis for APP vs WT while controlling for Age
dge_dds_hippo_gene_se_10filtered <- DESeq(dds_hippo_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_hippo_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
APP_hippo_age_controlled_results <- results(dge_dds_hippo_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_hippo_age_controlled_results

## see how many genes reach threshold
summary(APP_hippo_age_controlled_results, alpha= 0.05)
## 2092 genes are upregulated and 2243 genes are downegulated in APP mice compared to WT with padj < 0.05
summary(APP_hippo_age_controlled_results)
## 2742 genes are upregulated and 2919 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_hippo_results_age_controlled_padj_ordered <- APP_hippo_age_controlled_results[order(APP_hippo_age_controlled_results$padj),]
APP_hippo_results_age_controlled_padj_ordered

APP_hippo_results_age_controlled_pval_ordered <- APP_hippo_age_controlled_results[order(APP_hippo_age_controlled_results$pvalue),]
APP_hippo_results_age_controlled_pval_ordered

APP_hippo_results_age_controlled_log2fc_ordered <- APP_hippo_age_controlled_results[order(APP_hippo_age_controlled_results$log2FoldChange),]
APP_hippo_results_age_controlled_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_hippo_results_age_controlled_padj_ordered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_APP_v_WT_DGE_Analysis_age_controlled_padj_sorted_7_11_2022.csv")






## Now Do APP vs WT while controlling for diet
dds_hippo_app_diet_control_gene_se <- DESeqDataSet(hippo_gene_se, design = ~ Diet + APP)

## Filter so only genes that have more than 10 reads total are included
dds_hippo_app_diet_control_gene_se_10filtered <- dds_hippo_app_diet_control_gene_se[rowSums(counts(dds_hippo_gene_se)) >= 10,]
## Filtered out ~15000 genes (35682 to 20480)

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_hippo_app_diet_control_gene_se_10filtered$APP <- relevel(dds_hippo_app_diet_control_gene_se_10filtered$APP, ref = "WT")

## Run DGE analysis for APP vs WT while controlling for Diet
dge_dds_hippo_app_diet_control_gene_se_10filtered <- DESeq(dds_hippo_app_diet_control_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_hippo_app_diet_control_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
APP_hippo_diet_controlled_results <- results(dge_dds_hippo_app_diet_control_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_hippo_diet_controlled_results

## see how many genes reach threshold
summary(APP_hippo_diet_controlled_results, alpha= 0.05)
## 1989 genes are upregulated and 2043 genes are downegulated in APP mice compared to WT with padj < 0.05
summary(APP_hippo_diet_controlled_results)
## 2614 genes are upregulated and 2714 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_hippo_results_diet_controlled_padj_ordered <- APP_hippo_diet_controlled_results[order(APP_hippo_diet_controlled_results$padj),]
APP_hippo_results_diet_controlled_padj_ordered

APP_hippo_results_diet_controlled_pval_ordered <- APP_hippo_diet_controlled_results[order(APP_hippo_diet_controlled_results$pvalue),]
APP_hippo_results_diet_controlled_pval_ordered

APP_hippo_results_diet_controlled_log2fc_ordered <- APP_hippo_diet_controlled_results[order(APP_hippo_diet_controlled_results$log2FoldChange),]
APP_hippo_results_diet_controlled_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_hippo_results_diet_controlled_padj_ordered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_APP_v_WT_DGE_Analysis_diet_controlled_padj_sorted_7_11_2022.csv")









































## Now Do APP vs WT while controlling for sex
dds_hippo_app_sex_control_gene_se <- DESeqDataSet(hippo_gene_se, design = ~ Sex + APP)

## Filter so only genes that have more than 10 reads total are included
dds_hippo_app_sex_control_gene_se_10filtered <- dds_hippo_app_sex_control_gene_se[rowSums(counts(dds_hippo_gene_se)) >= 10,]
## Filtered out ~15000 genes (35682 to 20480)

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_hippo_app_sex_control_gene_se_10filtered$APP <- relevel(dds_hippo_app_sex_control_gene_se_10filtered$APP, ref = "WT")

## Run DGE analysis for APP vs WT while controlling for Diet
dge_dds_hippo_app_sex_control_gene_se_10filtered <- DESeq(dds_hippo_app_sex_control_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_hippo_app_sex_control_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
APP_hippo_sex_controlled_results <- results(dge_dds_hippo_app_sex_control_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_hippo_sex_controlled_results

## see how many genes reach threshold
summary(APP_hippo_sex_controlled_results, alpha= 0.05)
## 2029 genes are upregulated and 2062 genes are downegulated in APP mice compared to WT with padj < 0.05
summary(APP_hippo_sex_controlled_results)
## 2661 genes are upregulated and 2739 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_hippo_results_sex_controlled_padj_ordered <- APP_hippo_sex_controlled_results[order(APP_hippo_sex_controlled_results$padj),]
APP_hippo_results_sex_controlled_padj_ordered

APP_hippo_results_sex_controlled_pval_ordered <- APP_hippo_sex_controlled_results[order(APP_hippo_sex_controlled_results$pvalue),]
APP_hippo_results_sex_controlled_pval_ordered

APP_hippo_results_sex_controlled_log2fc_ordered <- APP_hippo_sex_controlled_results[order(APP_hippo_sex_controlled_results$log2FoldChange),]
APP_hippo_results_sex_controlled_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_hippo_results_sex_controlled_padj_ordered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_APP_v_WT_DGE_Analysis_sex_controlled_padj_sorted_7_11_2022.csv")

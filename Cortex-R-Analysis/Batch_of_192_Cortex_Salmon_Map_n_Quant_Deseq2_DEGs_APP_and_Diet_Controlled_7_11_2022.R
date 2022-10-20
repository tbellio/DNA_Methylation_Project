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
library(clusterProfiler)

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




## Now use DESeq2 to look for DEGs in APP vs WT while controlling for nothing
dds_cortex_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ APP) ## this will create a DESEqDataSet with the design saying look for the effects of APP while controlling for age
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

## Look at results and create gene expression table as DESeqResults
APP_cortex_results <- results(dge_dds_cortex_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_cortex_results

## see how many genes reach threshold
summary(APP_cortex_results, alpha= 0.05)
## 3320 genes are upregulated and 3457 genes are downegulated in APP mice compared to WT with padj < 0.05; this is controlling for age, sex, diet
summary(APP_cortex_results)
## 3931 genes are upregulated and 4077 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_cortex_age_controlled_results_padj_ordered <- APP_cortex_age_controlled_results[order(APP_cortex_age_controlled_results$padj),]
APP_cortex_age_controlled_results_padj_ordered

APP_cortex_age_controlled_results_pval_ordered <- APP_cortex_age_controlled_results[order(APP_cortex_age_controlled_results$pvalue),]
APP_cortex_age_controlled_results_pval_ordered

APP_cortex_age_controlled_results_log2fc_ordered <- APP_cortex_age_controlled_results[order(APP_cortex_age_controlled_results$log2FoldChange),]
APP_cortex_age_controlled_results_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_cortex_age_controlled_results_padj_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_age_controlled_padj_sorted_7_11_2022.csv")





## Now use DESeq2 to look for DEGs in APP vs WT while controlling for Age
library('DESeq2')
dds_cortex_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ Age + APP) ## this will create a DESEqDataSet with the design saying look for the effects of APP while controlling for age
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

## Look at results and create gene expression table as DESeqResults
APP_cortex_age_controlled_results <- results(dge_dds_cortex_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_cortex_age_controlled_results

## see how many genes reach threshold
summary(APP_cortex_age_controlled_results, alpha= 0.05)
## 3342 genes are upregulated and 3516 genes are downegulated in APP mice compared to WT with padj < 0.05; this is controlling for age, sex, diet
summary(APP_cortex_age_controlled_results)
## 3959 genes are upregulated and 4120 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_cortex_age_controlled_results_padj_ordered <- APP_cortex_age_controlled_results[order(APP_cortex_age_controlled_results$padj),]
APP_cortex_age_controlled_results_padj_ordered

APP_cortex_age_controlled_results_pval_ordered <- APP_cortex_age_controlled_results[order(APP_cortex_age_controlled_results$pvalue),]
APP_cortex_age_controlled_results_pval_ordered

APP_cortex_age_controlled_results_log2fc_ordered <- APP_cortex_age_controlled_results[order(APP_cortex_age_controlled_results$log2FoldChange),]
APP_cortex_age_controlled_results_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_cortex_age_controlled_results_padj_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_age_controlled_padj_sorted_7_11_2022.csv")





## Now use DESeq2 to look for DEGs in APP vs WT while controlling for Diet
dds_cortex_app_diet_control_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ Diet + APP) ## this will create a DESEqDataSet with the design saying look for the effects of APP while controlling for age
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_app_diet_control_gene_se_10filtered <- dds_cortex_app_diet_control_gene_se[rowSums(counts(dds_cortex_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 21234

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_cortex_app_diet_control_gene_se_10filtered$APP <- relevel(dds_cortex_app_diet_control_gene_se_10filtered$APP, ref = "WT")

## Run DGE analysis
dge_dds_cortex_app_diet_control_gene_se_10filtered <- DESeq(dds_cortex_app_diet_control_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_cortex_app_diet_control_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
APP_cortex_diet_controlled_results <- results(dge_dds_cortex_app_diet_control_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_cortex_diet_controlled_results

## see how many genes reach threshold
summary(APP_cortex_diet_controlled_results, alpha= 0.05)
## 3321 genes are upregulated and 3453 genes are downegulated in APP mice compared to WT with padj < 0.05; this is controlling for age, sex, diet
summary(APP_cortex_diet_controlled_results)
## 3928 genes are upregulated and 4079 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_cortex_diet_controlled_results_padj_ordered <- APP_cortex_diet_controlled_results[order(APP_cortex_diet_controlled_results$padj),]
APP_cortex_diet_controlled_results_padj_ordered

APP_cortex_diet_controlled_results_pval_ordered <- APP_cortex_diet_controlled_results[order(APP_cortex_diet_controlled_results$pvalue),]
APP_cortex_diet_controlled_results_pval_ordered

APP_cortex_diet_controlled_results_log2fc_ordered <- APP_cortex_diet_controlled_results[order(APP_cortex_diet_controlled_results$log2FoldChange),]
APP_cortex_diet_controlled_results_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_cortex_diet_controlled_results_padj_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_diet_controlled_padj_sorted_7_11_2022.csv")





## Now use DESeq2 to look for DEGs in APP vs WT while controlling for Sex
dds_cortex_app_sex_control_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ Sex + APP) ## this will create a DESEqDataSet with the design saying look for the effects of APP while controlling for age
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_app_sex_control_gene_se_10filtered <- dds_cortex_app_sex_control_gene_se[rowSums(counts(dds_cortex_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 21234

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_cortex_app_sex_control_gene_se_10filtered$APP <- relevel(dds_cortex_app_sex_control_gene_se_10filtered$APP, ref = "WT")

## Run DGE analysis
dge_dds_cortex_app_sex_control_gene_se_10filtered <- DESeq(dds_cortex_app_sex_control_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_cortex_app_sex_control_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
APP_cortex_sex_controlled_results <- results(dge_dds_cortex_app_sex_control_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_cortex_sex_controlled_results

## see how many genes reach threshold
summary(APP_cortex_sex_controlled_results, alpha= 0.05)
## 3349 genes are upregulated and 3485 genes are downegulated in APP mice compared to WT with padj < 0.05; this is controlling for age, sex, diet
summary(APP_cortex_sex_controlled_results)
## 3965 genes are upregulated and 4104 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_cortex_diet_controlled_results_padj_ordered <- APP_cortex_diet_controlled_results[order(APP_cortex_diet_controlled_results$padj),]
APP_cortex_diet_controlled_results_padj_ordered

APP_cortex_diet_controlled_results_pval_ordered <- APP_cortex_diet_controlled_results[order(APP_cortex_diet_controlled_results$pvalue),]
APP_cortex_diet_controlled_results_pval_ordered

APP_cortex_diet_controlled_results_log2fc_ordered <- APP_cortex_diet_controlled_results[order(APP_cortex_diet_controlled_results$log2FoldChange),]
APP_cortex_diet_controlled_results_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_cortex_diet_controlled_results_padj_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_diet_controlled_padj_sorted_7_11_2022.csv")































## Now use DESeq2 to look for DEGs in APP vs WT while controlling for Sex, age, and diet
dds_cortex_app_all_control_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ Sex + Age + Diet+ APP) ## this will create a DESEqDataSet with the design saying look for the effects of APP while controlling for age
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_app_all_control_gene_se_10filtered <- dds_cortex_app_all_control_gene_se[rowSums(counts(dds_cortex_app_all_control_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 21234

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_cortex_app_all_control_gene_se_10filtered$APP <- relevel(dds_cortex_app_all_control_gene_se_10filtered$APP, ref = "WT")

## Run DGE analysis
dge_dds_cortex_app_all_control_gene_se_10filtered <- DESeq(dds_cortex_app_all_control_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just wt vs app
resultsNames(dge_dds_cortex_app_all_control_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
APP_cortex_all_controlled_results <- results(dge_dds_cortex_app_all_control_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so App vs wt for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
APP_cortex_all_controlled_results

## see how many genes reach threshold
summary(APP_cortex_all_controlled_results, alpha= 0.05)
## 3367 genes are upregulated and 3530 genes are downegulated in APP mice compared to WT with padj < 0.05; this is controlling for age, sex, diet
summary(APP_cortex_all_controlled_results)
## 3972 genes are upregulated and 4120 genes are downregulated in APP mice compared to WT with padj < 0.1; this is controlling for age, sex, diet

## Sort the gene expression data by FDR(adjusted p-value)
APP_cortex_all_controlled_results_padj_ordered <- APP_cortex_all_controlled_results[order(APP_cortex_all_controlled_results$padj),]
APP_cortex_all_controlled_results_padj_ordered

APP_cortex_all_controlled_results_pval_ordered <- APP_cortex_all_controlled_results[order(APP_cortex_all_controlled_results$pvalue),]
APP_cortex_all_controlled_results_pval_ordered

APP_cortex_all_controlled_results_log2fc_ordered <- APP_cortex_all_controlled_results[order(APP_cortex_all_controlled_results$log2FoldChange),]
APP_cortex_all_controlled_results_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(APP_cortex_all_controlled_results_padj_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_7_11_2022.csv")


## Changing Ensembl Ids to Gene Symbols so easier to Read
symbols <- bitr(rownames(APP_cortex_all_controlled_results_padj_ordered), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
APP_cortex_all_controlled_results_padj_ordered_df <- as.data.frame(APP_cortex_all_controlled_results_padj_ordered)
APP_cortex_all_controlled_results_padj_ordered_df$ENSEMBL <- rownames(APP_cortex_all_controlled_results_padj_ordered_df)
APP_cortex_all_controlled_results_padj_ordered_df
merged_symbols_APP_cortex_all_controlled_results_padj_ordered <- merge(APP_cortex_all_controlled_results_padj_ordered_df, symbols, by.x="ENSEMBL", by.y="ENSEMBL")
merged_symbols_APP_cortex_all_controlled_results_padj_ordered <- merged_symbols_APP_cortex_all_controlled_results_padj_ordered[order(merged_symbols_APP_cortex_all_controlled_results_padj_ordered$padj),]
merged_symbols_APP_cortex_all_controlled_results_padj_ordered
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(merged_symbols_APP_cortex_all_controlled_results_padj_ordered, file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_with_symbols_8_10_2022.csv")




























## Effect of DIET

## Now use DESeq2 to look for DEGs for Control vs Supp while controlling for nothing
dds_cortex_diet_gene_se <- DESeqDataSet(cortex_gene_se, design = ~Diet) ## this will create a DESEqDataSet with the design saying look for the effects of diet while controlling for nothing


## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_diet_gene_se_10filtered <- dds_cortex_diet_gene_se[rowSums(counts(dds_cortex_diet_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 21234

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_cortex_diet_gene_se_10filtered$Diet <- relevel(dds_cortex_diet_gene_se_10filtered$Diet, ref = "Control")

## Run DGE analysis for Control vs Supp; no controlling
dge_dds_cortex_diet_gene_se_10filtered <- DESeq(dds_cortex_diet_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just control vs supp
resultsNames(dge_dds_cortex_diet_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
Diet_cortex_results <- results(dge_dds_cortex_diet_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so Control vs Supp for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
Diet_cortex_results

## see how many genes reach threshold
summary(Diet_cortex_results, alpha= 0.05)
## 0 genes are upregulated and 0 genes are downegulated in supplemented mice compared to control with padj < 0.05
summary(Diet_cortex_results)
## 0 genes are upregulated and 0 genes are downregulated in supplemented mice compared to control with padj < 0.1 

## Sort the gene expression data by FDR(adjusted p-value)
Diet_cortex_results_padj_ordered <- Diet_cortex_results[order(Diet_cortex_results$padj),]
Diet_cortex_results_padj_ordered

Diet_cortex_results_pval_ordered <- Diet_cortex_results[order(Diet_cortex_results$pvalue),]
Diet_cortex_results_pval_ordered

Diet_cortex_results_log2fc_ordered <- Diet_cortex_results[order(Diet_cortex_results$log2FoldChange),]
Diet_cortex_results_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(Diet_cortex_results_padj_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Supp_v_Control_DGE_Analysis_no_control_padj_sorted_7_12_2022.csv")




















## Now use DESeq2 to look for DEGs for Control vs Supp while controlling for all other variables
dds_cortex_diet_all_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ Age + Sex + APP + Diet) ## this will create a DESEqDataSet with the design saying look for the effects of diet while controlling for all other variables


## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_diet_all_gene_se_10filtered <- dds_cortex_diet_all_gene_se[rowSums(counts(dds_cortex_diet_all_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 20480

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_cortex_diet_all_gene_se_10filtered$Diet <- relevel(dds_cortex_diet_all_gene_se_10filtered$Diet, ref = "Control")

## Run DGE analysis for Control vs Supp; no controlling
dge_dds_cortex_diet_all_gene_se_10filtered <- DESeq(dds_cortex_diet_all_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just control vs supp
resultsNames(dge_dds_cortex_diet_all_gene_se_10filtered)

## Look at results and create gene expression table as DESeqResults
Diet_all_controlled_cortex_results <- results(dge_dds_cortex_diet_all_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so Control vs Supp for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
Diet_all_controlled_cortex_results

## see how many genes reach threshold
summary(Diet_all_controlled_cortex_results, alpha= 0.05)
## 2 genes are upregulated and 45 genes are downegulated in supplemented mice compared to control with padj < 0.05
summary(Diet_all_controlled_cortex_results)
## 5 genes are upregulated and 48 genes are downregulated in supplemented mice compared to control with padj < 0.1 

## Sort the gene expression data by FDR(adjusted p-value)
Diet_all_controlled_cortex_results_padj_ordered <- Diet_all_controlled_cortex_results[order(Diet_all_controlled_cortex_results$padj),]
Diet_all_controlled_cortex_results_padj_ordered

Diet_all_controlled_cortex_results_pval_ordered <- Diet_all_controlled_cortex_results[order(Diet_all_controlled_cortex_results$pvalue),]
Diet_all_controlled_cortex_results_pval_ordered

Diet_all_controlled_cortex_results_log2fc_ordered <- Diet_all_controlled_cortex_results[order(Diet_all_controlled_cortex_results$log2FoldChange),]
Diet_all_controlled_cortex_results_log2fc_ordered

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(Diet_all_controlled_cortex_results_padj_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Supp_v_Control_DGE_Analysis_all_controlled_padj_sorted_7_12_2022.csv")

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(Diet_all_controlled_cortex_results_pval_ordered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Supp_v_Control_DGE_Analysis_all_controlled_pval_sorted_7_12_2022.csv")




















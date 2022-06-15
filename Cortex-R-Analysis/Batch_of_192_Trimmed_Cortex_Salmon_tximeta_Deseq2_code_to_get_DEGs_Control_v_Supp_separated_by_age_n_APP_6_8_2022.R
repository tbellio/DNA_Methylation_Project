## Batch of 192- Trimmed Cortex Data 
## This is the data that was run on Salmon 6-6-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 6/8/2022
## Will import data using tximeta to create a summarized experiment First
## Using https://www.reneshbedre.com/blog/deseq2.html and https://lashlock.github.io/compbio/R_presentation.html as a guide
## Trying to Use DESeq2 to get DEGs for Control vs Supp Diet at all time points in separated APP and WT (so a comparison of control vs supp at 3, 6, 9, and 12 months old in separately APP and WT mice)
## Adding export to excel files 6-10-22




library(BiocManager)
library(tximport)
library(tximportData)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)
library(ggplot2)

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

## Create dataframe that has the location of the quant files and all sample metadata
cortex_coldata <- data.frame(files, names= sample_metadata$names, APP= sample_metadata$APP, Age= sample_metadata$Age, Diet= sample_metadata$Diet, Sex= sample_metadata$Sex, stringsAsFactors = F)
cortex_coldata

## Subset data to only include 3, 6, 9, and 12 month animals
## Make sure only 24 rows for each
cortex_three_only_coldata <- cortex_coldata[cortex_coldata$Age== c("3"), ]
cortex_six_only_coldata <- cortex_coldata[cortex_coldata$Age== c("6"),]
cortex_nine_only_coldata <- cortex_coldata[cortex_coldata$Age== c("9"),]
cortex_twelve_only_coldata <- cortex_coldata[cortex_coldata$Age== c("12"),]

##Subset for APP Only and WT Only
APP_cortex_coldata <- cortex_coldata[cortex_coldata$APP==c("APP"),]
APP_cortex_three_only_coldata <- cortex_three_only_coldata[cortex_three_only_coldata$APP== c("APP"),]
APP_cortex_six_only_coldata <- cortex_six_only_coldata[cortex_six_only_coldata$APP== c("APP"),]
APP_cortex_nine_only_coldata <- cortex_nine_only_coldata[cortex_nine_only_coldata$APP== c("APP"),]
APP_cortex_twelve_only_coldata <- cortex_twelve_only_coldata[cortex_twelve_only_coldata$APP== c("APP"),]

WT_cortex_coldata <- cortex_coldata[cortex_coldata$APP==c("WT"),]
WT_cortex_three_only_coldata <- cortex_three_only_coldata[cortex_three_only_coldata$APP== c("WT"),]
WT_cortex_six_only_coldata <- cortex_six_only_coldata[cortex_six_only_coldata$APP== c("WT"),]
WT_cortex_nine_only_coldata <- cortex_nine_only_coldata[cortex_nine_only_coldata$APP== c("WT"),]
WT_cortex_twelve_only_coldata <- cortex_twelve_only_coldata[cortex_twelve_only_coldata$APP== c("WT"),]


## Now we have the information for tximeta, lets load it in to create a summarized experiment with all the sample metadata
## This should also find all the matching transcriptome information-- Ensembl Mus musculus release 104
APP_cortex_se <- tximeta(APP_cortex_coldata)
APP_cortex__three_se <- tximeta(APP_cortex_three_only_coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe
APP_cortex_six_se <- tximeta(APP_cortex_six_only_coldata)
APP_cortex_nine_se <- tximeta(APP_cortex_nine_only_coldata)
APP_cortex_twelve_se <- tximeta(APP_cortex_twelve_only_coldata)

WT_cortex_se <- tximeta(WT_cortex_coldata)
WT_cortex__three_se <- tximeta(WT_cortex_three_only_coldata)
WT_cortex_six_se <- tximeta(WT_cortex_six_only_coldata)
WT_cortex_nine_se <- tximeta(WT_cortex_nine_only_coldata)
WT_cortex_twelve_se <- tximeta(WT_cortex_twelve_only_coldata)



suppressPackageStartupMessages(library(SummarizedExperiment))
##adding exons to the summarized experiment (hippo__three_se, etc.)
APP_cortex_se.exons <- addExons(APP_cortex_se)
APP_cortex_three_se.exons <- addExons(APP_cortex__three_se)
APP_cortex_six_se.exons <- addExons(APP_cortex_six_se)
APP_cortex_nine_se.exons <- addExons(APP_cortex_nine_se)
APP_cortex_twelve_se.exons <- addExons(APP_cortex_twelve_se)

WT_cortex_se.exons <- addExons(WT_cortex_se)
WT_cortex_three_se.exons <- addExons(WT_cortex__three_se)
WT_cortex_six_se.exons <- addExons(WT_cortex_six_se)
WT_cortex_nine_se.exons <- addExons(WT_cortex_nine_se)
WT_cortex_twelve_se.exons <- addExons(WT_cortex_twelve_se)

##Check rowranges--shouldnt have to make sure it works for all four as same code
rowRanges(APP_cortex_three_se.exons)

## Summarize read counts to the gene level
APP_cortex_gene_se <- summarizeToGene(APP_cortex_se.exons)
APP_cortex_three_gene_se <- summarizeToGene(APP_cortex_three_se.exons)
APP_cortex_six_gene_se <- summarizeToGene(APP_cortex_six_se.exons)
APP_cortex_nine_gene_se <- summarizeToGene(APP_cortex_nine_se.exons)
APP_cortex_twelve_gene_se <- summarizeToGene(APP_cortex_twelve_se.exons)

WT_cortex_gene_se <- summarizeToGene(WT_cortex_se.exons)
WT_cortex_three_gene_se <- summarizeToGene(WT_cortex_three_se.exons)
WT_cortex_six_gene_se <- summarizeToGene(WT_cortex_six_se.exons)
WT_cortex_nine_gene_se <- summarizeToGene(WT_cortex_nine_se.exons)
WT_cortex_twelve_gene_se <- summarizeToGene(WT_cortex_twelve_se.exons)

rowRanges(APP_cortex_three_gene_se)


## Now use DESeq2 to look for DEGs at each time point
library('DESeq2')
dds_APP_cortex_gene_se <- DESeqDataSet(APP_cortex_gene_se, design = ~ Diet)
dds_APP_cortex_three_gene_se <- DESeqDataSet(APP_cortex_three_gene_se, design = ~ Diet)
dds_APP_cortex_six_gene_se <- DESeqDataSet(APP_cortex_six_gene_se, design = ~ Diet)
dds_APP_cortex_nine_gene_se <- DESeqDataSet(APP_cortex_nine_gene_se, design = ~ Diet)
dds_APP_cortex_twelve_gene_se <- DESeqDataSet(APP_cortex_twelve_gene_se, design = ~ Diet)

dds_WT_cortex_gene_se <- DESeqDataSet(WT_cortex_gene_se, design = ~ Diet)
dds_WT_cortex_three_gene_se <- DESeqDataSet(WT_cortex_three_gene_se, design = ~ Diet)
dds_WT_cortex_six_gene_se <- DESeqDataSet(WT_cortex_six_gene_se, design = ~ Diet)
dds_WT_cortex_nine_gene_se <- DESeqDataSet(WT_cortex_nine_gene_se, design = ~ Diet)
dds_WT_cortex_twelve_gene_se <- DESeqDataSet(WT_cortex_twelve_gene_se, design = ~ Diet)
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; if not it will change it to a factor (this code has it as a character so it knows to change it)

## Filter so that only genes that have more than 10 reads are inculuded
dds_APP_cortex_gene_se_10filtered <- dds_APP_cortex_gene_se[rowSums(counts(dds_APP_cortex_gene_se)) >= 10,]
dds_APP_cortex_three_gene_se_10filtered <- dds_APP_cortex_three_gene_se[rowSums(counts(dds_APP_cortex_three_gene_se)) >= 10,]
dds_APP_cortex_six_gene_se_10filtered <- dds_APP_cortex_six_gene_se[rowSums(counts(dds_APP_cortex_six_gene_se)) >= 10,]
dds_APP_cortex_nine_gene_se_10filtered <- dds_APP_cortex_nine_gene_se[rowSums(counts(dds_APP_cortex_nine_gene_se)) >= 10,]
dds_APP_cortex_twelve_gene_se_10filtered <- dds_APP_cortex_twelve_gene_se[rowSums(counts(dds_APP_cortex_twelve_gene_se)) >= 10,]

dds_WT_cortex_gene_se_10filtered <- dds_WT_cortex_gene_se[rowSums(counts(dds_WT_cortex_gene_se)) >= 10,]
dds_WT_cortex_three_gene_se_10filtered <- dds_WT_cortex_three_gene_se[rowSums(counts(dds_WT_cortex_three_gene_se)) >= 10,]
dds_WT_cortex_six_gene_se_10filtered <- dds_WT_cortex_six_gene_se[rowSums(counts(dds_WT_cortex_six_gene_se)) >= 10,]
dds_WT_cortex_nine_gene_se_10filtered <- dds_WT_cortex_nine_gene_se[rowSums(counts(dds_WT_cortex_nine_gene_se)) >= 10,]
dds_WT_cortex_twelve_gene_se_10filtered <- dds_WT_cortex_twelve_gene_se[rowSums(counts(dds_WT_cortex_twelve_gene_se)) >= 10,]
## this filtered out ~18000 genes; went from 35682 to 18000-19000

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the contorl diet mice
dds_APP_cortex_gene_se_10filtered$Diet <- relevel(dds_APP_cortex_gene_se_10filtered$Diet, ref = "Control")
dds_APP_cortex_three_gene_se_10filtered$Diet <- relevel(dds_APP_cortex_three_gene_se_10filtered$Diet, ref = "Control")
dds_APP_cortex_six_gene_se_10filtered$Diet <- relevel(dds_APP_cortex_six_gene_se_10filtered$Diet, ref = "Control")
dds_APP_cortex_nine_gene_se_10filtered$Diet <- relevel(dds_APP_cortex_nine_gene_se_10filtered$Diet, ref = "Control")
dds_APP_cortex_twelve_gene_se_10filtered$Diet <- relevel(dds_APP_cortex_twelve_gene_se_10filtered$Diet, ref = "Control")

dds_WT_cortex_gene_se_10filtered$Diet <- relevel(dds_WT_cortex_gene_se_10filtered$Diet, ref = "Control" )
dds_WT_cortex_three_gene_se_10filtered$Diet <- relevel(dds_WT_cortex_three_gene_se_10filtered$Diet, ref = "Control")
dds_WT_cortex_six_gene_se_10filtered$Diet <- relevel(dds_WT_cortex_six_gene_se_10filtered$Diet, ref = "Control")
dds_WT_cortex_nine_gene_se_10filtered$Diet <- relevel(dds_WT_cortex_nine_gene_se_10filtered$Diet, ref = "Control")
dds_WT_cortex_twelve_gene_se_10filtered$Diet <- relevel(dds_WT_cortex_twelve_gene_se_10filtered$Diet, ref = "Control")

## Run DGE analysis
dge_dds_APP_cortex_gene_se_10filtered <- DESeq(dds_APP_cortex_gene_se_10filtered)
dge_dds_APP_cortex_three_gene_se_10filtered <- DESeq(dds_APP_cortex_three_gene_se_10filtered)
dge_dds_APP_cortex_six_gene_se_10filtered <- DESeq(dds_APP_cortex_six_gene_se_10filtered)
dge_dds_APP_cortex_nine_gene_se_10filtered <- DESeq(dds_APP_cortex_nine_gene_se_10filtered)
dge_dds_APP_cortex_twelve_gene_se_10filtered <- DESeq(dds_APP_cortex_twelve_gene_se_10filtered)

dge_dds_WT_cortex_gene_se_10filtered <- DESeq(dds_WT_cortex_gene_se_10filtered)
dge_dds_WT_cortex_three_gene_se_10filtered <- DESeq(dds_WT_cortex_three_gene_se_10filtered)
dge_dds_WT_cortex_six_gene_se_10filtered <- DESeq(dds_WT_cortex_six_gene_se_10filtered)
dge_dds_WT_cortex_nine_gene_se_10filtered <- DESeq(dds_WT_cortex_nine_gene_se_10filtered)
dge_dds_WT_cortex_twelve_gene_se_10filtered <- DESeq(dds_WT_cortex_twelve_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just Control vs Supp
resultsNames(dge_dds_APP_cortex_gene_se_10filtered)
resultsNames(dge_dds_APP_cortex_three_gene_se_10filtered)
resultsNames(dge_dds_APP_cortex_six_gene_se_10filtered)
resultsNames(dge_dds_APP_cortex_nine_gene_se_10filtered)
resultsNames(dge_dds_APP_cortex_twelve_gene_se_10filtered)

resultsNames(dge_dds_WT_cortex_gene_se_10filtered)
resultsNames(dge_dds_WT_cortex_three_gene_se_10filtered)
resultsNames(dge_dds_WT_cortex_six_gene_se_10filtered)
resultsNames(dge_dds_WT_cortex_nine_gene_se_10filtered)
resultsNames(dge_dds_WT_cortex_twelve_gene_se_10filtered)

## Get gene expression table as DESeqResults
gene_expression_from_dge_dds_APP_cortex_gene_se_10filtered <- results(dge_dds_APP_cortex_gene_se_10filtered) 
gene_expression_from_dge_dds_APP_cortex_three_gene_se_10filtered <- results(dge_dds_APP_cortex_three_gene_se_10filtered)
gene_expression_from_dge_dds_APP_cortex_six_gene_se_10filtered <- results(dge_dds_APP_cortex_six_gene_se_10filtered)
gene_expression_from_dge_dds_APP_cortex_nine_gene_se_10filtered <- results(dge_dds_APP_cortex_nine_gene_se_10filtered)
gene_expression_from_dge_dds_APP_cortex_twelve_gene_se_10filtered <- results(dge_dds_APP_cortex_twelve_gene_se_10filtered)

gene_expression_from_dge_dds_WT_cortex_gene_se_10filtered <- results(dge_dds_WT_cortex_gene_se_10filtered)
gene_expression_from_dge_dds_WT_cortex_three_gene_se_10filtered <- results(dge_dds_WT_cortex_three_gene_se_10filtered)
gene_expression_from_dge_dds_WT_cortex_six_gene_se_10filtered <- results(dge_dds_WT_cortex_six_gene_se_10filtered)
gene_expression_from_dge_dds_WT_cortex_nine_gene_se_10filtered <- results(dge_dds_WT_cortex_nine_gene_se_10filtered)
gene_expression_from_dge_dds_WT_cortex_twelve_gene_se_10filtered <- results(dge_dds_WT_cortex_twelve_gene_se_10filtered)

## See what the gene expression table says
## Will tell you what it was comparing (so APP vs WT for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
gene_expression_from_dge_dds_APP_cortex_gene_se_10filtered
gene_expression_from_dge_dds_APP_cortex_three_gene_se_10filtered
gene_expression_from_dge_dds_APP_cortex_six_gene_se_10filtered
gene_expression_from_dge_dds_APP_cortex_nine_gene_se_10filtered
gene_expression_from_dge_dds_APP_cortex_twelve_gene_se_10filtered

## Sort the gene expression data by FDR(adjusted p-value)
gene_expression_padj_ordered_from_dge_dds_APP_cortex_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_APP_cortex_three_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_three_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_three_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_APP_cortex_six_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_six_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_six_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_APP_cortex_nine_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_nine_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_APP_cortex_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_twelve_gene_se_10filtered$padj),]

gene_expression_padj_ordered_from_dge_dds_WT_cortex_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_WT_cortex_three_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_three_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_three_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_WT_cortex_six_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_six_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_six_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_WT_cortex_nine_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_nine_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_WT_cortex_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_twelve_gene_se_10filtered$padj),]

gene_expression_padj_ordered_from_dge_dds_APP_cortex_gene_se_10filtered 
gene_expression_padj_ordered_from_dge_dds_WT_cortex_gene_se_10filtered

## Sort by p-value
gene_expression_pval_ordered_from_dge_dds_APP_cortex_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_APP_cortex_three_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_three_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_three_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_APP_cortex_six_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_six_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_six_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_APP_cortex_nine_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_nine_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_APP_cortex_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_APP_cortex_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_APP_cortex_twelve_gene_se_10filtered$pvalue),]

gene_expression_pval_ordered_from_dge_dds_WT_cortex_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_WT_cortex_three_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_three_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_three_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_WT_cortex_six_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_six_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_six_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_WT_cortex_nine_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_nine_gene_se_10filtered$pvalue),]
gene_expression_pval_ordered_from_dge_dds_WT_cortex_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_WT_cortex_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_WT_cortex_twelve_gene_se_10filtered$pvalue),]


## Get a summary of DEGs with FDR (padj) of <0.05 for all ages combined (separate WT and APP)
summary(results(dge_dds_WT_cortex_gene_se_10filtered, alpha = 0.05))
## no genes differentially expressed in supplemented WT mice
summary(results(dge_dds_APP_cortex_gene_se_10filtered, alpha = 0.05))
## 1 gene was downregulated in supplemented APP mice

## Get a summary of DEGs with FDR (p adj) of <0.05
summary(results(dge_dds_APP_cortex_three_gene_se_10filtered, alpha = 0.05))
## Output is that of the 18837 genes,10 genes were upregulated and 15 were downregulated in supplemented 3 motnh old APP mice
summary(results(dge_dds_WT_cortex_three_gene_se_10filtered, alpha = 0.05))
## Output is that of the 18837 genes,1 genes were upregulated and 2 were downregulated in supplemented 3 motnh old WT mice


summary(results(dge_dds_APP_cortex_six_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19242 genes, 1 gene was upregulated and 4 genes were downregulated in supp mice 6mo APP mice
summary(results(dge_dds_WT_cortex_six_gene_se_10filtered, alpha = 0.05))
## Output is that of the 19242 genes,3 genes were upregulated and 2 were downregulated in supplemented 6 motnh old WT mice


summary(results(dge_dds_APP_cortex_nine_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19225 genes, no genes were differentially expressed in supplemented 9mo APP mice
summary(results(dge_dds_WT_cortex_nine_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19225 genes, no genes were differentially expressed in supplemented 9mo WT mice


summary(results(dge_dds_APP_cortex_twelve_gene_se_10filtered, alpha = 0.05))
##Output is that of the 18741 genes, 1 genes was upregulated and 1 gene was downregulated in supp APP mice at 12 months
summary(results(dge_dds_WT_cortex_twelve_gene_se_10filtered, alpha = 0.05))
##Output is that of the 18741 genes, 2 genes were upregulated in supp WT mice at 12 months


##export DGE Analysis to csv File
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_cortex_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp__APP_DGE_Analysis_padj_sorted_6_13_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_cortex_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_3mo_APP_DGE_Analysis_padj_sorted_6_10_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_cortex_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_6mo_APP_DGE_Analysis_padj_sorted_6_10_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_cortex_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_9mo_APP_DGE_Analysis_padj_sorted_6_10_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_cortex_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_12mo_APP_DGE_Analysis_padj_sorted_6_10_2022.csv")

write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_cortex_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_WT_DGE_Analysis_padj_sorted_6_13_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_cortex_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_3mo_WT_DGE_Analysis_padj_sorted_6_10_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_cortex_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_6mo_WT_DGE_Analysis_padj_sorted_6_10_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_cortex_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_9mo_WT_DGE_Analysis_padj_sorted_6_10_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_cortex_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_12mo_WT_DGE_Analysis_padj_sorted_6_10_2022.csv")


write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_cortex_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_APP_DGE_Analysis_pval_sorted_6_13_2022.csv')
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_cortex_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_3mo_APP_DGE_Analysis_pval_sorted_6_10_2022.csv')
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_cortex_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_6mo_APP_DGE_Analysis_pval_sorted_6_10_2022.csv")
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_cortex_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_9mo_APP_DGE_Analysis_pval_sorted_6_10_2022.csv")
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_APP_cortex_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_12mo_APP_DGE_Analysis_pval_sorted_6_10_2022.csv")

write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_cortex_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_WT_DGE_Analysis_pval_sorted_6_13_2022.csv')
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_cortex_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_3mo_WT_DGE_Analysis_pval_sorted_6_10_2022.csv')
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_cortex_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_6mo_WT_DGE_Analysis_pval_sorted_6_10_2022.csv")
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_cortex_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_9mo_WT_DGE_Analysis_pval_sorted_6_10_2022.csv")
write.csv(as.data.frame(gene_expression_pval_ordered_from_dge_dds_WT_cortex_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_Control_v_Supp_12mo_WT_DGE_Analysis_pval_sorted_6_10_2022.csv")
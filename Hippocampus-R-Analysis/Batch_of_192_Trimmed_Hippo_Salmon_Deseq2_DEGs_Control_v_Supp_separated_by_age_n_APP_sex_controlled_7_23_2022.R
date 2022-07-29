## Batch of 192- Trimmed Hippocampus Data 
## This is the data that was run on Salmon 6-6-2022 and 6-7-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 6/10/2022
## Will import data using tximeta to create a summarized experiment First
## Using https://www.reneshbedre.com/blog/deseq2.html and https://lashlock.github.io/compbio/R_presentation.html as a guide
## Trying to Use DESeq2 to get DEGs for Control vs Supp diet at all time points but separated into just APP and just WT (so a comparison of Control vs Supplemented at 3, 6, 9, and 12 months old)




library(BiocManager)
library(tximport)
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

dir <- setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Trimmed_FASTQ_Aligned_and_Quant_Salmon_6_6_2022")
files <- file.path(dir, sample_metadata$names,  "quant.sf" )
file.exists(files)

## Create dataframe that has the location of the quant files and all sample metadata
hippo_coldata <- data.frame(files, names= sample_metadata$names, APP= sample_metadata$APP, Age= sample_metadata$Age, Diet= sample_metadata$Diet, Sex= sample_metadata$Sex, stringsAsFactors = F)
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
dds_APP_hippo_gene_se <- DESeqDataSet(APP_hippo_gene_se, design = ~ Age + Sex + Diet)
dds_APP_hippo_three_gene_se <- DESeqDataSet(APP_hippo_three_gene_se, design = ~ Sex + Diet)
dds_APP_hippo_six_gene_se <- DESeqDataSet(APP_hippo_six_gene_se, design = ~ Sex + Diet)
dds_APP_hippo_nine_gene_se <- DESeqDataSet(APP_hippo_nine_gene_se, design = ~ Sex + Diet)
dds_APP_hippo_twelve_gene_se <- DESeqDataSet(APP_hippo_twelve_gene_se, design = ~  Sex + Diet)

dds_WT_hippo_gene_se <- DESeqDataSet(WT_hippo_gene_se, design = ~ Age + Sex + Diet)
dds_WT_hippo_three_gene_se <- DESeqDataSet(WT_hippo_three_gene_se, design = ~ Sex + Diet)
dds_WT_hippo_six_gene_se <- DESeqDataSet(WT_hippo_six_gene_se, design = ~ Sex + Diet)
dds_WT_hippo_nine_gene_se <- DESeqDataSet(WT_hippo_nine_gene_se, design = ~ Sex + Diet)
dds_WT_hippo_twelve_gene_se <- DESeqDataSet(WT_hippo_twelve_gene_se, design = ~ Sex + Diet)


dds_hippo_gene_se <- DESeqDataSet(hippo_gene_se, design = ~ Age + Sex + APP + Diet)
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
## 5 genes were upregulated and 9 genes were downrefgulated in supplemented mice
summary(results(dge_dds_hippo_gene_se_10filtered))
## 7 genes were upregulated and 11 genes were downregulated in supplemented mice

## Get a summary of DEGs with FDR (padj) of <0.05
summary(results(dge_dds_APP_hippo_gene_se_10filtered, alpha = 0.05))
## 0 genes were upregulated, 1 genes were upregulated in supp APP mice
summary(results(dge_dds_WT_hippo_gene_se_10filtered, alpha = 0.05))
## 0 gene was upregulated and 2 gene downregulated in supp WT mice
## Get a summary of DEGs with FDR (padj) of <0.05
summary(results(dge_dds_APP_hippo_gene_se_10filtered))
## 0 genes were upregulated, 2 genes were upregulated in supp APP mice
summary(results(dge_dds_WT_hippo_gene_se_10filtered))
## 0 gene was upregulated and 2 gene downregulated in supp WT mice




## 3 Months

## Get a summary of DEGs with FDR (p adj) of <0.05
summary(results(dge_dds_APP_hippo_three_gene_se_10filtered, alpha = 0.05))
## Output is that of the 18837 genes,0 genes were upregulated and 0 were downregulated in supplemented 3 motnh old APP mice
summary(results(dge_dds_WT_hippo_three_gene_se_10filtered, alpha = 0.05))
## Output is that of the 18837 genes,0 genes were upregulated and 1 was downregulated in supplemented 3 motnh old WT mice
## Get a summary of DEGs with FDR (p adj) of <0.05
summary(results(dge_dds_APP_hippo_three_gene_se_10filtered))
## Output is that of the 18837 genes,0 genes were upregulated and 0 were downregulated in supplemented 3 motnh old APP mice
summary(results(dge_dds_WT_hippo_three_gene_se_10filtered))
## Output is that of the 18837 genes,2 genes were upregulated and 1 was downregulated in supplemented 3 motnh old WT mice


## 6 Months

summary(results(dge_dds_APP_hippo_six_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19242 genes, 8 genes were upregulated and 3 genes were downregulated in supp mice 6mo APP mice
summary(results(dge_dds_WT_hippo_six_gene_se_10filtered, alpha = 0.05))
## Output is that of the 19242 genes,0 genes were upregulated and 1 were downregulated in supplemented 6 motnh old WT mice
summary(results(dge_dds_APP_hippo_six_gene_se_10filtered))
##Output is that of the 19242 genes, 15 genes were upregulated and 8 genes were downregulated in supp mice 6mo APP mice
summary(results(dge_dds_WT_hippo_six_gene_se_10filtered))
## Output is that of the 19242 genes,0 genes were upregulated and 1 were downregulated in supplemented 6 motnh old WT mice



## 9 Months

summary(results(dge_dds_APP_hippo_nine_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19225 genes, 0 genes were upregulated and 1 genes were downregulated in supplemented 9mo APP mice
summary(results(dge_dds_WT_hippo_nine_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19225 genes, no genes were differentially expressed in supplemented 9mo WT mice
summary(results(dge_dds_APP_hippo_nine_gene_se_10filtered))
##Output is that of the 19225 genes, 0 genes were upregulated and 1 genes were downregulated in supplemented 9mo APP mice
summary(results(dge_dds_WT_hippo_nine_gene_se_10filtered))
##Output is that of the 19225 genes, no genes were differentially expressed in supplemented 9mo WT mice



## 12 Months

summary(results(dge_dds_APP_hippo_twelve_gene_se_10filtered, alpha = 0.05))
##Output is that of the 18741 genes, 172 genes was upregulated and 151 gene was downregulated in supp APP mice at 12 months
summary(results(dge_dds_WT_hippo_twelve_gene_se_10filtered, alpha = 0.05))
##Output is that of the 18741 genes, 0 genes were upregulated and 0 were downregulated in supp WT mice at 12 months
summary(results(dge_dds_APP_hippo_twelve_gene_se_10filtered))
##Output is that of the 18741 genes, 262 genes was upregulated and 259 gene was downregulated in supp APP mice at 12 months
summary(results(dge_dds_WT_hippo_twelve_gene_se_10filtered))
##Output is that of the 18741 genes, 1 genes were upregulated and 1 were downregulated in supp WT mice at 12 months



##export DGE Analysis to csv File
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_APP_DGE_Analysis_Age_sex_controlled_padj_sorted_7_23_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_3mo_APP_DGE_Analysis_sex_controlled_padj_sorted_7_23_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_6mo_APP_DGE_Analysis_sex_controlled_padj_sorted_7_23_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_9mo_APP_DGE_Analysis_sex_controlled_padj_sorted_7_23_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_APP_hippo_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_12mo_APP_DGE_Analysis_sex_controlled_padj_sorted_7_23_2022.csv")

write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_WT_DGE_Analysis_age_sex_controlled_padj_sorted_7_23_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_3mo_WT_DGE_Analysis_sex_controlled_padj_sorted_7_23_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_6mo_WT_DGE_Analysis_sex_controlled_padj_sorted_7_23_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_9mo_WT_DGE_Analysis_sex_controlled_padj_sorted_7_23_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_WT_hippo_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_12mo_WT_DGE_Analysis_sex_controlled_padj_sorted_7_23_2022.csv")

write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_hippo_gene_se_10filtered), file = "Batch_of_192_Trimmed_Hippo_DESeq2_Control_v_Supp_All_animals_DGE_Analysis_sex_age_app_controlled_padj_sorted_7_23_2022.csv")

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

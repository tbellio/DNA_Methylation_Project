## Batch of 192- Trimmed Cortex Data 
## This is the data that was run on Salmon 6-6-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 6/8/2022
## Will import data using tximeta to create a summarized experiment First
## Using https://www.reneshbedre.com/blog/deseq2.html and https://lashlock.github.io/compbio/R_presentation.html as a guide
## Trying to Use DESeq2 to get DEGs for WT vs APP at all time points (so a comparison of WT vs APP at 3, 6, 9, and 12 months old)




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

## Now we have the information for tximeta, lets load it in to create a summarized experiment with all the sample metadata
## This should also find all the matching transcriptome information-- Ensembl Mus musculus release 104
cortex__three_se <- tximeta(cortex_three_only_coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe
cortex_six_se <- tximeta(cortex_six_only_coldata)
cortex_nine_se <- tximeta(cortex_nine_only_coldata)
cortex_twelve_se <- tximeta(cortex_twelve_only_coldata)


suppressPackageStartupMessages(library(SummarizedExperiment))
## check the column data to see if everything loaded properly; should have all the names, APP, Age, Diet, and Sex information
colData(cortex__three_se)
colData(cortex_six_se)
colData(cortex_nine_se)
colData(cortex_twelve_se)
## Check the assay names; should be 3- counts, abundance, length
assayNames(cortex__three_se)
assayNames(cortex_six_se)
assayNames(cortex_nine_se)
assayNames(cortex_twelve_se)
## Check the row ranges
rowRanges(cortex__three_se)
rowRanges(cortex_six_se)
rowRanges(cortex_nine_se)
rowRanges(cortex_twelve_se)
## Check the seqinfo
seqinfo(cortex__three_se)
seqinfo(cortex_six_se)
seqinfo(cortex_nine_se)
seqinfo(cortex_twelve_se)
## get the database that was used to match the Salmon input quants
edb <- retrieveDb(cortex__three_se)
class(edb)

##adding exons to the summarized experiment (hippo__three_se, etc.)
cortex_three_se.exons <- addExons(cortex__three_se)
cortex_six_se.exons <- addExons(cortex_six_se)
cortex_nine_se.exons <- addExons(cortex_nine_se)
cortex_twelve_se.exons <- addExons(cortex_twelve_se)

##Check rowranges--shouldnt have to make sure it works for all four as same code
rowRanges(cortex_three_se.exons)

## Summarize read counts to the gene level
cortex_three_gene_se <- summarizeToGene(cortex_three_se.exons)
cortex_six_gene_se <- summarizeToGene(cortex_six_se.exons)
cortex_nine_gene_se <- summarizeToGene(cortex_nine_se.exons)
cortex_twelve_gene_se <- summarizeToGene(cortex_twelve_se.exons)

rowRanges(cortex_three_gene_se)


## Now use DESeq2 to look for DEGs at each time point
library('DESeq2')
dds_cortex_three_gene_se <- DESeqDataSet(cortex_three_gene_se, design = ~ APP)
dds_cortex_six_gene_se <- DESeqDataSet(cortex_six_gene_se, design = ~ APP)
dds_cortex_nine_gene_se <- DESeqDataSet(cortex_nine_gene_se, design = ~ APP)
dds_cortex_twelve_gene_se <- DESeqDataSet(cortex_twelve_gene_se, design = ~ APP)
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; if not it will change it to a factor (this code has it as a character so it knows to change it)

## Filter so that only genes that have more than 10 reads are inculuded
dds_cortex_three_gene_se_10filtered <- dds_cortex_three_gene_se[rowSums(counts(dds_cortex_three_gene_se)) >= 10,]
dds_cortex_six_gene_se_10filtered <- dds_cortex_six_gene_se[rowSums(counts(dds_cortex_six_gene_se)) >= 10,]
dds_cortex_nine_gene_se_10filtered <- dds_cortex_nine_gene_se[rowSums(counts(dds_cortex_nine_gene_se)) >= 10,]
dds_cortex_twelve_gene_se_10filtered <- dds_cortex_twelve_gene_se[rowSums(counts(dds_cortex_twelve_gene_se)) >= 10,]
## this filtered out ~18000 genes; went from 35682 to 18837, 19242, 19225, 19428, respectively

## Since DESeq uses alphabetical as reference when comparing want to change the reference to the wildtype mice
dds_cortex_three_gene_se_10filtered$APP <- relevel(dds_cortex_three_gene_se_10filtered$APP, ref = "WT")
dds_cortex_six_gene_se_10filtered$APP <- relevel(dds_cortex_six_gene_se_10filtered$APP, ref = "WT")
dds_cortex_nine_gene_se_10filtered$APP <- relevel(dds_cortex_nine_gene_se_10filtered$APP, ref = "WT")
dds_cortex_twelve_gene_se_10filtered$APP <- relevel(dds_cortex_twelve_gene_se_10filtered$APP, ref = "WT")
## Run DGE analysis
dge_dds_cortex_three_gene_se_10filtered <- DESeq(dds_cortex_three_gene_se_10filtered)
dge_dds_cortex_six_gene_se_10filtered <- DESeq(dds_cortex_six_gene_se_10filtered)
dge_dds_cortex_nine_gene_se_10filtered <- DESeq(dds_cortex_nine_gene_se_10filtered)
dge_dds_cortex_twelve_gene_se_10filtered <- DESeq(dds_cortex_twelve_gene_se_10filtered)

## See which comparisons were run; should only be one as model was just APP vs WT
resultsNames(dge_dds_cortex_three_gene_se_10filtered)
resultsNames(dge_dds_cortex_six_gene_se_10filtered)
resultsNames(dge_dds_cortex_nine_gene_se_10filtered)
resultsNames(dge_dds_cortex_twelve_gene_se_10filtered)

## Get gene expression table as DESeqResults
gene_expression_from_dge_dds_cortex_three_gene_se_10filtered <- results(dge_dds_cortex_three_gene_se_10filtered)
gene_expression_from_dge_dds_cortex_six_gene_se_10filtered <- results(dge_dds_cortex_six_gene_se_10filtered)
gene_expression_from_dge_dds_cortex_nine_gene_se_10filtered <- results(dge_dds_cortex_nine_gene_se_10filtered)
gene_expression_from_dge_dds_cortex_twelve_gene_se_10filtered <- results(dge_dds_cortex_twelve_gene_se_10filtered)
## See what the gene expression table says
## Will tell you what it was comparing (so APP vs WT for this) and then the log2 fold change, log fold change, stat, pvalue, and adjusted p value (this is Benjamini-Hochberg FDR method)
gene_expression_from_dge_dds_cortex_three_gene_se_10filtered
gene_expression_from_dge_dds_cortex_six_gene_se_10filtered
gene_expression_from_dge_dds_cortex_nine_gene_se_10filtered
gene_expression_from_dge_dds_cortex_twelve_gene_se_10filtered

## Sort the gene expression data by FDR(adjusted p-value)
gene_expression_padj_ordered_from_dge_dds_cortex_three_gene_se_10filtered <- gene_expression_from_dge_dds_cortex_three_gene_se_10filtered[order(gene_expression_from_dge_dds_cortex_three_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_cortex_six_gene_se_10filtered <- gene_expression_from_dge_dds_cortex_six_gene_se_10filtered[order(gene_expression_from_dge_dds_cortex_six_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_cortex_nine_gene_se_10filtered <- gene_expression_from_dge_dds_cortex_nine_gene_se_10filtered[order(gene_expression_from_dge_dds_cortex_nine_gene_se_10filtered$padj),]
gene_expression_padj_ordered_from_dge_dds_cortex_twelve_gene_se_10filtered <- gene_expression_from_dge_dds_cortex_twelve_gene_se_10filtered[order(gene_expression_from_dge_dds_cortex_twelve_gene_se_10filtered$padj),]


## Get a summary of DEGs with FDR (p adj) of <0.05
summary(results(dge_dds_cortex_three_gene_se_10filtered, alpha = 0.05))
## Output is that of the 18837 genes, 171 genes were upregulated in APP and 13 was downregulated at 3 months
summary(results(dge_dds_cortex_six_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19242 genes, 1075 genes were upregulated and 338 were downregulated in APP mice at 6 months
summary(results(dge_dds_cortex_nine_gene_se_10filtered, alpha = 0.05))
##Output is that of the 19225 genes, 1730 genes were upregulated and 932 were downregulated in APP mice at 9 months
summary(results(dge_dds_cortex_twelve_gene_se_10filtered, alpha = 0.05))
##Output is that of the 18741 genes, 2037 genes were upregulated and 1741 were downregulated in APP mice at 12 months

##export DGE Analysis to csv File
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis/Trimmed_Cortex_R_Analysis_Starting_6_8_2022/")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_cortex_three_gene_se_10filtered), file = 'Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_3mo_DGE_Analysis_padj_sorted_6_8_2022.csv')
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_cortex_six_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_6mo_DGE_Analysis_padj_sorted_6_8_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_cortex_nine_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_9mo_DGE_Analysis_padj_sorted_6_8_2022.csv")
write.csv(as.data.frame(gene_expression_padj_ordered_from_dge_dds_cortex_twelve_gene_se_10filtered), file = "Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_12mo_DGE_Analysis_padj_sorted_6_8_2022.csv")

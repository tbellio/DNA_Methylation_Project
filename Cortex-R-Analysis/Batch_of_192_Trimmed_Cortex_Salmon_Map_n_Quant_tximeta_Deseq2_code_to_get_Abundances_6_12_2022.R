## Batch of 192- Trimmed Cortex Data
## This is the data that was run on Salmon 6-6-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 6/12/2022
##trying to use tximeta to create abundances table with gene symbol

library(BiocManager)
library(tximport)
library(tximportData)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192")
samples <- read_excel("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Batch_of_192_metadata.xlsx")  #reads in excel file with data about each sample
samples
samples <- as.data.frame(samples)  ##turns the samples into a data frame
samples
class(samples$Sample)  #3checking class of the samples column in the samples data frame
class(samples$Sex)
samples <- transform(samples, Age = as.character(Age))  ##changes the class of the Age column in samples to character
class(samples$Age)
rownames(samples) <- samples$Sample   ##turns the row names of samples data frame into what is said in the sample column
samples$names <- samples$Sample
samples

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
dds_cortex_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ APP)
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_gene_se_10filtered <- dds_cortex_gene_se[rowSums(counts(dds_cortex_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 20480

##create abundances data frame
filtered_cortex_gene_abundances <- as.data.frame(dds_cortex_gene_se_10filtered@assays@data@listData[["abundance"]])

##want to change the names of the rows from the ENSEMBL gene names to the Symbol
##first need to create annotation that matches the rownames of filtered_gene_abundances data frame with the gene symbol
gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, row.names(filtered_cortex_gene_abundances), "SYMBOL", keytype = "ENSEMBL")


##this creates gene_symbols matrix that has 2 columns (ensembl gene name and gene symbol); however there are more genes in this annotation than we have in our filtered_gene_abundances matrix
length(unique(gene_symbols$ENSEMBL))  #21234 unique emsembl gene names
length(unique(gene_symbols$SYMBOL))  ##19372 unique gene symbols
## so we need to adjust gene_symbols so we dont have any overlapping genes
gene_symbols[duplicated(gene_symbols$ENSEMBL),]  #gives you the duplicated ensembl gene names
##use dplyr to select just the genes that are not duplicated and call it nonduplicated_gene_symbols
library(dplyr)
nonduplicated_gene_symbols <- distinct(gene_symbols, ENSEMBL, .keep_all = T)
##now make sure that the nonduplicated_gene_symbols is the same length as your filtered_gene_abundances
nrow(filtered_cortex_gene_abundances)
nrow(nonduplicated_gene_symbols)
## or 
nrow(filtered_cortex_gene_abundances) ==  nrow(nonduplicated_gene_symbols) ##make sure this is TRUE

##now you need to make sure that the row names of your abundance data set matches the data in your annotation
test <- row.names(filtered_cortex_gene_abundances) %in% nonduplicated_gene_symbols$ENSEMBL
sum(test) #make sure that what this returns matches your number of rows in both your abundance data set and your annotation (gene symbols)

##now you can set the row.names of your filtered_gene_abundancs to the gene symbols from nonduplicated_gene_symbols since you know they match
row.names(filtered_cortex_gene_abundances) <- nonduplicated_gene_symbols$SYMBOL

##will have to create column within the filtered_gene_abundances that has gene symbols
filtered_cortex_gene_abundances$symbol <- nonduplicated_gene_symbols$SYMBOL
dup_symbols <- which(duplicated(filtered_cortex_gene_abundances$symbol))
sum(dup_symbols) ##get over 1947 duplicate gene names
## wait just kidding, maybe if i transform filtered_gene_abundances to a matrix I can rename the columns and have duplicates
## but want to make sure I still have the Ensembl name in this so make a new column in the dataframe for the esnsembl name before converting to a matrix
filtered_cortex_gene_abundances$ENSEMBL <- row.names(filtered_cortex_gene_abundances)
filtered_cortex_gene_abundances_matrix <- as.matrix(filtered_cortex_gene_abundances)
row.names(filtered_cortex_gene_abundances_matrix) <- filtered_cortex_gene_abundances$symbol


## export to csv
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis/")
write.csv(filtered_cortex_gene_abundances_matrix,
          file = "Gene_Abundances(TPM)_for_filtered_genes_with_gene_symbols_Batch_of_192_Trimmed_Cortex_Salmon_map_and_quant_6_12_2022.csv")

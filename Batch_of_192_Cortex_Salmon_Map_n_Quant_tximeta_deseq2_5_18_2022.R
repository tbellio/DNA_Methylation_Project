##trying to use tximeta to create abundances table with gene symbol

library(BiocManager)
library(tximport)
library(tximportData)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

getwd()
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

dir <- setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_aligned_Against_Mouse_cdna_Output_5_2022")
getwd()
files <- file.path(dir, samples$names,  "quant.sf" )
file.exists(files)

coldata <- data.frame(files, names= samples$names, APP= samples$APP, Age= samples$Age, Diet= samples$Diet, Sex= samples$Sex, stringsAsFactors = F)
coldata

se <- tximeta(coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe


suppressPackageStartupMessages(library(SummarizedExperiment))
colData(se)

assayNames(se)
rowRanges(se)
seqinfo(se)


edb <- retrieveDb(se)
class(edb)

##adding exons to the summarized experiment (se)
se.exons <- addExons(se)

rowRanges(se.exons)


gse <- summarizeToGene(se)
rowRanges(gse)




suppressPackageStartupMessages(library(DESeq2))
suppressWarnings({dds <- DESeqDataSet(gse, design = ~ APP)})
dds

##Filter the dds DESeqDataSet so that there are only reads with more than 1 read
keep <- rowSums(counts(dds)) > 1
filtered_dds <- dds[keep,]  
##This reduced total genes from 35682 to 24642

##create abundances data frame
filtered_gene_abundances <- as.data.frame(filtered_dds@assays@data@listData[["abundance"]])

##want to change the names of the rows from the ENSEMBL gene names to the Symbol
##first need to create annotation that matches the rownames of filtered_gene_abundances data frame with the gene symbol
gene_symbols <- select(org.Mm.eg.db, row.names(filtered_gene_abundances), "SYMBOL", keytype = "ENSEMBL")


##this creates gene_symbols matrix that has 2 columns (ensembl gene name and gene symbol); however there are more genes in this annotation than we have in our filtered_gene_abundances matrix
length(unique(gene_symbols$ENSEMBL))  #18681 unique emsembl gene names
length(unique(gene_symbols$SYMBOL))  ##17593 unique gene symbols
## so we need to adjust gene_symbols so we dont have any overlapping genes
gene_symbols[duplicated(gene_symbols$ENSEMBL),]  #gives you the duplicated ensembl gene names
##use dplyr to select just the genes that are not duplicated and call it nonduplicated_gene_symbols
library(dplyr)
nonduplicated_gene_symbols <- distinct(gene_symbols, ENSEMBL, .keep_all = T)
##now make sure that the nonduplicated_gene_symbols is the same length as your filtered_gene_abundances
nrow(filtered_gene_abundances)
nrow(nonduplicated_gene_symbols)
## or 
nrow(filtered_gene_abundances) ==  nrow(nonduplicated_gene_symbols) ##make sure this is TRUE

##now you need to make sure that the row names of your abundance data set matches the data in your annotation
test <- row.names(filtered_gene_abundances) %in% nonduplicated_gene_symbols$ENSEMBL
sum(test) #make sure that what this returns matches your number of rows in both your abundance data set and your annotation (gene symbols)

##now you can set the row.names of your filtered_gene_abundancs to the gene symbols from nonduplicated_gene_symbols since you know they match
row.names(filtered_gene_abundances) <- nonduplicated_gene_symbols$SYMBOL
##will have to create column within the filtered_gene_abundances that has gene symbols
filtered_gene_abundances$symbol <- nonduplicated_gene_symbols$SYMBOL
dup_symbols <- which(duplicated(filtered_gene_abundances$symbol))
sum(dup_symbols) ##get over 3418 duplicate gene names
## wait just kidding, maybe if i transform filtered_gene_abundances to a matrix I can rename the columns and have duplicates
## but want to make sure I still have the Ensembl name in this so make a new column in the dataframe for the esnsembl name before converting to a matrix
filtered_gene_abundances$ENSEMBL <- row.names(filtered_gene_abundances)
filtered_gene_abundances_matrix <- as.matrix(filtered_gene_abundances)
row.names(filtered_gene_abundances_matrix) <- filtered_gene_abundances$symbol
##delete the symbol column as now rows are names with the symbol
filtered_gene_abundances_matrix <- filtered_gene_abundances_matrix[,-97]

## now transpose the filtered gene abundances matrix so you can export it to excel for jmp file
transposed_filtered_gene_abundances <- t(filtered_gene_abundances_matrix)

## export to csv
write.csv(transposed_filtered_gene_abundances,
          file = "Gene_Abundances(TPM)_for_filtered_genes_with_gene_symbols_as_column_names_Batch_of_192_Cortex_Salmon_map_and_quant_5_18_2022.csv")

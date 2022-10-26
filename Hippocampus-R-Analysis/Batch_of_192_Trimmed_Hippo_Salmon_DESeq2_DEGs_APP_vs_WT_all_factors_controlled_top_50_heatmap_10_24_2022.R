## Batch of 192- Trimmed Hippocampus Data
## This is the data that was run on Salmon 6-6-2022 and 6-7-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 6/22/2022 for DEGs
## Using DESeq2 to get DEGs for WT vs APP
## Trying to make heatmap from top 50 DEGs
## Will import data using tximeta to create a summarized experiment First
## https://igordot.github.io/tutorials/heatmaps-2017-07.nb.html as a guide


library(BiocManager)
library(tximport)
library(tximportData)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)
library(ggplot2)
library('DESeq2')
library(pheatmap)
library(dendextend)
library(RColorBrewer)

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

## Now use DESeq2 to look for DEGs
library('DESeq2')
dds_hippo_gene_se <- DESeqDataSet(hippo_gene_se, design = ~ Age + Sex + Diet + APP)
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_hippo_gene_se_10filtered <- dds_hippo_gene_se[rowSums(counts(dds_hippo_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 20480

##create abundances data frame
filtered_hippo_gene_abundances <- as.data.frame(dds_hippo_gene_se_10filtered@assays@data@listData[["abundance"]])

## turn cortex gene abundances into matrix
filtered_hippo_gene_abundances_matrix <- as.matrix(filtered_hippo_gene_abundances)

## Subset filtered cortex abund so only genes that have more than 500000 total reads
filtered_hippo_gene_abundances_matrix_subset <- filtered_hippo_gene_abundances_matrix[rowSums(filtered_hippo_gene_abundances_matrix)>500000,]
head(filtered_hippo_gene_abundances_matrix_subset)

## Create a heatmap of genes that have more than 500000 total reads
pheatmap(filtered_hippo_gene_abundances_matrix_subset)

## can now scale the rows via z-score
cal_z_score <- function(x) {
  (x- mean(x)) / sd(x)
}
filtered_hippo_gene_abundances_matrix_subset_zscore <- t(apply(filtered_hippo_gene_abundances_matrix_subset, 1, cal_z_score))
pheatmap(filtered_hippo_gene_abundances_matrix_subset_zscore)

## Ok so now I know that I can create a heatmap
## Lets try to get only the top 100 DEGs for APP vs WT in this heatmap
## So lets get the top 50 genes
## Need to do actual DEG analysis with DESeq2
## First need to make sure WT used as the control
dds_hippo_gene_se_10filtered$APP <- relevel(dds_hippo_gene_se_10filtered$APP, ref = "WT")

## Run DEG analysis with DESeq2
dge_dds_hippo_gene_se_10filtered <- DESeq(dds_hippo_gene_se_10filtered) 

## Look at results from DEG analysis
results_dge_dds_hippo_gene_se_10filtered <- results(dge_dds_hippo_gene_se_10filtered)
results_dge_dds_hippo_gene_se_10filtered

## Sort the results of DEG analysis by p-adj
padj_sorted_results_dge_dds_hippo_gene_se_10filtered <- results_dge_dds_hippo_gene_se_10filtered[order(results_dge_dds_hippo_gene_se_10filtered$padj),]
padj_sorted_results_dge_dds_hippo_gene_se_10filtered

## Create list of top 50 DEGs
top50 <-row.names(padj_sorted_results_dge_dds_hippo_gene_se_10filtered[1:50,])



## Making the rownames gene symbols instead of ENSEMBL names
gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, row.names(filtered_hippo_gene_abundances_matrix), "SYMBOL", keytype = "ENSEMBL")
##this creates gene_symbols matrix that has 2 columns (ensembl gene name and gene symbol); however there are more genes in this annotation than we have in our filtered_gene_abundances matrix
length(unique(gene_symbols$ENSEMBL))  #20480 unique emsembl gene names
length(unique(gene_symbols$SYMBOL))  ##18824 unique gene symbols
## so we need to adjust gene_symbols so we dont have any overlapping genes
gene_symbols[duplicated(gene_symbols$ENSEMBL),]  #gives you the duplicated ensembl gene names
##use dplyr to select just the genes that are not duplicated and call it nonduplicated_gene_symbols
library(dplyr)
nonduplicated_gene_symbols <- distinct(gene_symbols, ENSEMBL, .keep_all = T)
##now make sure that the nonduplicated_gene_symbols is the same length as your filtered_gene_abundances
nrow(filtered_hippo_gene_abundances_matrix)
nrow(nonduplicated_gene_symbols)
## or 
nrow(filtered_hippo_gene_abundances_matrix) ==  nrow(nonduplicated_gene_symbols) ##make sure this is TRUE

##now you need to make sure that the row names of your abundance data set matches the data in your annotation
test <- row.names(filtered_hippo_gene_abundances_matrix) %in% nonduplicated_gene_symbols$ENSEMBL
sum(test) #make sure that what this returns matches your number of rows in both your abundance data set and your annotation (gene symbols)

##now you can set the row.names of your filtered_gene_abundancs to the gene symbols from nonduplicated_gene_symbols since you know they match
row.names(filtered_hippo_gene_abundances_matrix) <- nonduplicated_gene_symbols$SYMBOL
## Check to see if row names have been renamed
filtered_hippo_gene_abundances_matrix ##they have


## Switch top 100 and top 50 lists from ENSEMBL names to Gene symbols
top50 <- as.data.frame(top50)
top50
top50_gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, top50$top50, "SYMBOL", keytype = "ENSEMBL")
row.names(top50_gene_symbols) <- top50_gene_symbols$SYMBOL
top50_gene_symbols
top50 <- row.names(top50_gene_symbols)
top50

## now need to subset the filtered_cortex_gene_abundances_matrix so that we only have the top 100 genes remaining
top50_filtered_hippo <- filtered_hippo_gene_abundances_matrix[rownames(filtered_hippo_gene_abundances_matrix) %in% top50,]
top50_filtered_hippo


## make heatmap
pheatmap(top50_filtered_hippo)
## this looks awful because things aren't normalized so lets normalize
## make function to z-score
cal_z_score <- function(x) {
  (x- mean(x)) / sd(x)
}

## apply z-score function to top 100 gene abundances
top50_filtered_hippo_gene_abundances_matrix_z_score <- t(apply(top50_filtered_hippo,1, cal_z_score))

## Make heatmap of top 100 DEGs APP vs WT (all ages) in z-score
pheatmap(top50_filtered_hippo_gene_abundances_matrix_z_score)


## Make annotation information with metadata
my_sample_col <- data.frame(APP = sample_metadata$APP, Diet= sample_metadata$Diet, Age= sample_metadata$Age)
rownames(my_sample_col) <- sample_metadata$Sample
my_sample_col

## Test heatmap
pheatmap(top50_filtered_hippo_gene_abundances_matrix_z_score, annotation_col = my_sample_col)


## Okay lets change colors
annot_colors_gradient= list(APP=c(APP="red1", WT= "royalblue1"), Diet= c(Control= "gold", Supplemented= "goldenrod3"), Age= c('3'='gray99', '6'='gray75', '9'= 'gray50', '12'= 'gray25'))
hippo_top50_heatmap <- pheatmap(top50_filtered_hippo_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors_gradient,  cluster_rows = F, color = brewer.pal(11, "PRGn"), border_color = "gray60")

## Get dendrogram of clustering so can figure out how to rotate branch
my_hclust_sample_top50_hippo <- hclust(dist(t(top50_filtered_hippo_gene_abundances_matrix_z_score)), method= "complete")
as.dendrogram(my_hclust_sample_top50_hippo) %>%
  plot()

## Rotate branches by adjusting order
col_dend <- hippo_top50_heatmap[[2]]
col_dend
col_dend <- dendextend::rotate(col_dend, order= c( '46', '103', '129', '202', '203', '161', '213', '130', '189', 
                                                   '31', '155', '67', '141', '72', '134', '33', '76',
                                                   '164', '79', '24', '25', '139', '57', '190', '152', '124', '125',
                                                   
                                                   '186', '58', '211', '149', '127', '216', '217', '183',
                                                   '9', '8', 
                                                   '236', '20', '226',
                                                   '204', '27', '36', 
                                                   '174', '61', '96', '74', 
                                                   '157', '106', 
                                                   '95', 
                                                   '287', '265', '268', '251', '249', '122',
                                                   '178', '242', '185', '18',
                                                   
                                                   '84', '41', '40', '102', '132', '90', '89', '194', '191', '252', '290', '192',
                                                   
                                                   '275', '225', '147', '117', '153', '228', '282', '160', '199', 
                                                   '262', '244', '168', '107', '263', '221', '207', '285', '272', '283', '277', '259', '223', '238',
                                                   
                                                   
                                                   '253'
                                                   ))
hippo_top50_heatmap <- pheatmap(top50_filtered_hippo_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors_gradient,  cluster_cols = as.hclust(col_dend),  cluster_rows = F, color = brewer.pal(11, "PRGn"), border_color = "gray60")
hippo_top50_heatmap                               

hippo_top50_heatmap_row_cluster <- pheatmap(top50_filtered_hippo_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors_gradient,  cluster_cols = as.hclust(col_dend),  cluster_rows = T, color = brewer.pal(11, "PRGn"), border_color = "black")


## Trying to adjust scale so it is the same as the cortex heatmap
## Using https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r as a guide
## Create vector called myColor that has pallette
myColor <- brewer.pal(11,"PRGn")
myColor
myColor2 <- c( "#762A83", "#9970AB", "#C2A5CF", '#E7D4E8', "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B")
# this shows you what colors are in this palette
## palette length is max of 11 with PRGn with brewer.pal() function
paletteLength <- 11

## length(breaks) == length(palettelength) +1
myBreaks <- c(-1.33,-1,-.66,-0.33,0,.75,1.5,2.25, 3)
myBreaks2 <- c(-1.33,-1,-.66,-0.33,0, 0.33, 0.66, 1, 1.33, 1.66, 2)
myBreaks3 <- c(-1.5,-1,-.5,0,.5,1,1.5,2, 2.5, 3)

hippo_top50_heatmap_scaled <- pheatmap(top50_filtered_hippo_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors_gradient,  cluster_cols = as.hclust(col_dend),  cluster_rows = F, 
                                       color = myColor, border_color = "gray60", breaks = myBreaks2)


## Ok lets try this again but using their methodology
## 10-26-2022
paletteLength <- 11
myColor_ex <- brewer.pal(11, "PRGn")
myColor_ex
palettelength2 <- 50
myColor_ex2 <- colorRampPalette(c("#762A83", "#9970AB", "#C2A5CF", '#E7D4E8', "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B"))(palettelength2)
myColor_ex2
## length(breaks) == length(palettelength) +1
## use floor and ceiling to deal with even/odd length pallettelenghts
myBreaks <- c(seq(min(top50_filtered_hippo_gene_abundances_matrix_z_score), 0, length.out=ceiling(palettelength2/2)+1),
              seq(max(top50_filtered_hippo_gene_abundances_matrix_z_score)/palettelength2, max(top50_filtered_hippo_gene_abundances_matrix_z_score), length.out=floor(palettelength2/2)))

hippo_top50_heatmap_overflow <- pheatmap(top50_filtered_hippo_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors = annot_colors_gradient, cluster_cols = as.hclust(col_dend), cluster_rows = F,
                                         color = myColor_ex2, border_color = "gray60", breaks = myBreaks,
                                         main = "Hippocampus Hierarchical Clustering of Top 50 DEGs",
                                         fontsize = 24,
                                         fontsize_row =  10,
                                         fontsize_col= 10)

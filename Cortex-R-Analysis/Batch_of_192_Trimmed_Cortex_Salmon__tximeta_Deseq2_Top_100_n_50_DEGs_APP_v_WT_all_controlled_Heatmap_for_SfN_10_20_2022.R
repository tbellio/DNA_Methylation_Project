## Batch of 192- Trimmed Cortex Data 
## This is the data that was run on Salmon 6-6-2022 with the Trimmed reads from Trimmomatic (used 1P and 2P.fastq files for salmon)
## Starting this on 6/8/2022
## Will import data using tximeta to create a summarized experiment First
## Using  https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/ as a guide
## Trying to visualize DEGs
## Have very nice Top 50 heatmap that is ready for poster as of 10-24-2022




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

## Now we have the information for tximeta, lets load it in to create a summarized experiment with all the sample metadata
## This should also find all the matching transcriptome information-- Ensembl Mus musculus release 104
cortex_se <- tximeta(cortex_coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe


suppressPackageStartupMessages(library(SummarizedExperiment))
## check the column data to see if everything loaded properly; should have all the names, APP, Age, Diet, and Sex information
colData(cortex_se)

##adding exons to the summarized experiment (hippo_se)
cortex_se.exons <- addExons(cortex_se)

rowRanges(cortex_se.exons)

## Summarize read counts to the gene level
cortex_gene_se <- summarizeToGene(cortex_se)
rowRanges(cortex_gene_se)

## Now use DESeq2 to look for DEGs
library('DESeq2')
dds_cortex_gene_se <- DESeqDataSet(cortex_gene_se, design = ~ Age + Sex + Diet + APP)
## Check to make sure that APP in the dds_hippo_gene_se DESeqDataSet is a factor; it will automatically convert if not

## Filter so that only genes that have more than 10 reads total are inculuded
dds_cortex_gene_se_10filtered <- dds_cortex_gene_se[rowSums(counts(dds_cortex_gene_se)) >= 10,]
## this filtered out ~15000 genes; went from 35682 to 20480

##create abundances data frame
filtered_cortex_gene_abundances <- as.data.frame(dds_cortex_gene_se_10filtered@assays@data@listData[["abundance"]])

## turn cortex gene abundances into matrix
filtered_cortex_gene_abundances_matrix <- as.matrix(filtered_cortex_gene_abundances)

## Subset filtered cortex abund so only genes that have more than 500000 total reads
filtered_cortex_gene_abundances_matrix_subset <- filtered_cortex_gene_abundances_matrix[rowSums(filtered_cortex_gene_abundances_matrix)>500000,]
head(filtered_cortex_gene_abundances_matrix_subset)

## Create a heatmap of genes that have more than 500000 total reads
pheatmap(filtered_cortex_gene_abundances_matrix_subset)

## can now scale the rows via z-score
cal_z_score <- function(x) {
        (x- mean(x)) / sd(x)
    }
filtered_cortex_gene_abundances_matrix_subset_zscore <- t(apply(filtered_cortex_gene_abundances_matrix_subset, 1, cal_z_score))
pheatmap(filtered_cortex_gene_abundances_matrix_subset_zscore)

## Ok so now I know that I can create a heatmap
## Lets try to get only the top 100 DEGs for APP vs WT in this heatmap
## So lets get the top 100 genes
## Need to do actual DEG analysis with DESeq2
## First need to make sure WT used as the control
dds_cortex_gene_se_10filtered$APP <- relevel(dds_cortex_gene_se_10filtered$APP, ref = "WT")

## Run DEG analysis with DESeq2
dge_dds_cortex_gene_se_10filtered <- DESeq(dds_cortex_gene_se_10filtered) 

## Look at results from DEG analysis
results_dge_dds_cortex_gene_se_10filtered <- results(dge_dds_cortex_gene_se_10filtered)
results_dge_dds_cortex_gene_se_10filtered

## Sort the results of DEG analysis by p-adj
padj_sorted_results_dge_dds_cortex_gene_se_10filtered <- results_dge_dds_cortex_gene_se_10filtered[order(results_dge_dds_cortex_gene_se_10filtered$padj),]
padj_sorted_results_dge_dds_cortex_gene_se_10filtered

## Create list of top 100 DEGs
top100 <- row.names(padj_sorted_results_dge_dds_cortex_gene_se_10filtered[1:100,])
top50 <-row.names(padj_sorted_results_dge_dds_cortex_gene_se_10filtered[1:50,])



## Making the rownames gene symbols instead of ENSEMBL names
gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, row.names(filtered_cortex_gene_abundances_matrix), "SYMBOL", keytype = "ENSEMBL")
##this creates gene_symbols matrix that has 2 columns (ensembl gene name and gene symbol); however there are more genes in this annotation than we have in our filtered_gene_abundances matrix
length(unique(gene_symbols$ENSEMBL))  #21234 unique emsembl gene names
length(unique(gene_symbols$SYMBOL))  ##19373 unique gene symbols
## so we need to adjust gene_symbols so we dont have any overlapping genes
gene_symbols[duplicated(gene_symbols$ENSEMBL),]  #gives you the duplicated ensembl gene names
##use dplyr to select just the genes that are not duplicated and call it nonduplicated_gene_symbols
library(dplyr)
nonduplicated_gene_symbols <- distinct(gene_symbols, ENSEMBL, .keep_all = T)
##now make sure that the nonduplicated_gene_symbols is the same length as your filtered_gene_abundances
nrow(filtered_cortex_gene_abundances_matrix)
nrow(nonduplicated_gene_symbols)
## or 
nrow(filtered_cortex_gene_abundances_matrix) ==  nrow(nonduplicated_gene_symbols) ##make sure this is TRUE

##now you need to make sure that the row names of your abundance data set matches the data in your annotation
test <- row.names(filtered_cortex_gene_abundances_matrix) %in% nonduplicated_gene_symbols$ENSEMBL
sum(test) #make sure that what this returns matches your number of rows in both your abundance data set and your annotation (gene symbols)

##now you can set the row.names of your filtered_gene_abundancs to the gene symbols from nonduplicated_gene_symbols since you know they match
row.names(filtered_cortex_gene_abundances_matrix) <- nonduplicated_gene_symbols$SYMBOL
## Check to see if row names have been renamed
filtered_cortex_gene_abundances_matrix ##they have


## Switch top 100 and top 50 lists from ENSEMBL names to Gene symbols
top100 <- as.data.frame(top100)
top100
top50 <- as.data.frame(top50)
top50
top100_gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, top100$top100, "SYMBOL", keytype = "ENSEMBL")
length(unique(top100_gene_symbols$SYMBOL))
row.names(top100_gene_symbols) <- top100_gene_symbols$SYMBOL
top100_gene_symbols
top50_gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, top50$top50, "SYMBOL", keytype = "ENSEMBL")
row.names(top50_gene_symbols) <- top50_gene_symbols$SYMBOL
top50_gene_symbols
top100 <- row.names(top100_gene_symbols)
top50 <- row.names(top50_gene_symbols)
top50

## now need to subset the filtered_cortex_gene_abundances_matrix so that we only have the top 100 genes remaining
top100_filtered_cortex_gene_abundances_matrix <- filtered_cortex_gene_abundances_matrix[rownames(filtered_cortex_gene_abundances_matrix) %in% top100, ]
top50_filtered_cortex <- filtered_cortex_gene_abundances_matrix[rownames(filtered_cortex_gene_abundances_matrix) %in% top50,]












## make heatmap
pheatmap(top100_filtered_cortex_gene_abundances_matrix)
pheatmap(top50_filtered_cortex)
## this looks awful because things aren't normalized so lets normalize
## make function to z-score
cal_z_score <- function(x) {
  (x- mean(x)) / sd(x)
}

## apply z-score function to top 100 gene abundances
top100_filtered_cortex_gene_abundances_matrix_z_score <- t(apply(top100_filtered_cortex_gene_abundances_matrix, 1, cal_z_score))
top50_filtered_cortex_gene_abundances_matrix_z_score <- t(apply(top50_filtered_cortex,1, cal_z_score))

## Make heatmap of top 100 DEGs APP vs WT (all ages) in z-score
pheatmap(top100_filtered_cortex_gene_abundances_matrix_z_score)
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score)

## Now perform hierarchical clustering to obtain gene cluseters
my_hclust_gene_100 <- hclust(dist(top100_filtered_cortex_gene_abundances_matrix_z_score), method = "complete")
my_hclust_gene_50 <- hclust(dist(top50_filtered_cortex_gene_abundances_matrix_z_score), method = "complete")

library(dendextend)
as.dendrogram(my_hclust_gene_100) %>%
    plot(horiz= T)
as.dendrogram(my_hclust_gene_50) %>%
    plot(horiz=T)
## Now we can form as many clusters as we want we cutree() function
## We'll do two here
## my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene_100), k=2 )
## my_gene_col
## my_gene_col_50 <- cutree(tree = as.dendrogram(my_hclust_gene_50), k=2)

## my_gene_col <- data.frame(cluster= ifelse(test = my_gene_col ==1, yes = "cluster 1", no= "cluster 2"))
## my_gene_col_50 <- data.frame(cluster= ifelse(test = my_gene_col_50 ==1, yes = "cluster 1", no= "cluster 2"))
## head(my_gene_col)
## head(my_gene_col_50)

my_sample_col <- data.frame(APP = sample_metadata$APP, Diet= sample_metadata$Diet, Age= sample_metadata$Age)
rownames(my_sample_col) <- sample_metadata$Sample
my_sample_col

pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col)
pheatmap(top100_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col)

## Trying to add different color scheme to the heatmap starting 10-18-2022
## usign https://biocorecrg.github.io/CRG_RIntroduction/pheatmap-function-from-the-pheatmap-package.html as a guide
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F)  #can do this as will get rid of weird gene clustering on left side of graph
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, cluster_cols = F) ##dont want to do this, you lose the clustering by genes

##change colors to rainbow in heatmpa
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = rainbow(50))  ##not a fan
## change colors to blue to red (middle as white)
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = colorRampPalette(c('navy', 'white', 'red'))(50))  ## this makes it looks red, white, and blue- also not a huge fan

## https://bookdown.org/rdpeng/exdata/plotting-and-color-in-r.html as a guide
library(RColorBrewer)
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = brewer.pal(11, "Spectral"))
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = brewer.pal(11, "PiYG"))
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = brewer.pal(11, "PRGn")) ## I like this one a lot
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = brewer.pal(11, "RdBu"))
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = brewer.pal(9, "Blues"))
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = brewer.pal(9, "Greens"))

top50_heatmap <- pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, cluster_rows = F, color = brewer.pal(11, "PRGn"))

col_dend <- top50_heatmap[[2]]
col_dend <- dendextend::rotate(col_dend, order= c('WT', 'APP'))
## trying to change annotation colors (ie APP, Age)
## https://bioinformatics.stackexchange.com/questions/18116/changing-the-order-of-colors-in-pheatmap as a guide to change annotation colors
annot_colors= list(APP=c(APP="red1", WT= "royalblue1"), Diet= c(Control= "olivedrab", Supplemented= "cyan"), Age= c('3'='darkgreen', '6'='darkorange', '9'= 'steelblue1', '12'= 'darkred'))
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors,  cluster_rows = F, color = brewer.pal(11, "PRGn"), border_color = "gray") ## this actually worked so need to play more with it



## This is combination I like and think I will use for poster
annot_colors_gradient= list(APP=c(APP="red1", WT= "royalblue1"), Diet= c(Control= "chartreuse2", Supplemented= "springgreen4"), Age= c('3'='gray99', '6'='gray75', '9'= 'gray50', '12'= 'gray25'))
my_sample_col$APP <- factor(my_sample_col$APP, levels = c('WT', 'APP'))
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors_gradient,  cluster_rows = F, color = brewer.pal(11, "PRGn"), border_color = "gray60")


## Trying to rotate the branch so that all APP on inside
my_hclust_sample_top50 <- hclust(dist(t(top50_filtered_cortex_gene_abundances_matrix_z_score)), method= "complete")
as.dendrogram(my_hclust_sample_top50) %>%
  plot()
## so this creates the dendrogram but now i just want the flip the first branch of this!! How do I do it???
annot_colors_gradient= list(APP=c(APP="red1", WT= "royalblue1"), Diet= c(Control= "gold", Supplemented= "goldenrod3"), Age= c('3'='gray99', '6'='gray75', '9'= 'gray50', '12'= 'gray25'))
top50_heatmap <- pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors_gradient,  cluster_rows = F, color = brewer.pal(11, "PRGn"), border_color = "gray60")
col_dend <- top50_heatmap[[2]]
col_dend
col_dend <- dendextend::rotate(col_dend, order= c('226', '189','217','287','216','202','190','155','204','27','76','67','31','183','46','211','125','72','79','139','25','24','236','33','124','57',
                                                  '58','268','203','242','185','265','103','213','174','186','149','130','129','106','161','157','152','164','134','141','178','122','74','18','249','20','251','8', '127',
                                                  '95','61','96','36','9', 
                                                  '275', '168', '84', '40', '194', '191', '132', '90', '41', '252', '277', 
                                                  '290', '160', '192', '199', '262', '89', '102', '272', '238', '285', '263', '153', '253', '283', '107', '147', '228',
                                                  '282', '207', '223', '117', '244', '221', '225', '259' ))
                                                  
                                                  
                          
pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors_gradient, cluster_cols = as.hclust(col_dend), cluster_rows = F, color = brewer.pal(11, "PRGn"), border_color = "gray50")
## You did it!!!!
top50_heatmap <- pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors=annot_colors_gradient, cluster_cols = as.hclust(col_dend), cluster_rows = F, color = brewer.pal(11, "PRGn"), border_color = "gray60")


## Need to adjust scale of z-scores to show that white is no change, purple is downregulated and greeen is upregulated
## Using https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r as a guide

palettelength2 <- 50
myColor_ex2 <- colorRampPalette(c("#762A83", "#9970AB", "#C2A5CF", '#E7D4E8', "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B"))(palettelength2)
myColor_ex2
## length(breaks) == length(palettelength) +1
## use floor and ceiling to deal with even/odd length pallettelenghts
myBreaks <- c(seq(min(top50_filtered_cortex_gene_abundances_matrix_z_score), 0, length.out=ceiling(palettelength2/2)+1),
              seq(max(top50_filtered_cortex_gene_abundances_matrix_z_score)/palettelength2, max(top50_filtered_cortex_gene_abundances_matrix_z_score), length.out=floor(palettelength2/2)))

cortex_top50_heatmap_overflow <- pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors = annot_colors_gradient, cluster_cols = as.hclust(col_dend), cluster_rows = F,
                                         color = myColor_ex2, border_color = "gray60", breaks = myBreaks,
                                         main = "Cortex Hierarchical Clustering of Top 50 DEGs",
                                         fontsize = 24,
                                         fontsize_row =  10,
                                         fontsize_col= 10)



## Now trying to get scale of green same as hippo graph
palettelength2 <- 50
myColor_ex2 <- colorRampPalette(c("#762A83", '#E7D4E8', "#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61"))(palettelength2)
myColor_ex2
## length(breaks) == length(palettelength) +1
## use floor and ceiling to deal with even/odd length pallettelenghts
myBreaks <- c(seq(min(top50_filtered_cortex_gene_abundances_matrix_z_score), 0, length.out=ceiling(palettelength2/2)+1),
              seq(max(top50_filtered_cortex_gene_abundances_matrix_z_score)/palettelength2, max(top50_filtered_cortex_gene_abundances_matrix_z_score), length.out=floor(palettelength2/2)))

cortex_top50_heatmap_overflow <- pheatmap(top50_filtered_cortex_gene_abundances_matrix_z_score, annotation_col = my_sample_col, annotation_colors = annot_colors_gradient, cluster_cols = as.hclust(col_dend), cluster_rows = F,
                                          color = myColor_ex2, border_color = "gray60", breaks = myBreaks,
                                          main = "Cortex Hierarchal Clustering of Top 50 DEGs",
                                          fontsize = 18,
                                          fontsize_row =  13,
                                          fontsize_col= 11)

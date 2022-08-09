## GSEA with Trimmed Cortex Data for APP vs WT (Diet, Sex, Age Adjusted)
## Using https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/



BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")

library(BiocManager)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Mm.eg.db)




## Read in DEG data
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
df <- read.csv("Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_7_11_2022.csv", header = T)

## Take the adjusted p-value
original_gene_list <- df$padj

## name the vector
names(original_gene_list) <- df$X

## Omit any NA values
gene_list <- na.omit(original_gene_list)

## Already sorted so don't really need to but will do it anyway; needs to be in decreasing order for 
gene_list <- sort(gene_list, decreasing = T)

## Check what key types are available for Org.Mm.eg.db
keytypes(org.Mm.eg.db)

## Run Gene Set Enrichment Analysis of Gene Ontology
gse <- gseGO( geneList = gene_list,
              ont = "ALL",
              keyType = "ENSEMBL",
              minGSSize = 5,
              maxGSSize = 10000,
              pvalueCutoff = 0.05,
              verbose = T,
              OrgDb = org.Mm.eg.db,
              pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

emapplot(gse, showCategory = 10)

ridgeplot(gse) + labs(x= "Enrichment Distribution")




gse_mf <- gseGO( geneList = gene_list,
              ont = "MF",
              keyType = "ENSEMBL",
              minGSSize = 5,
              maxGSSize = 10000,
              pvalueCutoff = 0.05,
              verbose = T,
              OrgDb = org.Mm.eg.db,
              pAdjustMethod = "none")
require(DOSE)
dotplot(gse_mf, showCategory=10, split=".sign") + facet_grid(.~.sign)

emapplot(gs_mf, showCategory = 10)

ridgeplot(gse_mf) + labs(x= "Enrichment Distribution")

cnetplot(gse_mf, categorySize= "pvalue", showCategory = 3, foldChange = gene_list)




gse_cc <- gseGO( geneList = gene_list,
                 ont = "CC",
                 keyType = "ENSEMBL",
                 minGSSize = 5,
                 maxGSSize = 10000,
                 pvalueCutoff = 0.05,
                 verbose = T,
                 OrgDb = org.Mm.eg.db,
                 pAdjustMethod = "none")
require(DOSE)
dotplot(gse_cc, showCategory=10, split=".sign") + facet_grid(.~.sign)

ridgeplot(gse_cc) + labs(x= "Enrichment Distribution")




gse_bp <- gseGO( geneList = gene_list,
                 ont = "BP",
                 keyType = "ENSEMBL",
                 minGSSize = 5,
                 maxGSSize = 10000,
                 pvalueCutoff = 0.05,
                 verbose = T,
                 OrgDb = org.Mm.eg.db,
                 pAdjustMethod = "none",
                 nproc=1)
require(DOSE)
dotplot(gse_bp, showCategory=10, split=".sign") + facet_grid(.~.sign)

ridgeplot(gse_bp) + labs(x= "Enrichment Distribution")

## Now trying as suggested by logfoldchange

## Read in DEG data
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
df <- read.csv("Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_7_11_2022.csv", header = T)

## Take the adjusted p-value
original_gene_list <- df$log2FoldChange

## name the vector
names(original_gene_list) <- df$X

## Omit any NA values
gene_list <- na.omit(original_gene_list)

## Already sorted so don't really need to but will do it anyway; needs to be in decreasing order for 
gene_list <- sort(gene_list, decreasing = T)

## Check what key types are available for Org.Mm.eg.db
keytypes(org.Mm.eg.db)

## Run Gene Set Enrichment Analysis of Gene Ontology
gse <- gseGO( geneList = gene_list,
              ont = "ALL",
              keyType = "ENSEMBL",
              minGSSize = 5,
              maxGSSize = 10000,
              pvalueCutoff = 0.05,
              verbose = T,
              OrgDb = org.Mm.eg.db,
              pAdjustMethod = "fdr")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

emapplot(gse, showCategory = 10)



## Convert GeneIDs to Gene Names
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
df <- read.csv("Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_7_11_2022.csv", header = T)
names(df)
Ensembl <- df$X
names(Ensembl) <- df$X
keytypes(org.Mm.eg.db)
symbols <- bitr(names(Ensembl), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)  ## 9.02% of Ensembl names dont have a symbol
unique_symbols <- symbols[!duplicated(symbols[c("ENSEMBL")]),]
duplicated(unique_symbols)
unique_symbols[duplicated(unique_symbols)]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$X %in% unique_symbols$ENSEMBL,]
## Create a new column in df2 with corresponding Gene symbols
df2$Symbol <- unique_symbols$SYMBOL
## Create a vector of the gene universe
kegg_gene_list <- df2$log2FoldChange
##Name vector with Symbols
names(kegg_gene_list) <- df2$Symbol
##imut any NA values
kegg_gene_list <- na.omit(kegg_gene_list)
## Sort in decreasing order
kegg_gene_list <- sort(kegg_gene_list, decreasing = T)

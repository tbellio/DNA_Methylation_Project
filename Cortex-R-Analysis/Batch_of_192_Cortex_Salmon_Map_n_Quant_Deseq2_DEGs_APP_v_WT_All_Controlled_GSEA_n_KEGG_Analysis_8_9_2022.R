## Looking for GSEA analysis with Cortex data

## Using https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html as a guide

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


## Tutorial
data(geneList, package = "DOSE")
gene <- names(geneList)[abs(geneList) >2]
ggo <-groupGO(gene = gene,
              OrgDb = org.Hs.eg.db,
              ont = "CC",
              level = 3,
              readable = T)
head(ggo)


## My Data
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
df <- read.csv("Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_7_11_2022.csv", header = T)
df_sig <- df[df$padj <= 0.05,]
df_sig <- na.omit(df_sig)
ids <- bitr(df_sig$X, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
df_sig <- ids$ENTREZID

df_sig_GO <- groupGO(gene = df_sig,
                     OrgDb = org.Mm.eg.db,
                     ont = "CC",
                     level = 3,
                     readable = T)
head(df_sig_GO)


## TUtorial GO Over-representation analysis
ego <- enrichGO(gene = gene,
                universe = names(geneList),
                OrgDb = org.Hs.eg.db,
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = T)
head(ego)


## My data GO Over-representation analysis
df_sig_EGO <- enrichGO(gene = df_sig,
                       OrgDb = org.Mm.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05,
                       readable = T)
head(df_sig_EGO)
head(df_sig_EGO, 5)
head(df_sig_EGO, 10)



## Tutorial GO Gene Set Enrichment Analaysis
ego3 <- gseGO(geneList = geneList,
              OrgDb = org.Hs.eg.db,
              ont = "CC",
              minGSSize = 100,
              maxGSSize = 500,
              pvalueCutoff = 0.05,
              verbose = T)
ridgeplot(ego3)
require(DOSE)
dotplot(ego3, showCategory=10, split=".sign") + facet_grid(.~.sign)

# My Data GO Gene Set Enrichment Analysis
## Need to make data so that first column is GeneID
## Second column is fold change
## and that it is sorted in decreasing order
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
df <- read.csv("Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_7_11_2022.csv", header = T)
df_sig <- df[df$padj <= 0.05,]
df_sig <- na.omit(df_sig)
ids <- bitr(df_sig$X, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
joint_df_sig <- merge(df_sig, ids, by.x= "X", by.y= "ENSEMBL")
joint_df_sig_genelist = joint_df_sig[,3]
names(joint_df_sig_genelist) = as.character(joint_df_sig[,"ENTREZID"])
joint_df_sig_genelist = sort(joint_df_sig_genelist, decreasing = T)
## ok all data is now in GeneList format called 'joint_df_sig_genelist
df_sig_GSEA <- gseGO(geneList = joint_df_sig_genelist,
                     OrgDb = org.Mm.eg.db,
                     ont = "CC",
                     minGSSize = 50,
                     maxGSSize = 1000,
                     pvalueCutoff = 0.05,
                     verbose = T)
ridgeplot(df_sig_GSEA)
require(DOSE)
dotplot(df_sig_GSEA, showCategory=10, split=".sign") + facet_grid(.~.sign)
goplot(df_sig_GSEA)

df_sig_GSEA2 <- gseGO(geneList = joint_df_sig_genelist,
                     OrgDb = org.Mm.eg.db,
                     ont = "CC",
                     minGSSize = 100,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = T)
ridgeplot(df_sig_GSEA2)
require(DOSE)
dotplot(df_sig_GSEA2, showCategory=10, split=".sign") + facet_grid(.~.sign)

goplot(df_sig_EGO)
barplot(df_sig_EGO)
require(DOSE)
dotplot(df_sig_GSEA2, showCategory=10, split=".sign") + facet_grid(.~.sign)




df_sig_GSEA3 <- gseGO(geneList = joint_df_sig_genelist,
                      OrgDb = org.Mm.eg.db,
                      ont = "MF",
                      minGSSize = 50,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      verbose = T)
ridgeplot(df_sig_GSEA3)
require(DOSE)
dotplot(df_sig_GSEA3, showCategory=10, split=".sign") + facet_grid(.~.sign)
goplot(df_sig_GSEA3)



df_sig_GSEA4 <- gseGO(geneList = joint_df_sig_genelist,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      minGSSize = 50,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      verbose = T)
ridgeplot(df_sig_GSEA4)
require(DOSE)
dotplot(df_sig_GSEA4, showCategory=10, split=".sign") + facet_grid(.~.sign)
goplot(df_sig_GSEA4)


search_kegg_organism('mus', by="scientific_name")
mouse <- search_kegg_organism('Mus musculus', by='scientific_name')
dim(mouse)


## KEGG Analysis
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
df <- read.csv("Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_7_11_2022.csv", header = T)
df_sig <- df[df$padj <= 0.05,]
df_sig <- na.omit(df_sig)
ids <- bitr(df_sig$X, fromType = "ENSEMBL", toType = "UNIPROT", OrgDb = "org.Mm.eg.db")
df_sig <- ids$UNIPROT
df_sig_enrichKEGG <- enrichKEGG(gene=df_sig,
                                organism = 'mmu',
                                keyType = 'uniprot',
                                pvalueCutoff = 0.05)
head(df_sig_enrichKEGG)  

## KEGG Pathway Gene Set Enrichment Analysis
## Need to change ensembl ids firt though
## Using http://bioconductor.riken.jp/packages/3.4/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#how-to-build-a-biomart-query as a guide
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Cortex_Salmon_Mapped_and_Quant_R_Analysis")
df <- read.csv("Batch_of_192_Trimmed_Cortex_DESeq2_APP_v_WT_DGE_Analysis_all_controlled_padj_sorted_7_11_2022.csv", header = T)
df_sig <- df[df$padj <= 0.05,]
df_sig <- na.omit(df_sig)
keytypes(org.Mm.eg.db)
library(biomaRt)
listMarts()
useMart("ENSEMBL_MART_MOUSE")
mouse = useMart("ENSEMBL_MART_MOUSE")
listDatasets(mouse)
mouse = useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
listFilters(mouse)
ncbi_id <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),
                 filters = 'ensembl_gene_id',
                 values = df_sig$X,
                 mart = mouse)

merged_ensembl_geneid <- merge(df_sig, ncbi_id, by.x="X", by.y= 'ensembl_gene_id')
merged_ensembl_geneid_geneList <- merged_ensembl_geneid[,3]
names(merged_ensembl_geneid_geneList) <- as.character(merged_ensembl_geneid[,"entrezgene_id"])
merged_ensembl_geneid_geneList = sort(merged_ensembl_geneid_geneList, decreasing=T)
## ok now we have a 'geneList' to use for the KEGG pathway gene set enrichment analysis


df_sig_gseKEGG <- gseKEGG(geneList = merged_ensembl_geneid_geneList,
                          organism = 'mmu',
                          minGSSize = 50,
                          pvalueCutoff = 0.05,
                          verbose = T)
head(df_sig_gseKEGG)
dotplot(df_sig_gseKEGG)

## KEGG Module gene set enrichment analysis
df_sig_gseMKEGG <- gseMKEGG(geneList = merged_ensembl_geneid_geneList,
                            organism = 'mmu',
                            pvalueCutoff = 1)
head(df_sig_gseMKEGG)
dotplot(df_sig_gseMKEGG)

head(df_sig_enrichKEGG)
browseKEGG(df_sig_enrichKEGG, 'mmu04062')
browseKEGG(df_sig_enrichKEGG, 'mmu05171')


library(pathview)
mmu04062 <- pathview(gene.data = merged_ensembl_geneid_geneList,
                     pathway.id = 'mmu04062',
                     species = 'mmu',
                     limit = list(gene=max(abs(merged_ensembl_geneid_geneList)), cpd=1))
view(mmu04062)
mmu04062
##that didnt work, try another pathway
mmu05171 <- pathview(gene.data = merged_ensembl_geneid_geneList,
                     pathway.id = 'mmu05171',
                     species = 'mmu',
                     limit = list(gene=max(abs(merged_ensembl_geneid_geneList)), cpd=1))
##this isnt workign either; saying it cant find the pathway even though when i use KEGG browser it can find it
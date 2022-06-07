##trying to use tximeta to do DGE with edgeR with correct gene id-transcript id matching
## Using RNA-Seq data from Batch of 192 Hippocampus Salmon Mapped and Quantified


library(BiocManager)
library(tximport)
library(tximportData)
library(readxl)
library(tximeta)
library(apeglm)
library(edgeR)
library(GO.db)
library(edgeR)

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


##sometimes have to run this bit twice for some reason to work properly
dir <- setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_aligned_Against_Mouse_cdna_Output_5_2022")
files <- file.path(dir, samples$names,  "quant.sf" )
file.exists(files)


## Use coldata to add in the metadata to the files data frame
coldata <- data.frame(files, names= samples$names, APP= samples$APP, Age= samples$Age, Diet= samples$Diet, Sex= samples$Sex, stringsAsFactors = F)
coldata

## create se which is large ranged summarized experiment that includes all transcripts read and aligned to reference 
se <- tximeta(coldata)  ##creates se which is a very large summarized experiment dataset....IT FUCKING WORKED...well maybe
colData(se)
assayNames(se)
rowRanges(se)
seqinfo(se)

##Retreive the database used for se
edb <- retrieveDb(se)
##check class of the database and it should be an ensembl package
class(edb)

## Take se and summarize it to gene level instead of transcript level
gse <- summarizeToGene(se)
rowRanges(gse)

##now using edgeR to look for DGE
##have to create matrices of the group information and the counts information from the tximeta output 'gse'
library(edgeR)
group <- gse@colData@listData[["APP"]]
counts <- gse@assays@data@listData[["counts"]]
##create a DGElist to be used for getting DGE in the samples
y <- DGEList(counts = counts, group = group)
y$samples
y$counts
##filter out lowly expressed genes
keep <- filterByExpr(y)
y_filtered <- y[keep, , keep.lib.sizes=F]  ##essentially goes from ~35000 down to ~15000 genes

##Need to Normalize the data for library size effects (ie need to make sure some genes arent being expressed more in a particular sample just because its library is bigger than another sample)
y_filtered_and_normalized <- calcNormFactors(y_filtered)
y_filtered_and_normalized$samples

##Estimate dispersions
##estimate common and tagwise dispersion in one run
y_filtered_and_normalized_dispersion <- estimateDisp(y_filtered_and_normalized)

##get DEGs using log counts per million
exact_test <- exactTest(y_filtered_and_normalized_dispersion)
sorted_exact_test <- exact_test[order(exact_test$table$PValue),]

summary(sorted_exact_test$table)

setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
write.csv(as.data.frame(sorted_exact_test$table),
          file = "Batch_of_192_Hippocampus_Salmon_Map_n_Quant_edgeR_Results_using_qCML_5_20_2022.csv")  

### Now trying using generalized linear models (GLMs)
group_matrix <- as.matrix(group)
group_matrix <- model.matrix(~group)
y_CR_disperion <- estimateDisp(y_filtered, group_matrix)

fit <- glmQLFit(y_CR_disperion, group_matrix)
APP_vs_WT <- glmQLFTest(fit, coef = 2)
topTags(APP_vs_WT)

sorted_APP_vs_WT_glm <- APP_vs_WT[order(APP_vs_WT$table$PValue),]

write.csv(as.data.frame(sorted_APP_vs_WT_glm),
          file = "Batch_of_192_Hippocampus_Salmon_Map_n_Quant_edgeR_results_using_glm_from_5_20_2022.csv")

##Gene ontology and pathway analysis
library(GO.db)
library(org.Mm.eg.db)
APP_vs_WT <- glmQLFTest(fit, coef = 2)

APP_vs_WT$table$genes <- rownames(APP_vs_WT$table)

z <- select(org.Mm.eg.db, keys(org.Mm.eg.db), "SYMBOL")

APP_vs_WT$table$keys <- keys(org.Mm.eg.db)  ##this doesnt work
##this creates a data.frame by using the row names from APP_vs_WT$table and matching the corresponding ENTREZID to the row.names (which are ensembl genes) from the org.Mm.eg.db (which is an online database of a bunch of annotations)
annot <- select(org.Mm.eg.db, row.names(APP_vs_WT$table), "SYMBOL", keytype = "ENSEMBL")  ##Holy shit it finally fucking worked!!!!!
##but there are more genes when creatign this annot than there are as rows in APP_vs_WT which doesnt make sense; means theres probably duplications
length(unique(annot$ENSEMBL))  ##15011 ensembl genes
length(unique(annot$SYMBOL))  ##14640 symbol which means more ensembl genes (ENSUMG000) than gene symbols
length(unique(annot$ENSEMBL)) == nrow(annot)  ##seeing if the number of rows in annot is the same as unique ensembl genes...it is not
## so we want to turn annot so that there is no duplicated gene symbols; it will keep the first one that shows up here
annot <- annot[!duplicated(annot$SYMBOL),]
library(dplyr)
nonduplicated_annot <- distinct(annot, ENSEMBL, .keep_all = T)
row.names(nonduplicated_annot) <- nonduplicated_annot$ENSEMBL
write.csv(as.data.frame(nonduplicated_annot),
          file = "Ensemble_matched_Entrezid_annotation.csv")
nonduplicated_annot$ENSEMBL <- NULL ##deletes the ensembl column since now column names are labeled with ensembl id
##run goana; must have column names be ENTREZID for geneid=
gene_ontology <- goana(APP_vs_WT, geneid = nonduplicated_annot$ENTREZID,  species= "Mm")
gene_ontology <- gene_ontology[order(gene_ontology$P.Up),]
top_gene_ontology <- topGO(gene_ontology, sort = "up")

write.csv(as.data.frame(gene_ontology),
          file = "Gene_Ontology_Analysis_of_App_v_WT_using_goana_v1_8_18_2021.csv")


##Running KEGG Pathway analysis
#can just use generalized linear model output (APP_vs_WT) from line of code 'APP_vs_WT <- glmQLFTest(fit, coef = 2)'
kegg_pathway <- kegga(APP_vs_WT, geneid = nonduplicated_annot$ENTREZID, species= "Mm")
kegg_pathway <- kegg_pathway[order(kegg_pathway$P.Up),]
topKEGG(kegg_pathway)

write.csv(as.data.frame(kegg_pathway),
          file = "Pathway_Analysis_of_APP_v_WT_using_kegga_v1_8_18_2021.csv")

##need to somehow get all of my ensembl gene ids; just kidding i didnt need to do this, just needed to figure out how to run the line of code that creates 'annot' properly
##but none the less, this is still useful as now i have created a nice index that aligns all the Refseq, Ensembl, Uniprot, NCBI, GO annotations together
library(org.Mm.eg.db)
columns(org.Mm.eg.db)
uniKeys <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
cols <- c("ENTREZID", "GENENAME", "GENETYPE", "GO", "ONTOLOGY", "REFSEQ", "SYMBOL", "UNIPROT")
ensm_ncbi_go_refseq_annot <- select(org.Mm.eg.db, keys = uniKeys, columns = cols, keytype = "ENSEMBL")

write.csv(as.data.frame(ensm_ncbi_go_refseq_annot),
          file = "Annotation_Combination_of_Ensembl_Refseq_GO_for_mouse.csv")

##plotMDS
plotMDS(APP_vs_WT$table)
plotMDS(y_filtered)


##counts per million with filtered and normalized data
logcpm <- cpm(y_filtered_and_normalized, log = T)


##heatmap using pheatmap
install.packages('pheatmap')
library(pheatmap)
pheatmap(filtered_dds_deseq) ##doesnt work; need a matrix
filtered_dds_deseq_matrix <- as.matrix(filtered_dds_deseq@assays@data@listData[["abundance"]])
subset_filtered_dds_deseq_matrix <- filtered_dds_deseq_matrix[rowSums(filtered_dds_deseq_matrix)>25000,]
pheatmap(subset_filtered_dds_deseq_matrix)


#heatmap using 

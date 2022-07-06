

library(BiocManager)
library(tximport)
library(readxl)
library(tximeta)
library(apeglm)
library(GO.db)
library(org.Mm.eg.db)
library(ggplot2)
library('DESeq2')
library(pheatmap)
library(dendextend)
library(readr)
library(AnnotationDbi)

Hippo_nlcRNAs <- read_excel("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/B192_Hippocampus_nlcRNAs_only_6_29_2022.xlsx")

geneids <- as.data.frame(Hippo_nlcRNAs$Geneid)
geneids <- as.matrix(geneids)

justgeneids <- gsub("\\..*", "", geneids)

justgeneids <- as.data.frame(justgeneids)

colnames(justgeneids) <- "GeneIDNew"

write_delim(justgeneids, "Mouse_nlcRNA_list.txt", delim = " ")

write.table(justgeneids, "Mouse_nlcRNAs_list.txt", append=F, sep= " ", row.names = F, col.names = F)

colnames(Hippo_nlcRNAs)

Hippo_nlcRNAs$GeneIDNew <- justgeneids
Hippo_nlcRNAs$GeneIDNew <- as.character(Hippo_nlcRNAs$GeneIDNew)

## Need to have nlcRNA ensembl data as vector for next line in code
nlcRNAgeneids <- as.vector(Hippo_nlcRNAs$GeneIDNew)

nlcRNA_symbols <- AnnotationDbi::select(org.Mm.eg.db, nlcRNAgeneids, "SYMBOL", keytype = "ENSEMBL")
## ok so this isnt really working lets try something else


## Will create two data frames then merge with the overlapping information which will be the Ensembl gene id for this
Gene_name_and_ensembl <- read_excel("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Mouse_nlcRNAs_gene_name_ensemblid_list.xls")
Gene_name_and_ensembl <- as.data.frame(Gene_name_and_ensembl)
merged <- merge(Hippo_nlcRNAs, Gene_name_and_ensembl, by.x= "GeneIDNew", by.y= "Gene stable ID", all.x= T)

charactervectorgeneidnew <- Hippo_nlcRNAs[['GeneIDNew']]
charactervectorgeneidnew2 <- as.vector(Hippo_nlcRNAs['GeneIDNew'])
class(charactervectorgeneidnew2)
Hippo_nlcRNAs$GeneIDNew2 <- charactervectorgeneidnew2
charactervectorgeneidnew3 <- as.character(Hippo_nlcRNAs['GeneIDNew'])
class(charactervectorgeneidnew3)
Hippo_nlcRNAs$GeneIDNew3 <- charactervectorgeneidnew3 ## dont do this it fucks it all up
class(Hippo_nlcRNAs$GeneIDNew2)
class(Hippo_nlcRNAs$Geneid)
Hippo_nlcRNAs$GeneIDNew <- as.character(Hippo_nlcRNAs$GeneIDNew)
class(Hippo_nlcRNAs$GeneIDNew)

merged <- subset(merged, select= -c(GeneIDNew2, GeneIDNew3))


test <- test[-c(2:6)]
test_matrix <- as.matrix(test[,-1])
rownames(test_matrix) <- test$Geneid
metadata <- read_excel("C:/Users/tbell/Downloads/B192_test_v2_meta.xlsx")
metadata <- as.matrix(metadata)

dds <- DESeqDataSetFromMatrix(countData = test_matrix,
                              colData = metadata,
                              design = ~dex)
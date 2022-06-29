



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

test <- read_excel("C:/Users/tbell/Downloads/B192_test_v2.xlsx")
test <- test[-c(2:6)]
test_matrix <- as.matrix(test[,-1])
rownames(test_matrix) <- test$Geneid
metadata <- read_excel("C:/Users/tbell/Downloads/B192_test_v2_meta.xlsx")
metadata <- as.matrix(metadata)

dds <- DESeqDataSetFromMatrix(countData = test_matrix,
                              colData = metadata,
                              design = ~dex)

## RNA-seq Analysis easy as 1-2-3 with limma, glimma, and edgeR tutorial (5-20-2022)
## Tutorial found https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/
## Using RNA-Seq data from Batch of 192 Hippocampus Salmon Mapped and Quantified

library(BiocManager)
library(tximport)
library(tximportData)
library(readxl)
library(tximeta)
library(apeglm)
library(edgeR)
library(GO.db)
library(DESeq2)
library(limma)
library(Glimma)


##tutorial uses readDGE but I can just use tximeta to import all my data from salmon

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

##make DGEList from gse counts called test
## remember that counts is not normalized at all, will do normalization later down in this code
test <- DGEList(counts = gse@assays@data@listData[["counts"]])
class(test)
dim(test)

##need to add sample information to the samples dataframe of test
## Add APP genotype to test DGElist
APP <- as.factor(samples$APP)
test$samples$APP <- APP
## Add Age to test DGElist
Age <- as.factor(samples$Age)
test$samples$Age <- Age
## Add Diet to test DGElist
Diet <- as.factor(samples$Diet)
test$samples$Diet <- Diet
## Ad Sex to test DGElist
Sex <- as.factor(samples$Sex)
test$samples$Sex <- Sex


## load in the ensembl mouse reference
library(ensembldb)
library(org.Mm.eg.db)
library(AnnotationFilter)
edb <- org.Mm.eg.db
edb

## create a vector called ensembl_gene_names that uses the rownames of tests 
ensembl_gene_names <- rownames(test)

##check the key types of the ensembl mouse reference- key types are the different kinds of variables within the org.Mm.eg.db
## key types shoudl be ensembl, go, entrezid, symbol, ontology, etc
keytypes(org.Mm.eg.db)
columns(org.Mm.eg.db)

##create genes data frame that uses the ensembl mouse reference genome, and uses our ensembl gene names vector (which is just the genes mapped), and takes the gene symbol, chromosome location, and entrezid (NCBI gene id) and puts it into a data frame called genes
genes <- select(org.Mm.eg.db, keys = ensembl_gene_names, columns = c("SYMBOL", "CHRLOC", "ENTREZID"), keytype = "ENSEMBL")
##this created an error message when creating genes 
genes
## however we get more mapped places than we have (40522 genes to 35682 in test)
## so we need to do something with the duplicate maps
## here we select one of the chromosome mappigns to represent a gene with duplicate annotation
## essentially what the code is saying is to keep all the genes that are present and the ones that are located in more than one place to take just the first occurence of that gene
genes <- genes[!duplicated(genes$ENSEMBL),]
##genes should now be down to 35682 from the 40522-which means there were roughly 5000 genes that were duplicated and we took just the first occurence for these

##add our genes data frame to our test DGEList
test$genes <- genes
##check to make sure that the genes are ordered properly in the test DGElist--this should be true
genes %in% test$genes

##now test has 3 elements within it (counts, samples, and genes)

## Now we use edgeR to get counts per million and log2 counts per million
library(edgeR)
cpm_test <- cpm(test)
log_cpm_test <- cpm(test, log = T)
## note that cpm and log_cpm do not take into effect the size of gene length so it is not good to compare "expression level" between two different genes in these but is okay to look for differential expression of the same genes within the samples since they are comparing the same thing within a gene 
## a log cpm of 1 means having 20 counts in a sample with 20 million reads, or 50 counts in a sample with 50 million reads etc

## to get rpkm can just use rpkm() function; however need to have gene size given; can create gene lengths by using Iranges from gse
gene_lengths <- gse@rowRanges@ranges@width
test$genes$Length <- gene_lengths
rpkm_test <- rpkm(test)


## we will use log-cpm here
## log-cpm uses the average library size to make sure that values of same cpm are the same log-cpm
L <- mean(test$samples$lib.size) *1e-6
M <- median(test$samples$lib.size) *1e-6
c(L,M)
## so mean library size of all the samples is 28.51573 million reads and median is 28.50237 million reads 

summary(log_cpm_test)

## look for lowly expressed genes, essentially asking how many genes have no counts in all 96 samples
table(rowSums(test$counts==0)==96)
## says that 10344 genes are not expressed, while 25338 genes have atleast one read 

## filter genes by filterByExpr() function
keep.exprs <- filterByExpr(test, group = APP)  ##this creates a logical where if a gene has 10 or more read counts in a minimum number of samples than it is true, if not it is false
test_filtered <- test[keep.exprs,,keep.lib.sizes=F] ##cuts DGElist down to 15011 genes
dim(test_filtered)
log_cpm_test_filtered <- cpm(test_filtered, log=T)

## graph before and after filtering
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
samplenames <- colnames(test)
nsamples <- ncol(test)
col <- brewer.pal(nsamples, "Paired")

par(mfrow=c(1,2))
plot(density(log_cpm_test[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(log_cpm_test[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

plot(density(log_cpm_test_filtered[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(log_cpm_test_filtered[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")



## Save the filtered and unfiltered log cpm graphs in a pdf
setwd("C:/Users/tbell/Documents/Boston University/DNA_Methylation/RNA_Seq/Batch_of_192/Hippocampus_Salmon_Mapped_and_Quant_R_Analysis")
pdf("Unfiltered_v_Filtered_LogCPM.pdf")
par(mfrow=c(1,2))
plot(density(log_cpm_test[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(log_cpm_test[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

plot(density(log_cpm_test_filtered[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(log_cpm_test_filtered[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.off()

## now going to normalize using calcNormFactors() function which takes into account the library size 
test_filtered_TMM_norm <- calcNormFactors(test_filtered, method="TMM")
test_filtered_TMM_norm$samples$norm.factors
log_cpm_test_filtered_TMM_norm <- cpm(test_filtered_TMM_norm, log=T)


## Graph differences between normalized and unnormalized logcpm 
##graph unnormalized filtered test (logcpm)
par(mfrow=c(1,1.5))
boxplot(log_cpm_test_filtered, las=2, col=col, main="")
title(main="A. Example: Unnormalised data", ylab="Log-cpm")
##graph normalized filtered test (logcpm)
boxplot(log_cpm_test_filtered_TMM_norm, las=2, col=col, main="")
title(main="B. Example: Normalised data", ylab="Log-cpm")

##Save the graph of TMM normalized vs unnormalized log cpm data
pdf("TMM_Normalization_of_Filtered_log_CPM.pdf")
par(mfrow=c(1,1.5))
boxplot(log_cpm_test_filtered, las=2, col=col, main="")
title(main="A. Example: Unnormalised data", ylab="Log-cpm")
##graph normalized filtered test (logcpm)
boxplot(log_cpm_test_filtered_TMM_norm, las=2, col=col, main="")
title(main="B. Example: Normalised data", ylab="Log-cpm")
dev.off()

## plot MDS to see which samples group together without any outside influence
plotMDS(log_cpm_test_filtered_TMM_norm)

## plot using plotMDS based on APP genotype
par(mfrow=c(1,2))
col.APP <- APP
levels(col.APP) <-  brewer.pal(nlevels(col.APP), "Set1")
col.APP <- as.character(col.APP)
col.Age <- Age
levels(col.Age) <-  brewer.pal(nlevels(col.Age), "Set2")
col.Age <- as.character(col.Age)
plotMDS(log_cpm_test_filtered_TMM_norm, labels=APP, col=col.APP)
title(main="A. Sample groups")
plotMDS(log_cpm_test_filtered_TMM_norm, labels = Age, col=col.Age, dim= c(3,4))
title(main = "B. Age")

## save these two graphs in pdf format
pdf("Unsupervised_Clustering_of_Filtered_and_TMM_Normalized_Data_APP_and_Age.pdf")
par(mfrow=c(1,1.5))
col.APP <- APP
levels(col.APP) <-  brewer.pal(nlevels(col.APP), "Set1")
col.APP <- as.character(col.APP)
col.Age <- Age
levels(col.Age) <-  brewer.pal(nlevels(col.Age), "Set2")
col.Age <- as.character(col.Age)
plotMDS(log_cpm_test_filtered_TMM_norm, labels=APP, col=col.APP)
title(main="A. Sample groups")
plotMDS(log_cpm_test_filtered_TMM_norm, labels = Age, col=col.Age, dim= c(3,4))
title(main = "B. Age")
dev.off()

## plot using plotMDS based on APP genotype, age, and diet
par(mfrow=c(1,3))
col.APP <- APP
levels(col.APP) <-  brewer.pal(nlevels(col.APP), "Set1")
col.APP <- as.character(col.APP)
col.Age <- Age
levels(col.Age) <-  brewer.pal(nlevels(col.Age), "Set2")
col.Age <- as.character(col.Age)
col.Diet <- Diet
levels(col.Diet) <- brewer.pal(nlevels(col.Diet), "Set3")
col.Diet <- as.character(col.Diet)
plotMDS(log_cpm_test_filtered_TMM_norm, labels=APP, col=col.APP)
title(main="A. Sample groups")
plotMDS(log_cpm_test_filtered_TMM_norm, labels = Age, col=col.Age)
title(main = "B. Age")
plotMDS(log_cpm_test_filtered_TMM_norm, labels = Diet, col= col.Diet)
title(main = "C. Diet")

## Print dimension 1 and 2 with APP, Age, Diet
pdf("Unsupervised_Clustering_of_Filtered_and_TMM_Normalized_Data_APP_Diet_and_Age.pdf")
par(mfrow=c(1,1))
col.APP <- APP
levels(col.APP) <-  brewer.pal(nlevels(col.APP), "Set1")
col.APP <- as.character(col.APP)
col.Age <- Age
levels(col.Age) <-  brewer.pal(nlevels(col.Age), "Set2")
col.Age <- as.character(col.Age)
col.Diet <- Diet
levels(col.Diet) <- brewer.pal(nlevels(col.Diet), "Set3")
col.Diet <- as.character(col.Diet)
plotMDS(log_cpm_test_filtered_TMM_norm, labels=APP, col=col.APP)
title(main="A. Sample groups")
plotMDS(log_cpm_test_filtered_TMM_norm, labels = Age, col=col.Age)
title(main = "B. Age")
plotMDS(log_cpm_test_filtered_TMM_norm, labels = Diet, col= col.Diet)
title(main = "C. Diet")
dev.off()

library(Glimma)
glMDSPlot(log_cpm_test_filtered_TMM_norm, labels = paste(APP, Age, Diet, sep = "_"), groups = test_filtered_TMM_norm$samples[,c(4,5,6)],
          launch = TRUE, path = getwd(), folder = "GLimma_Plots_v1_5_20_2022", html = "MDS-Plot")


## Differential Gene Expression Analysis
## First need to create a design matrix
design <- model.matrix(~0+APP)
colnames(design) <- c("APP", "WT")
design

## Make Contrast
contr.matrix <- makeContrasts(
  APPvsWT = APP- WT,
  levels = colnames(design)
)

##voom filtering; When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation factors from x itself.
## It has been shown that for RNA-seq count data, the variance is not independent of the mean 13 - this is true of raw counts or when transformed to log-CPM values. Methods that model counts using a Negative Binomial distribution assume a quadratic mean-variance relationship. In limma, linear modelling is carried out on the log-CPM values which are assumed to be normally distributed and the mean-variance relationship is accommodated using precision weights calculated by the voom function
voom_test_filtered <- voom(test_filtered, design, plot = T)
voom_test_filtered
voom_lmfit_test_filtered <- lmFit(voom_test_filtered, design)
contrasts_voom_fit_test_filtered <- contrasts.fit(voom_lmfit_test_filtered, contrasts =  contr.matrix)
Bayes_contrasts_voom_fit_test_filtered <- eBayes(contrasts_voom_fit_test_filtered)
Bayes_voom_fit_test_filtered <- eBayes(voom_fit_test_filtered)
plotSA(Bayes_voom_fit_test_filtered)
plotSA(Bayes_contrasts_voom_fit_test_filtered)

#3plot both on same screen
par(mfrow=c(1,2))
voom_test_filtered <- voom(test_filtered, design, plot = T)
plotSA(Bayes_voom_fit_test_filtered, main= "Bayes-Voom:Mean-variance trend")

## save pdf of both voom and bayes-voom
pdf("Voom_Bayes_Log2_stdev_v1_5_20_2022.pdf")
par(mfrow=c(1,2))
voom_test_filtered <- voom(test_filtered, design, plot = T)
plotSA(Bayes_voom_fit_test_filtered, main= "Bayes-Voom:Mean-variance trend")
dev.off()


## THIS IS WHERE I DONT KNOW WHATS GOING ON



## find DEGs
summary(decideTests(voom_fit_test_filtered))  ##well fucked this up LOL im tired
summary(decideTests(Bayes_voom_fit_test_filtered))
summary(decideTests(Bayes_contrasts_voom_fit_test_filtered, adjustment.method= "fdr", p.value=0.000003))
summary(decideTests(Bayes_contrasts_voom_fit_test_filtered, adjustment.method= "fdr", p.value=0.9))
summary(decideTests(Bayes_contrasts_voom_fit_test_filtered, adjustment.method= "fdr", method = "nestedF"))
et <- exactTest(test_filtered_TMM_norm, pair = c("WT", "APP"))

design2 <- model.matrix(~0+APP+Diet+Age, data = test_filtered_TMM_norm$samples)
colnames(design2) <- c('APP', 'WT', 'Supplemented', '3', '6', '9')
design2
Group <- c('6SuppApp', '6ControlAPP', '6SuppWT', '3ControlAPP', '3SuppAPP',  )


design3 <- model.matrix(~APP + APP:Age, data = test_filtered_TMM_norm)
dispersions_test_filtered_TMM_norm <- estimateDisp(test_filtered_TMM_norm, design3)
glmwlfit_test_filtered_TMM_norm <- glmQLFit(dispersions_test_filtered_TMM_norm, design3)
colnames(glmwlfit_test_filtered_TMM_norm)
qlfAPPvsWT <- glmQLFTest(glmwlfit_test_filtered_TMM_norm, coef = 2)
topTags(qlfAPPvsWT)

qlfAPPvsWT9mo <- glmQLFTest(glmwlfit_test_filtered_TMM_norm, coef = 8)
topTags(qlfAPPvsWT9mo)


##design matrix with APP, APPand age, and Age
design4 <- model.matrix(~APP+ APP:Age + Age, data = test_filtered_TMM_norm)
dispersions_test_filtered_TMM_norm <- estimateDisp(test_filtered_TMM_norm, design4)
glmwlfit_test_filtered_TMM_norm <- glmQLFit(dispersions_test_filtered_TMM_norm, design4)
colnames(glmwlfit_test_filtered_TMM_norm)
qlfAPPvsWT <- glmQLFTest(glmwlfit_test_filtered_TMM_norm, coef = 2)
topTags(qlfAPPvsWT)



design5 <- model.matrix(~0+Age, test_filtered_TMM_norm$samples)
dispersions_test_filtered_TMM_norm <- estimateDisp(test_filtered_TMM_norm, design5)
glmwlfit_test_filtered_TMM_norm <- glmQLFit(dispersions_test_filtered_TMM_norm, design5)
colnames(glmwlfit_test_filtered_TMM_norm)
qlfAPPvsWT <- glmQLFTest(glmwlfit_test_filtered_TMM_norm, coef = 2)
topTags(qlfAPPvsWT)

design6 <- model.matrix(~0+APP, test_filtered_TMM_norm$samples)
dispersions_test_filtered_TMM_norm <- estimateDisp(test_filtered_TMM_norm, design6)
glmwlfit_test_filtered_TMM_norm <- glmQLFit(dispersions_test_filtered_TMM_norm, design6)
colnames(glmwlfit_test_filtered_TMM_norm)
qlfAPPvsWT <- glmQLFTest(glmwlfit_test_filtered_TMM_norm, coef = 2)
topTags(qlfAPPvsWT)

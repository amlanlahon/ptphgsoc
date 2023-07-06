setwd("C:/Users/Chocobunbhai/Documents")
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("hugene10stv1cdf")
library(makecdfenv)
library(affy)
library(hugene10stv1cdf)
library(limma)
library(hursta2a520709cdf)
library(ComplexHeatmap)
library(circlize)
library(preprocessCore)

#Create CDF package in temp directory 

pkgpath <- tempdir()

make.cdf.package("GPL10379_HuRSTA-2a520709_custom_MMPM.cdf", 
                 cdf.path <- "C:/Users/Chocobunbhai/Documents", 
                 compress <- FALSE, 
                 species = "Homo_sapiens",
                 packagename = "hursta2a520709cdf",
                 package.path = pkgpath)

# Note the temporary location used for the source package

dir(pkgpath)

# USE TERMINAL to run :- sudo R CMD INSTALL C:\Users\CHOCOB~1\AppData\Local\Temp\RtmpY9MkiL/hursta2a520709cdf

#RMA
ovar.data <- ReadAffy()
data.rma <- affy::rma(ovar.data)
write.exprs(data.rma, "expression_values_logFC1.csv")
str(ovar.data)
pheno <- pData(data.rma)
#pheno$batch = c(rep(1,62), rep(2,18))
edata <- exprs(data.rma)
#mod <- model.matrix(~as.factor(treatment), data <- pheno)
#mod0 <- model.matrix(~1,data <- pheno)
batch <- pheno$batch
#batch
#<-modcombat <- model.matrix(~1, data <- pheno)
#modcombat # equal to mod0
#combat_edata <- ComBat(dat <- edata, batch<-batch, mod <- modcombat, par.prior <- TRUE, prior.plots <- FALSE)
data.rma2 <- exprs(data.rma)

#Limma
#http://genomicsclass.github.io/book/pages/using_limma.html

cov <- read.table("covdesc.txt")
fac2 <- cov$Treatment
fac2
fit <- lmFit(data.rma2, design<-model.matrix(~ fac2))
colnames(coef(fit))
fit <- eBayes(fit)
#tt <- topTable(fit, coef=2)
res_table <- topTable(fit, coef=2, number=Inf, sort.by="none")
#write.table(res_table, "limma_results.txt")
#dim(res_table)
res_new <- res_table[order(res_table$adj.P.Val),]
write.csv(res_new, "limma_results_logFC1.csv")

#Annotation
#Download the SOFT Annotation file from https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10379/annot/ and process it into a CSV file 

array_annot = read.csv("annotations.csv")
res_new$annotation = array_annot[match(rownames(res_new), array_annot$ID),10]
res_new$description = array_annot[match(rownames(res_new), array_annot$ID),2]
write.csv(res_new, "result__with annotation_logFC1.csv")

#Threshold (p value > 0.01 & |logFC| > 1)

res_new$threshold = as.logical[(res_new$adj.P.Val < 0.01) & (abs(res_new$logFC)>1)]
degs_2fc =  res_new[which(res_new$threshold),]
write.csv(degs_2fc, "two_fold_change_degs_logFC1.csv")
write.csv(data.rma2[match(rownames(degs_2fc), rownames(data.rma2)),], "expression_values_of_degs_logFC1.csv")
exp_data_logFC1 = read.csv("expression_values_of_degs_logFC1.csv", row.names = 1)

#Normalization & Scaling

exp_data.norm_logFC1 = normalize.quantiles(as.matrix(exp_data_logFC1))
max(exp_data.norm_logFC1)
min(exp_data.norm_logFC1)
exp_data.norm.scaled_logFC1 <- scale(exp_data.norm_logFC1, center <- T, scale <- T)
max(exp_data.norm.scaled_logFC1)
min(exp_data.norm.scaled_logFC1)
colnames(exp_data.norm.scaled_logFC1) <- colnames(exp_data_logFC1)
rownames(exp_data.norm.scaled_logFC1) <- rownames(exp_data_logFC1)
ids <- array_annot[match(rownames(exp_data.norm.scaled_logFC1), array_annot$ID),]

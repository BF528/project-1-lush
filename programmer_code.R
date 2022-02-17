#part 1
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))

#part2 
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi) 
library(hgu133plus2.db)

#part 3
data<-ReadAffy(celfile.path = '/projectnb/bf528/users/lush_2022/project_1/samples')
normal<-rma(data)

#part 4
pset<-fitPLM(data,normalize = TRUE,background=TRUE)

nuse<-NUSE(pset,type = "stats")
nuse_medians<-nuse[1,]
jpeg(file = 'nuse_medians.jpeg')
hist(nuse_medians, main ="NUSE Medians", c="purple", xlab = "Medians(NUSE)")
dev.off()

rle<-RLE(pset,type = "stats")
rle_medians<-rle[1,]
jpeg(file = 'rle_medians.jpeg')
hist(rle_medians, main = "RLE Medians", c="pink", xlab = "Medians(RLE)")
dev.off()

#part 5
getwd()
cobat_data <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")

normbatch<-as.numeric(factor(cobat_data$normalizationcombatbatch))
normmod<-as.numeric(factor(cobat_data$normalizationcombatmod))

ComBat(normal, batch=normbatch, mod=normmod)

df<-exprs(normal)

write.csv(df,"/projectnb/bf528/users/lush_2022/project_1/expression_data.csv",row.names=TRUE)

#part 6
tscaled<-scale(t(df))
scaled <- t(tscaled)
principal<-prcomp(scaled,center= FALSE, scale. = FALSE)
principal$rotation
jpeg(file='pca_plot.jpeg')
plot(principal$rotation[,1], principal$rotation[,2], main = "PCA Plot", xlab="PC1 (14.5%)", ylab= "PC2 (9.5%)")
dev.off()
summary(principal)


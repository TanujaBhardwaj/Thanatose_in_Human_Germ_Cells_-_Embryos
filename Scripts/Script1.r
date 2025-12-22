# getting all the datasets

library(patchwork)
library(Seurat)
library(tidyverse)
library(ggplot2)

dataset1 <- Read10X(data.dir = "Spermatogonia/")
dataset2 <- Read10X(data.dir = "Spermatocytes/")
dataset3 <- Read10X(data.dir = "Spermatids/")
dataset4 <- Read10X(data.dir = "Spermatogenesis/")
dataset5 <- readRDS("Spermatogenesis2.RDS")
dataset6 <- read.table("spermatogonia_samples/GSE153947_RAW/counts_1.tsv.gz", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dataset6 <- as.matrix(dataset6)
dataset7 <- read.table("spermatogonia_samples/GSE153947_RAW/counts_2.tsv.gz", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dataset7 <- as.matrix(dataset7)
dataset8 <- read.table("spermatogonia_samples/GSE153947_RAW/counts_3.tsv.gz", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dataset8 <- as.matrix(dataset8)


dim(dataset1)   # 33694 11104
dim(dataset2)   # 33694  4884
dim(dataset3)   # 33694  7434
dim(dataset4)   # 33694  7134
dim(dataset5)   # 21761  3046
dim(dataset6)   # 56981  5437
dim(dataset7)   # 56981  5438
dim(dataset8)   # 56981  4672

head(dataset5@meta.data,5)


dataset1 <- CreateSeuratObject(counts = dataset1, project = "spermatogonia", min.cells = 3, min.features = 200)
dataset2 <- CreateSeuratObject(counts = dataset2, project = "Spermatocytes", min.cells = 3, min.features = 200)
dataset3 <- CreateSeuratObject(counts = dataset3, project = "Spermatids", min.cells = 3, min.features = 200)
dataset4 <- CreateSeuratObject(counts = dataset4 , project = "Spermatogenesis1", min.cells = 3, min.features = 200)
dataset6 <- CreateSeuratObject(counts = dataset6 , project = "Spermatogenesis3", min.cells = 3, min.features = 200)
dataset7 <- CreateSeuratObject(counts = dataset7 , project = "Spermatogenesis4", min.cells = 3, min.features = 200)
dataset8 <- CreateSeuratObject(counts = dataset8 , project = "Spermatogenesis5", min.cells = 3, min.features = 200)


dim(dataset1)   # 28320 11104
dim(dataset2)   # 27405  4884
dim(dataset3)   # 27898  7434
dim(dataset4)   # 29119  7132
dim(dataset5)   # 21761  3046
dim(dataset6)   # 30008  5437
dim(dataset7)   # 32876  5437
dim(dataset8)   # 32357  4672


dataset1[["percent.mt"]] <- PercentageFeatureSet(dataset1, pattern = "^MT-")
dataset2[["percent.mt"]] <- PercentageFeatureSet(dataset2, pattern = "^MT-")
dataset3[["percent.mt"]] <- PercentageFeatureSet(dataset3, pattern = "^MT-")
dataset4[["percent.mt"]] <- PercentageFeatureSet(dataset4, pattern = "^MT-")
dataset5[["percent.mt"]] <- PercentageFeatureSet(dataset5, pattern = "^MT-")
dataset6[["percent.mt"]] <- PercentageFeatureSet(dataset6, pattern = "^MT-")
dataset7[["percent.mt"]] <- PercentageFeatureSet(dataset7, pattern = "^MT-")
dataset8[["percent.mt"]] <- PercentageFeatureSet(dataset8, pattern = "^MT-")

head(dataset1@meta.data, 5)
head(dataset2@meta.data, 5)
head(dataset3@meta.data, 5)
head(dataset4@meta.data, 5)
head(dataset5@meta.data, 5)
head(dataset6@meta.data, 5)
head(dataset7@meta.data, 5)
head(dataset8@meta.data, 5)


pdf("./quality_vlnplot.pdf", width=20)
lapply(c(dataset1,dataset2,dataset3,dataset4,dataset5,dataset6,dataset7,dataset8),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

VlnPlot(dataset1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(dataset2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(dataset3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(dataset4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(dataset5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(dataset6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(dataset7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(dataset8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



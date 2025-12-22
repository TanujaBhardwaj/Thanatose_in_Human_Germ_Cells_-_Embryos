######################################     ***Normalization***    ######################################
dim(dataset1)   # 28320 10162
dim(dataset2)   # 27405  4628
dim(dataset3)   # 27898  7029
dim(dataset4)   # 29119  6665
dim(dataset5)   # 21761  2893
dim(dataset6)   # 30008  5162
dim(dataset7)   # 32876  5148
dim(dataset8)   # 32357  4340


all.list <- c(dataset1,dataset2,dataset3, dataset4,dataset5,dataset6,dataset7,dataset8)
names(all.list) <- c("dataset1", "dataset2","dataset3", "dataset4", "dataset5", "dataset6", "dataset7", "dataset8")

for (i in names(all.list)) {
  all.list[[i]] <- NormalizeData(all.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}

######################################    ***Variable***   ############################################

for (i in names(all.list)) {
  all.list[[i]] <- FindVariableFeatures(all.list[[i]], selection.method = "vst", nfeatures = 2000)
}

######################################   ***INTEGRATION***  ###########################################

n=35
all.anchors <- FindIntegrationAnchors(object.list=all.list,dims = 1:n)
all.integrated <- IntegrateData(anchorset = all.anchors, dims = 1:n)
dim(all.integrated ) # 2000 46027


######################################Saving###############################################
saveRDS(all.integrated, file = "integratedData.RDS")
#save.image(file = "integratedData.RData")
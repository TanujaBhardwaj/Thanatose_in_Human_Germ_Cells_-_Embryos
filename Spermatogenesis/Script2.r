###   ******************** Filtering cells **************************  #####

#  dataset1 ************
lb <- quantile(dataset1[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02)  # 1742
ub <- quantile(dataset1[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97)  # 5739
dataset1 <- dataset1[, dataset1[["nFeature_RNA"]] > lb & 
                       dataset1[["nFeature_RNA"]] < ub & 
                       dataset1[["percent.mt"]] < 15]

dim(dataset1) # 28320 10162
VlnPlot(dataset1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#  dataset2 ************
lb <- quantile(dataset2[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02) # 2950
ub <- quantile(dataset2[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97) # 8461
dataset2 <- dataset2[, dataset2[["nFeature_RNA"]] > lb & 
                       dataset2[["nFeature_RNA"]] < ub & 
                       dataset2[["percent.mt"]] < 15]

dim(dataset2) # 27405  4628
VlnPlot(dataset2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#  dataset3 ************
lb <- quantile(dataset3[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02) # 1710
ub <- quantile(dataset3[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97) # 7166
dataset3 <- dataset3[, dataset3[["nFeature_RNA"]] > lb & 
                       dataset3[["nFeature_RNA"]] < ub & 
                       dataset3[["percent.mt"]] < 15]

dim(dataset3) # 27898  7029
VlnPlot(dataset3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#  dataset4 ************
lb <- quantile(dataset4[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02) # 1262
ub <- quantile(dataset4[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97) # 8482
dataset4 <- dataset4[, dataset4[["nFeature_RNA"]] > lb & 
                       dataset4[["nFeature_RNA"]] < ub & 
                       dataset4[["percent.mt"]] < 15]

dim(dataset4) # 29119  6665
VlnPlot(dataset4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#  dataset5 ************
lb <- quantile(dataset5[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02) # 975
ub <- quantile(dataset5[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97) # 10882
dataset5 <- dataset5[, dataset5[["nFeature_RNA"]] > lb & 
                       dataset5[["nFeature_RNA"]] < ub & 
                       dataset5[["percent.mt"]] < 15]

dim(dataset5) # 21761  2893
VlnPlot(dataset5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#  dataset6 ************
lb <- quantile(dataset6[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02) # 1347
ub <- quantile(dataset6[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97) # 5137
dataset6 <- dataset6[, dataset6[["nFeature_RNA"]] > lb & 
                       dataset6[["nFeature_RNA"]] < ub & 
                       dataset6[["percent.mt"]] < 15]

dim(dataset6) # 30008  5162
VlnPlot(dataset6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#  dataset7 ************
lb <- quantile(dataset7[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02) # 1399
ub <- quantile(dataset7[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97) # 6932
dataset7 <- dataset7[, dataset7[["nFeature_RNA"]] > lb & 
                       dataset7[["nFeature_RNA"]] < ub & 
                       dataset7[["percent.mt"]] < 15]

dim(dataset7) # 32876  5148
VlnPlot(dataset7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



#  dataset8 ************
lb <- quantile(dataset8[["nFeature_RNA"]]$nFeature_RNA, probs = 0.02) # 1426
ub <- quantile(dataset8[["nFeature_RNA"]]$nFeature_RNA, probs = 0.97) # 7070
dataset8 <- dataset8[, dataset8[["nFeature_RNA"]] > lb & 
                       dataset8[["nFeature_RNA"]] < ub & 
                       dataset8[["percent.mt"]] < 15]

dim(dataset8) # 32357  4340
VlnPlot(dataset8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


pdf("./quality_vlnplot_filtered.pdf", width=20)
lapply(c(dataset1,dataset2,dataset3,dataset4,dataset5,dataset6,dataset7,dataset8),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()





# *************************** ANALYSIS ******************************************** #
library(Seurat)
library(SeuratObject)

#all.integrated <- readRDS("integratedData.RDS")

######################################Dimension Reduction & Clustering##########################################

head(all.integrated@meta.data, 5)
all.integrated@meta.data$nCount_RNA


###   ****DR
DefaultAssay(all.integrated) <- "integrated"
all.integrated <- ScaleData(all.integrated, verbose = FALSE)
all.integrated <- RunPCA(all.integrated, npcs = 50, verbose = FALSE)

Elbowplot <- ElbowPlot(all.integrated, ndims = 35, reduction = "pca")

all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:30)
dim(all.integrated) # 2000 46027

DimPlot(all.integrated, reduction = 'umap', group.by = 'orig.ident')
DimPlot(all.integrated, reduction = 'umap', split.by = 'orig.ident')


###Clustering
all.integrated <- FindNeighbors(all.integrated, dims = 1:30)
all.integrated <- FindClusters(all.integrated, resolution = 0.3)

table(Idents(all.integrated))
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 5869 4850 4158 3723 3583 3124 3021 2766 2761 2686 2309 2085 1972 1561  992  567 


pdf("./integrated_UMAP.pdf", width=20)
DimPlot(all.integrated, reduction = 'umap', group.by = 'orig.ident')
DimPlot(all.integrated, reduction = "umap", group.by = 'ident')
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
DimPlot(all.integrated, reduction = 'umap', split.by ='orig.ident')
dev.off()


head(all.integrated@meta.data,5)
table(all.integrated@meta.data$orig.ident)

# Spermatids    Spermatocytes Spermatogenesis1 Spermatogenesis2 Spermatogenesis3 Spermatogenesis4 
#   7029             4628             6665             2893             5162             5148 
# Spermatogenesis5    spermatogonia 
#     4340            10162 



DefaultAssay(all.integrated) <- "RNA"


pdf("./germ_cells.pdf", width=10)
FeaturePlot(all.integrated, features = c("DDX4"),cols=c("lightgrey","darkviolet")) + RotatedAxis()
VlnPlot(all.integrated, c("DDX4"))
DotPlot(all.integrated, features =  c("DDX4"),cols=c("lightgrey","darkviolet")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


##   *************************************************  MARKER  *********************************   ###

pdf("./Markers_spermatogenesis.pdf", width=15)
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
DotPlot(all.integrated, features = c("CYP26B1","INSL3","MYH11","ALDH1A1","CD68","CD163" ),
        cols=c("lightgrey","black")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))  +
  ggtitle("Somatic cell markers") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
FeaturePlot(all.integrated, features = c("CYP26B1","INSL3","MYH11","ALDH1A1","CD68","CD163" ),
        cols=c("lightgrey","darkviolet")) + RotatedAxis()
DotPlot(all.integrated, features = c("NANOS2","PIWIL4","GFRA1","SALL4","MAGEA4","HMGA1" ),
        cols=c("lightgrey","black")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))  +
  ggtitle("Spermatogonia cell markers") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
FeaturePlot(all.integrated, features = c("NANOS2","PIWIL4","GFRA1","SALL4","MAGEA4","HMGA1" ),
        cols=c("lightgrey","darkviolet")) + RotatedAxis()
DotPlot(all.integrated, features = c("DMC1","RAD51AP2","PIWIL1","SYCP3","OVOL2","ACR", "PGK2" ),
        cols=c("lightgrey","black")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))  +
  ggtitle("Spermatocyte cell markers") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
FeaturePlot(all.integrated, features = c("DMC1","RAD51AP2","PIWIL1","SYCP3","OVOL2","ACR", "PGK2" ),
        cols=c("lightgrey","darkviolet")) 
DotPlot(all.integrated, features = c("TEX29","SUN5","SPEM1","SYCP3","OVOL2","ACR", "PGK2"),
        cols=c("lightgrey","black")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))  +
  ggtitle("Spermatid cell markers") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
FeaturePlot(all.integrated, features = c("TEX29","SUN5","SPEM1","SYCP3","OVOL2","ACR", "PGK2"),
        cols=c("lightgrey","darkviolet")) 
dev.off()


##  ************************************
pdf("./THAP9_in_cells.pdf", width=10)
FeaturePlot(all.integrated, features = c("THAP9"),cols=c("lightgrey","darkviolet")) + RotatedAxis()
VlnPlot(all.integrated, c("THAP9"))
DotPlot(all.integrated, features =  c("THAP9"),cols=c("lightgrey","darkviolet")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


#################            **********DEGs***********            #########################

for_DE <- all.integrated
for_DE <- JoinLayers(for_DE)  # 49139 features across 46027 samples
DefaultAssay(for_DE) <- "RNA"
all.markers <- FindAllMarkers(for_DE, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./DEGs_integrated.csv")
#all_markers <- read.csv("DEGs_integrated.csv")

library(openxlsx)
write.xlsx(all.markers, 
           file = "./DEGs_integrated_spermatogenesis.xlsx", 
           rowNames = TRUE)





"THAP9" %in% rownames(for_DE)
for_DE["THAP9", ]

subset(all_markers, gene == "THAP9")

grep("THAP9", rownames(for_DE), value = TRUE)

write.table(subset(all.markers, gene == "THAP9"),
            file = "THAP9_makers_spermatogenesis.txt",
            sep = "\t",
            quote = FALSE)



#####   ****** Heatmap ******  ######
top_markers <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC))

genes.use <- unique(top_markers$gene)
for_DE <- ScaleData(for_DE, features = unique(top_markers$gene), assay = "RNA")

pdf("./heatmap_spermatogonia.pdf", width=20, height = 15)
DoHeatmap(
  for_DE,
  features = unique(top_markers$gene),
  assay = "RNA",
  slot = "scale.data"
)
dev.off()




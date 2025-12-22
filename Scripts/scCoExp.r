#  --------------------------------------------------------------------  #
       #   ****************** LOAD libraries ****************   #
#  --------------------------------------------------------------------  #
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

#  --------------------------------------------------------------------  #
theme_set(theme_cowplot()) #cowplot theme
set.seed(12345) #random seed for reproducibility
enableWGCNAThreads(nThreads = 18)
#enableWGCNAThreads()

#  --------------------------------------------------------------------  #
#                           load the dataset
#  --------------------------------------------------------------------  #
seurat_obj <- readRDS("labelled_integrated.RDS")
seurat_obj <- JoinLayers(seurat_obj)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)
view(seurat_obj@meta.data)
head(seurat_obj@meta.data)

seurat_obj$cell_type <- Idents(seurat_obj)
head(seurat_obj@meta.data[, c("seurat_clusters", "cell_type")])
seurat_obj$sample <- seurat_obj@meta.data$orig.ident
head(seurat_obj@meta.data[, c("seurat_clusters", "cell_type", "sample")])

p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
  umap_theme() + ggtitle('spermatogenesis') + NoLegend()
p

#  --------------------------------------------------------------------  #
#              Setting up the WGCNA
#  --------------------------------------------------------------------  #
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, 
  wgcna_name = "spermatogenesis" # the name of the experiment
)

grep("THAP9", rownames(seurat_obj), value = TRUE)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = "cell_type", # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

grep("THAP9", rownames(seurat_obj), value = TRUE)



#  --------------------------------------------------------------------  #
# Process the Metacell Seurat Object
metacell_obj <- GetMetacellObject(seurat_obj)

seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='orig.ident')
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='pca', dims=1:30)
#seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:30)
#  --------------------------------------------------------------------  #

p1 <- DimPlotMetacells(seurat_obj, group.by='cell_type') + umap_theme() + ggtitle("Cell Type")
p2 <- DimPlotMetacells(seurat_obj, group.by="orig.ident") + umap_theme() + ggtitle("sample")
p1 | p2

#  --------------------------------------------------------------------  #

#  --------------------------------------------------------------------  #
#####   *****Co-expression network analysis****   #####
#  --------------------------------------------------------------------  #

#  --------------------------------------------------------------------  #

# Set up the expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("round SPT1","round SPT2", "elongated STD"), # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

#  --------------------------------------------------------------------  #
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)


plot_list <- PlotSoftPowers(seurat_obj) # plot the results
wrap_plots(plot_list, ncol=2) # assemble with patchwork


# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'spermatid')  # Using soft_power = 10
PlotDendrogram(seurat_obj, main='spermatids hdWGCNA Dendrogram')

#  --------------------------------------------------------------------  #
#  --------------------------------------------------------------------  #
TOM <- GetTOM(seurat_obj)
saveRDS(TOM, "TOM_spermatids.rds")
#  --------------------------------------------------------------------  #
#  --------------------------------------------------------------------  #


# Module Eigengenes and Connectivity
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj)) # need to run ScaleData first or else harmony throws an error

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)

hMEs <- GetMEs(seurat_obj) # harmonized module eigengenes
MEs <- GetMEs(seurat_obj, harmonized=FALSE) # module eigengenes
names(seurat_obj@misc$hdWGCNA$MEs) # compute eigengene-based connectivity (kME)


seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "cell_type",
  group_name = c("round SPT1","round SPT2", "elongated STD")
)
  
#  --------------------------------------------------------------------  #
#  --------------------------------------------------------------------  #

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "spermatids-M") # rename the modules

p <- PlotKMEs(seurat_obj, ncol=5) # plot genes ranked by kME for each module
p

# get the module assignment table
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
head(modules[,1:6])

#  --------------------------------------------------------------------  #
saveRDS(modules, file='Modules_hdWGCNA_Spermatids.rds')
#  --------------------------------------------------------------------  #


hub_df <- GetHubGenes(seurat_obj, n_hubs = 10) # get hub genes
head(hub_df)

#  --------------------------------------------------------------------  #
saveRDS(seurat_obj, file='hdWGCNA_Spermatids_object.rds')
#  --------------------------------------------------------------------  #


#  --------------------------------------------------------------------  #
###  *************for THAP9 gene***********************   ####
#  --------------------------------------------------------------------  #

modules[rownames(modules) == "THAP9", ]  # Check which module THAP9 belongs to
target_module <- modules$module[modules$gene == "THAP9"]

modules %>%
  filter(module == target_module) %>%
  arrange(desc(!!sym(paste0("kME_", target_module)))) %>%
  head(20)
modules[modules$gene == "THAP9", ]

#  --------------------------------------------------------------------  #
saveRDS(modules, "modules_spermatids_withTHAP9.RDS")
write.csv(modules, "modules_hdWGCNA_withTHAP9.csv")
#  --------------------------------------------------------------------  #

# Computing hub gene signature scores
library(UCell)

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

#featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
wrap_plots(plot_list, ncol=6)
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE)
wrap_plots(plot_list, ncol=6)

seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$sample), ' '))[,1]

#  --------------------------------------------------------------------  #
library(ggradar)
ModuleRadarPlot(
  seurat_obj,
  group.by = 'cluster',
  barcodes = seurat_obj@meta.data %>% subset(cell_type == 'elongated STD') %>% rownames(),
  axis.label.size=4,
  grid.label.size=4)
#  --------------------------------------------------------------------  #
#  --------------------------------------------------------------------  #

library(corrplot)
ModuleCorrelogram(seurat_obj) # plot module correlagram
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']


seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
p <- DotPlot(seurat_obj, features=mods, group.by = 'cell_type')
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
p

#  --------------------------------------------------------------------  #
#  **************** Saving files ****************  #
#  --------------------------------------------------------------------  #
saveRDS(MEs, "MEs.rds")
saveRDS(hMEs, "hMEs.rds")
saveRDS(TOM_M7, "TOM_M7.rds")
#  --------------------------------------------------------------------  #
#  --------------------------------------------------------------------  #

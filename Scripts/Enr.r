# ----------------------------------------------------- #
#          ************Libraries ***********
# ----------------------------------------------------- #

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)

# ----------------------------------------------------- #
theme_set(theme_cowplot())
set.seed(12345)
# ----------------------------------------------------- #


# ----------------------------------------------------- #
# enrichr databases
dbs <- c('GO_Biological_Process_2025','GO_Cellular_Component_2025','GO_Molecular_Function_2025')
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs,
  max_genes = Inf # use max_genes = Inf to choose all genes
)

enrich_df <- GetEnrichrTable(seurat_obj)
head(enrich_df)

# ----------------------------------------------------- #


# ----------------------------------------------------- #
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 20, # number of enriched terms to show
  plot_size = c(5,7), 
  logscale=TRUE 
)

EnrichrDotPlot(
  seurat_obj,
  mods = "all",
  database = "GO_Biological_Process_2025", 
  n_terms=2, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = TRUE 
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

EnrichrDotPlot(
  seurat_obj,
  mods = "spermatids-M7", 
  database = "GO_Biological_Process_2025", 
  n_terms=2, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = TRUE # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))


EnrichrDotPlot(
  seurat_obj,
  mods = "spermatids-M7",
  database = "GO_Cellular_Component_2025", 
  n_terms=10, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = TRUE 
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

# ----------------------------------------------------- #
# ----------------------------------------------------- #


# ----------------------------------------------------- #
#   ***Gene Set Enrichment Analysis (GSEA)***
# ----------------------------------------------------- #
library(fgsea)
pathways <- fgsea::gmtPathways('GO_Biological_Process_2025')
names(pathways) <- stringr::str_replace(names(pathways), " \\s*\\([^\\)]+\\)", "")
modules <- GetModules(seurat_obj) %>% subset(module != 'grey') # get the modules table and remove grey genes
# ----------------------------------------------------- #

###  **** using only the genes in one module ****
# rank by spermatids-M7 genes only by kME
cur_mod <- 'spermatids-M7'
modules <- GetModules(seurat_obj) %>% subset(module == cur_mod)
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]

gsea_df2 <- fgsea::fgsea(
  pathways = pathways, 
  stats = ranks,
  minSize = 3,
  maxSize = 500
)
n_signif1 <- subset(gsea_df, pval < 0.05) %>% nrow 
n_signif2 <- subset(gsea_df2, pval < 0.05) %>% nrow
print(paste0("All genes: ", n_signif1, ", ", cur_mod, " only: ", n_signif2))

top_pathways <- gsea_df2 %>% 
  subset(pval < 0.05) %>% 
  slice_max(order_by=NES, n=25) %>% 
  .$pathway

plotGseaTable(
  pathways[top_pathways], 
  ranks, 
  gsea_df2, 
  gseaParam=0.5,
  colwidths = c(10, 4, 1, 1, 1)
)

# ----------------------------------------------------- #
# ----------------------------------------------------- #
# rank by spermatids-M7 genes only by kME
cur_mod <- 'spermatids-M7'
modules <- GetModules(seurat_obj) %>% subset(module == cur_mod)
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]

gsea_df3 <- fgsea::fgsea(
  pathways = pathways, 
  stats = ranks,
  minSize = 3,
  maxSize = 500
)

n_signif3 <- subset(gsea_df3, pval < 0.05) %>% nrow 
print(n_signif3)
# ----------------------------------------------------- #

# ------------------------------------------------------------------------  #

    # ************ Libraries for trajectory analysis  *********** #

# ------------------------------------------------------------------------  #

library(monocle3)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(SeuratWrappers)

# ------------------------------------------------------------------------  #
                 #  ************ Load data **********  #
# ------------------------------------------------------------------------  #
data <- readRDS("labelled_integrated.RDS")

DefaultAssay(data) <- "RNA"
DimPlot(data, reduction = 'umap', group.by = 'orig.ident')
table(Idents(data))



# ------------------------------------------------------------------------  #
#############      *******MONOCLE3 cds********        ####################
# ------------------------------------------------------------------------  #

DefaultAssay(data) <- "RNA"
data <- JoinLayers(object = data, assays = "RNA")
cds <- SeuratWrappers::as.cell_data_set(data) #change to cds here
cds <- cluster_cells(cds)


colData(cds)  ## Gives information about the cells metadata
fData(cds)  ## Gene metadata
rownames(cds)  ### gene names

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name


# ------------------------------------------------------------------------  #

###########  ****RETRIVING THE CLUSTER INFORMATION***   ###########

# ------------------------------------------------------------------------  #

##   ***********assign partitions
recreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)


cds@clusters$UMAP$partitions <- recreate.partition

###   ****assigning the cluster information
list_cluster <- data@active.ident
cds@clusters$UMAP$clusters <- list_cluster


### ****assigning UMAP coordinate 
cds@int_colData@listData$reducedDims$UMAP <- data@reductions$umap@cell.embeddings

# ------------------------------------------------------------------------  #
#                     *********plot***************
# ------------------------------------------------------------------------  #
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = 'right')


# ------------------------------------------------------------------------  #
#########    *** STEPS ***   ###########
# ------------------------------------------------------------------------  #

plot_cells(cds, show_trajectory_graph = FALSE)

plot_cells(cds, show_trajectory_graph = FALSE,
           color_cells_by = "partition")

cds <- learn_graph(cds, use_partition = FALSE) # graph learned acress all partitions
## Use TRUE unless you know your partitions should be together

a1 <- plot_cells(cds,
                 color_cells_by = "cluster",
                 label_groups_by_cluster = FALSE,
                 label_branch_points = FALSE,
                 label_roots = FALSE,
                 label_leaves = FALSE,
                 group_label_size = 5)



cds <- order_cells(cds,reduction_method = 'UMAP')

a2 <-plot_cells(cds, color_cells_by = "pseudotime",
                label_groups_by_cluster = FALSE,
                label_branch_points = FALSE,
                label_roots = FALSE,
                label_leaves = FALSE,
                group_label_size = 5)

a1 + a2

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=3)


saveRDS(cds, "spermatogenesis_cds.RDS")

# ------------------------------------------------------------------------  #
# ------------------------------------------------------------------------  #

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, ident, fill = ident)) +
  geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime, median), fill = ident)) +
  geom_boxplot()


###  ****Finding genes that change as function of pseudotime  ***********
cds_pt_res <- graph_test(cds, neighbor_graph = "principal_graph")
cds_pt_res

cds_pt_res %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head()



####   **************     *******pseudotime in seurat clusters******   *************   ####
data$pseudotime <- pseudotime(cds)

FeaturePlot(data, 
            features = "pseudotime", 
            label = T,
            label.size = 4,
            label.color = 'black')

####   **************     *************   *************   ####

#cds_pt_res <- na.omit(cds_pt_res)
#cds_pt_res <- cds_pt_res[cds_pt_res$p_value < 0.05 & cds_pt_res$status == "OK",]
#cds_pt_res

cds_pt_res[order(-cds_pt_res$morans_test_statistic),]

plot_cells(cds,
           label_cell_groups = FALSE,
           show_trajectory_graph = TRUE,
           label_branch_points=TRUE,
           graph_label_size=3)




########## Subsetting one cell population ***************
cds_subset <- choose_cells(cds)
cds_subset_pt_res <- graph_test(cds_subset, neighbor_graph = "principal_graph")
cds_subset_pt_res <- na.omit(cds_subset_pt_res)
cds_subset_pt_res <- cds_subset_pt_res[cds_subset_pt_res$p_value < 0.05 & cds_subset_pt_res$status == "OK", ]
cds_subset_pt_res




plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by="monocle3_pseudotime",
                         min_expr=0.5)


plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by="ident",
                         min_expr=0.5)





################   *****************************************    ################
cds_3d <- preprocess_cds(cds, num_dim = 50, method = "PCA")
cds_3d <- reduce_dimension(cds_3d, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)

cds_3d <- order_cells(cds_3d, root_pr_nodes = get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")

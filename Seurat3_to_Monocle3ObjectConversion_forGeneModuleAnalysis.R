# Load all the required Libraries
library(Biobase)
library(Seurat)
library(plotly)
library(dplyr)
library(monocle3)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)

# Setting the default directory
setwd('/Users/sha6hg/Desktop/IPF_scRNA/')

########################################################### Running Monocle3 ############################################################ 
#  Create Monocle3 CDS object to get gene modules for all clusters
load("Mouse_all_integrated.Robj")
cds <- as.cell_data_set(integrated_all_mouse)
save(cds, file = "Mouse_all_integrated_Monocle3CDS.Robj")
load('Mouse_all_integrated_Monocle3_Step1.Robj')

# Retrieve Counts from cds object 
counts <- counts(cds)

# Read cell Meta
cell_meta <- as.matrix(read.table("cell_meta.txt", sep="\t", row.names=1, header=T))
cell_meta <- data.frame(cell_meta)

# Read gene Meta
genes_meta <- as.matrix(read.table("genes.txt", sep="\t", row.names=1, header=T))
genes_meta <- data.frame(genes_meta)

# Create new Monocle object
cds_Step2 <- new_cell_data_set(expression_data = counts, cell_metadata = cell_meta, gene_metadata = genes_meta)

# Remove genes and cells with no expression
cds_Step2 <- cds_Step2[,Matrix::rowSums(exprs(cds_Step2)) != 0]
cds_Step2 <- cds_Step2[,Matrix::colSums(exprs(cds_Step2)) != 0]

# Estimate size factors
cds_Step2 <- estimate_size_factors(cds_Step2)

# Preprocess monocle object
cds_Step2 <- preprocess_cds(cds_Step2, dim=40, method = "PCA")

# Reduce the dimensions using UMAP
cds_Step2 <- align_cds(cds_Step2)
cds_Step2 <- reduce_dimension(cds_Step2, reduction_method = "UMAP")

# Cluster the cells
cds_Step2 <- cluster_cells(cds_Step2, resolution=1e-5)

plot_cells(cds_Step2, label_groups_by_cluster=TRUE, label_leaves=FALSE, label_branch_points=FALSE)
#cds_Step2 <- learn_graph(cds_Step2)

pr_graph_test_res <- graph_test(cds_Step2, neighbor_graph="knn")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(cds_Step2[pr_deg_ids,], resolution=1e-2)
#cell_group_df <- tibble::tibble(cell=row.names(colData(cds_Step2)), cell_group=partitions(cds)[colnames(neurons_cds)])
#cell_group_df <- tibble::tibble(cell=row.names(colData(cds_Step2)), cell_group=cds_Step2@colData$integrated_snn_res.0.3)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_Step2)), cell_group=cds_Step2@colData$MouseGeneral_Annotation)
agg_mat <- aggregate_gene_expression(cds_Step2, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Clusters ", colnames(agg_mat))

write.table(agg_mat, file = "Mouse_all_integrated_Heatmap_byCluster_res0.3.txt", row.names = T, col.names = T, quote = F, sep = "\t")
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE, scale="column", clustering_method="ward.D2", fontsize=6)

# Save the final Monocle3 CDS object
save(cds_Step2, file="Mouse_all_integrated_Monocle3_Step2.Robj")
########################################################################################################################################

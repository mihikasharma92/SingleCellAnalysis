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

####################################################### Pre Processing #######################################################
# Get a list of all the folders to be processed
path_main = "/Users/sha6hg/Desktop/IPF_scRNA/"

# List of MOUSE datasets
mouse_GSE = c("GSE111664", "GSE121611", "GSE124872", "GSE127803", "GSE129605")

# Create a list containing path to each GSE dataset conatining Raw Counts Matrix
mouse_GSE_Path <- paste(path_main, mouse_GSE, "/Seurat", sep="")

for(GSE_dir in  mouse_GSE_Path)
{
  GSE_id <- unlist(strsplit(GSE_dir, "/"))[6]
  # Get a list of all GSM folders in the GSE folder 
  iter_file = 1
  list_obj_to_integrate = list()
  all_dir_list = list.dirs(GSE_dir, recursive = F)
  #print(all_dir_list)
  
  # Iterate through each GSE folder and add the Robj of each GSM to a list
  for(directory in all_dir_list)
  {
    all_file_list = list.files(directory, recursive = F)
    obj <- grep('_Stage2.Robj', all_file_list, value=TRUE)
    obj_path = paste0(directory,"/", obj)
    #print(obj_path)
    list_obj_to_integrate[iter_file] = obj_path
    iter_file = iter_file + 1
  }
  print1 <- paste0("Path for all the objects added to the list for: ", GSE_id)
  print(print1)
  
  # Iterate through path to load all the Robj for one GSE
  list1 <- unlist(list_obj_to_integrate)
  for(object in list1)
  {
    load(object)
  }
  print2 <- paste0("All objects loaded for: ", GSE_id)
  print(print2)
  
  # Get a list of all the Seurat Objects in the environment for a particular GSE
  obj_list <- objects()
  iter2 = 1
  list_obj_to_integrate1 = list()
  for(obj in obj_list)
  {
    if(class(get(obj))[1] == "Seurat")
    {
      obj1 = obj
      list_obj_to_integrate1[iter2] = get(obj)
      iter2 = iter2 + 1
    }
  }  
  print3 <- paste0("Obtained the list of all the objects to integrate")
  print(print3)
  integrated.anchors <- FindIntegrationAnchors(object.list = list_obj_to_integrate1, dims = 1:30)
  print4 <- paste0("Integration anchors found for all objects for: ", GSE_id)
  print(print4)
  # Integrate all Seurat Objects
  integrated <- IntegrateData(anchorset = integrated.anchors, dims = 1:20)
  DefaultAssay(integrated) <- "integrated"
  # Run the standard workflow for visualization and clustering
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
  # t-SNE and Clustering
  integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:20, check_duplicates = F)
  integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
  integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:20)
  # Calculate clusters at multiple resolutions for seurat object
  resolution_list <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2)
  for (i in resolution_list)
  {
    integrated <- FindClusters(integrated, resolution = i, reduction = "tsne")
  }
  print5 <- paste0("Integration completed for: ", GSE_id)
  print(print5)
  # Save integrated object
  filename = paste0(unlist(strsplit(GSE_dir, "/"))[6], ".Robj")
  save(integrated, file = filename)
  
  print6 <- paste0("Seurat R Object saved for: ", GSE_id)
  print(print6)
  
  # Clear the environment of all Seurat objects 
  print("Clearing the environment for new GSE integration")
  obj_list <- objects()
  for(obj in obj_list)
  {
    if(class(get(obj))[1] == "Seurat")
    {
      print(obj)
      rm(list = obj, envir = environment())
    }
  }
}


# Iterate over the folder to load all the R objects for each GSE 
all_file_list = list.files(path_main, recursive = F)
Robj_file_list <- paste(path_main, all_file_list, sep="")
for(file in Robj_file_list)
{
  if(length(grep('.Robj', file, value=TRUE)) > 0)
  {
    load(file)
  }  
}

# Create a list of all the loaded R objects
obj_list1 <- objects()
iter3 = 1
list_obj_to_integrate2 = list()
for(obj in obj_list1)
{
  if(class(get(obj))[1] == "Seurat")
  {
    list_obj_to_integrate2[iter3] = get(obj)
    iter3 = iter3 + 1
  }
}  

# Integrate all Seurat Objects
integrated.anchors <- FindIntegrationAnchors(object.list = list_obj_to_integrate2, dims = 1:30)
integrated_all_mouse <- IntegrateData(anchorset = integrated.anchors, dims = 1:20)
DefaultAssay(integrated_all_mouse) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated_all_mouse <- ScaleData(integrated_all_mouse, verbose = FALSE)
integrated_all_mouse <- RunPCA(integrated_all_mouse, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated_all_mouse <- RunTSNE(integrated_all_mouse, reduction = "pca", dims = 1:20, check_duplicates = F)
integrated_all_mouse <- RunUMAP(integrated_all_mouse, reduction = "pca", dims = 1:20)
integrated_all_mouse <- FindNeighbors(integrated_all_mouse, reduction = "pca", dims = 1:20)
# Calculate clusters at multiple resolutions for seurat object
resolution_list <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2)
for (i in resolution_list)
{
  integrated_all_mouse <- FindClusters(integrated_all_mouse, resolution = i, reduction = "umap")
}

# Save the R object
save(integrated_all_mouse, file = "Mouse_all_integrated.Robj")

# Save the metadata file for adding additional annotations
metadata <- as.matrix(integrated_all_mouse@meta.data)
write.table(metadata, file = "Integrated_all_mouse_Metadata.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

# Re-upload the altered metadata file to add new annotations
metafile1 <- "Integrated_all_mouse_Metadata.txt"
col1.data <- as.matrix(read.table(file = metafile1, sep = "\t", header=T, row.names = 1, check.names = FALSE))
col1_df <- data.frame(col1.data)
col1_df1 <- new("AnnotatedDataFrame", data = col1_df)
rownames(col1_df1) <- col1_df1$Cells

metadata <- integrated_all_mouse@meta.data
head(metadata)
metadata["Sample_type"] <- col1_df1$Sample_type
head(metadata)
CellsMetaTrim <- subset(metadata, select = c("Sample_type"))
integrated_all_mouse <- AddMetaData(integrated_all_mouse, CellsMetaTrim)

# Save the R object
save(integrated_all_mouse, file = "Mouse_all_integrated.Robj")

load("Mouse_Integration_3_Final_AfterSingleR.Robj")

# Save the metadata file for adding additional annotations
metadata <- as.matrix(Final.integrated@meta.data)
write.table(metadata, file = "Integrated_all_mouse_Metadata1.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

# Re-upload the altered metadata file to add new annotations
metafile1 <- "Integrated_all_mouse_Metadata1.txt"
col1.data <- as.matrix(read.table(file = metafile1, sep = "\t", header=T, row.names = 1, check.names = FALSE))
col1_df <- data.frame(col1.data)
col1_df1 <- new("AnnotatedDataFrame", data = col1_df)
rownames(col1_df1) <- col1_df1$Cells

# Save the R object
save(integrated_all_mouse, file = "Mouse_all_integrated.Robj")

# Create UMAP Cluster plots for Visualization
cols1 <- c("deeppink4", "limegreen", "darkorange2", "violet", "plum4", "dodgerblue1", "darkturquoise", "coral3", "powderblue", "goldenrod1", "darkcyan", "red2", "thistle4","mediumseagreen", "blue4", "ivory", "pink", "cornflowerblue", "indianred2", "lightblue3", "darksalmon", "olivedrab2", "violet", "darkslategray2", "#8B4513")
cols2 <- c("#FF4500", "#87CEEB", "#FFD700", "#D8BFD8", "#6A5ACD", "#000000", "#FFA500", "#6B8E23")

DimPlot(integrated_all_mouse, reduction = "umap", pt.size = 0.025, group.by = "MouseImmGen_Annotation", cols = cols1) + DarkTheme() + guides(col = guide_legend(ncol=1, override.aes = list(size = 9))) + theme(legend.text=element_text(size=20))
ggsave('AllMouse_UMAP_DimPlot_byImmGenAnnotations.png', width=20, height=15, dpi= 600)

DimPlot(integrated_all_mouse, reduction = "umap", pt.size = 0.025, group.by = "MouseGeneral_Annotation", cols = cols1) + DarkTheme()+ guides(col = guide_legend(ncol=1, override.aes = list(size = 9))) + theme(legend.text=element_text(size=20))
ggsave('AllMouse_UMAP_DimPlot_byMouseGeneral_Annotation.png', width=20, height=15, dpi = 600)

DimPlot(integrated_all_mouse, reduction = "umap", pt.size = 0.025, group.by = "integrated_snn_res.0.3", cols = cols1) + DarkTheme() + guides(col = guide_legend(ncol=1, override.aes = list(size = 9))) + theme(legend.text=element_text(size=20))
ggsave('AllMouse_UMAP_DimPlot_byCluster_res0.3.png', width=20, height=15, dpi = 600)

DimPlot(integrated_all_mouse, reduction = "umap", pt.size = 0.025, group.by = "orig.ident") + DarkTheme() + guides(col = guide_legend(ncol=1, override.aes = list(size = 5))) + theme(legend.text=element_text(size=10))
ggsave('AllMouse_UMAP_DimPlot_byGSMids.png', width=20, height=15, dpi = 600)

DimPlot(integrated_all_mouse, reduction = "umap", pt.size = 0.01, group.by = "Sample_type", cols = c("Red", "Blue", "Green")) + DarkTheme() + guides(col = guide_legend(ncol=1, override.aes = list(size = 9))) + theme(legend.text=element_text(size=20))
ggsave('AllMouse_UMAP_DimPlot_bySampleType.png', width=20, height=15, dpi = 600)
########################################################################################################################################

########################################################### Create Dot Plot ############################################################ 
# Create dot plot to help annotations
dir = "/Users/sha6hg/Desktop/IPF_scRNA/AnnotationFiles/KrasnowAnnotationFiles/Human_to_Mouse_Conversion/"
output_dir = "/Users/sha6hg/Desktop/IPF_scRNA/KrasnowAnnotationResults/"
all_file_list = list.files(dir, recursive = F)[grepl(".txt", list.files(dir, recursive = F))]

# Iterate through the file list to get a list of genes
object_name = integrated_all_mouse
for(file in all_file_list)
{
  file_path <- paste(dir, file, sep = "")
  DefaultAssay(integrated_all_mouse) <- "RNA"
  file_path <- paste0(dir, file)
  gene_list1 <- as.data.frame(read.table(file_path, sep = "\n"))$V1
  gene_list2 <- markers$gene
  gene_list_final <- unique(intersect(gene_list1, gene_list2))
  #print(gene_list1)
  filename1 = paste0(output_dir, strsplit(file, ".txt")[[1]], ".png")
  DotPlot(integrated_all_mouse, group.by = "integrated_snn_res.0.3", gene_list_final, assay = "RNA") + guides(shape = guide_legend(override.aes = list(size = 0.5))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(file = filename1, width = 30, height = 15)
}
########################################################################################################################################

library(crayon)
library(backports)
library(withr)
library(rstudioapi)
library(cli)

library(Seurat)
library(Signac)

library(labeling)
library(tidyverse)
library(gridExtra)
library(hdf5)


my_save_path <- ".../integrated_data/"

## Load NR2E3 Data
load(".../cluster_annotation/organoid_all_tp.RData")
Idents(organoid_all_tp) <- "draft_celltype_v3"

## Load Kallman reprocessed data
load(".../Kallman_NRL/plots/final_object_label_transfer.RData")

### Split lists into individual runs
NRL_list <- SplitObject(IVR_object.combined, split.by = "orig.ident")
NR2E3_list <- SplitObject(organoid_all_tp, split.by = "orig.ident")

combined_list <- c(NRL_list, NR2E3_list)


for(i in 1:length(combined_list)){
  my_initial_object <- combined_list[[i]]
  DefaultAssay(my_initial_object) <- "RNA"
  my_initial_object <- NormalizeData(my_initial_object, normalization.method = "LogNormalize", scale.factor = 10000)
  my_initial_object <- FindVariableFeatures(my_initial_object, selection.method = "vst", nfeatures = 2000)
  combined_list[[i]] <- my_initial_object
}

#########################################
########## PERFORM THE CCA ##############
#########################################
IVR_object.anchors <- FindIntegrationAnchors(object.list = combined_list, dims=1:25)
IVR_object.combined <- IntegrateData(anchorset = IVR_object.anchors, dims = 1:25)

DefaultAssay(IVR_object.combined) <- "integrated"
IVR_object.combined <- ScaleData(IVR_object.combined, verbose = FALSE)
IVR_object.combined <- RunPCA(IVR_object.combined, npcs = 30, verbose = FALSE)

ElbowPlot(IVR_object.combined, ndims = 25)
ggsave(paste0(my_save_path, "elbow_plot.pdf"), height = 8, width = 8)

#4. Perform the CCA

#load(file = paste0(my_save_path, "intermediate_object.RData"))

#### READ OBJECT
final.dim = 25

IVR_object.combined <- RunUMAP(IVR_object.combined, reduction = "pca", dims = 1:final.dim)
IVR_object.combined <- FindNeighbors(IVR_object.combined, reduction = "pca", dims = 1:final.dim)
IVR_object.combined <- FindClusters(IVR_object.combined, resolution = 0.5)

save(IVR_object.combined, file = paste0(my_save_path, "nr2e3_nrl_all_timepoints_final_object.RData"))
my_labels <- c("libraries", "filtering_criteria", "dimensions")
my_list <- list(my_filepaths, my_filter_csv, final.dim)
save(my_list, file = paste0(my_save_path, "documented_criteria_for_final_object.RData"))


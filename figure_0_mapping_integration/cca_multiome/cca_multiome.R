library(crayon, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(backports, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(withr, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(rstudioapi, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(cli, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(Seurat, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(Signac, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(labeling, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(tidyverse, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(gridExtra, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(hdf5r, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

args <- commandArgs(trailingOnly = TRUE)

csv_args <- args[2]

#### 1. Read in important files

#### FIRST ARGUMENT: CSV file specifying filepaths of cellranger output and library ID for each file to me CCA-ed.
#### SECOND ARGUMENT: filepath where data (plots, itermediate objects) can be saved. RECOMMEND CREATING OUTPUT FOLDER WITH DATE.
#### THIRD ARGUMENT: filepath of filtering criteria (INCLUDE even if at step 1).
#### FOURTH ARGUMENT: Step of analysis (either 1, 2, or 3):
      #### 1: read in data. Outputs: plots to determine filtering criteria.
      #### 2. filter data. Outputs: merged intermediate RData object, elbow plot.
      #### 3. finalize data object. Outputs: finalized RData object for manual analysis.
#### FIFTH ARGUMENT (IF STEP 3): number of dimensions to include.

my_cellranger_csv <- read.csv(as.character(args[2]))
my_filepaths <- as.character(my_cellranger_csv[,1])
my_h5 <- paste0(my_filepaths, ".h5")
my_library_ids <- as.character(my_cellranger_csv[,2])

my_save_path <- as.character(args[3])
my_filter_csv <- read.csv(as.character(args[4]))
#my_filter_csv <- read.csv("~/Desktop/filtering_criteria.csv")
my_step <- as.numeric(args[5])


#### 2. Read in data, store 10X data in a list
if(my_step %in% c(1,2)) {
  my_list <- vector(mode = "list", length = length(my_filepaths))

  for(i in 1:length(my_filepaths)){
    my_initial_data <- Read10X_h5(my_h5[i])
    my_inital_counts <- my_initial_data$`Gene Expression`
    my_initial_object <- CreateSeuratObject(counts = my_inital_counts, project = my_library_ids[i], assay = "RNA")
    my_initial_object[["mito.genes"]] <- PercentageFeatureSet(my_initial_object, pattern = "^MT-")
    my_initial_object[["percent.ribo.large"]] <- PercentageFeatureSet(my_initial_object, pattern = "^RPL")
    my_initial_object[["percent.ribo.small"]] <- PercentageFeatureSet(my_initial_object, pattern = "^RPS")
    my_initial_object[["ribo.genes"]] <- my_initial_object[["percent.ribo.small"]] + my_initial_object[["percent.ribo.large"]]

    p1 <- VlnPlot(object = my_initial_object, features = c("nCount_RNA")) + theme(legend.position = "none")
    p2 <- VlnPlot(object = my_initial_object, features = c("nFeature_RNA")) + theme(legend.position = "none")
    p3 <- VlnPlot(object = my_initial_object, features = c("mito.genes")) + theme(legend.position = "none")

    pdf(paste0(my_save_path, my_library_ids[i], ".pdf"), height = 8, width=12)
    grid.arrange(p1, p2, p3, ncol = 3)
    dev.off()

    my_list[[i]] <- my_initial_object

  }

  input.df <- data.frame(my_filepaths, my_library_ids, rep(my_save_path, length(my_filepaths)))
  write.csv(input.df, paste0(my_save_path, "input_1.csv"))
  }

##### STOP 1


# 3. Print VlnPlots and Specify Cutoffs
if(my_step == 2) {

###############################################
########## PERFORM THE FILTERING ##############
###############################################

 lower.cutoff = my_filter_csv[,1]
  higher.cutoff = my_filter_csv[,2]
  mito.cutoff = my_filter_csv[,3]
  ribo.cutoff = my_filter_csv[,4]
  selection_method = as.character(my_filter_csv[,5])
  variable_genes = my_filter_csv[,6]

  for(i in 1:length(my_filepaths)){
    my_initial_object <- my_list[[i]]
    my_intial_object <- subset(my_initial_object, subset = nFeature_RNA < higher.cutoff & nFeature_RNA > lower.cutoff & mito.genes < mito.cutoff & ribo.genes < ribo.cutoff)
    my_intial_object <- NormalizeData(my_intial_object, normalization.method = "LogNormalize", scale.factor = 10000)
    my_intial_object <- FindVariableFeatures(my_intial_object, selection.method = selection_method, nfeatures = variable_genes)
    my_list[[i]] <- my_intial_object
  }

#########################################
########## PERFORM THE CCA ##############
#########################################
  IVR_object.anchors <- FindIntegrationAnchors(object.list = my_list, dims=1:25)
  IVR_object.combined <- IntegrateData(anchorset = IVR_object.anchors, dims = 1:25)

  DefaultAssay(IVR_object.combined) <- "integrated"
  IVR_object.combined <- ScaleData(IVR_object.combined, verbose = FALSE)
  IVR_object.combined <- RunPCA(IVR_object.combined, npcs = 30, verbose = FALSE)

  ElbowPlot(IVR_object.combined, ndims = 25)
  ggsave(paste0(my_save_path, "elbow_plot.pdf"), height = 8, width = 8)

  save(IVR_object.combined, file = paste0(my_save_path, "intermediate_object.RData"))
}

### STOP 2



# 4. Perform the CCA
if(my_step == c(3)) {

  load(file = paste0(my_save_path, "intermediate_object.RData"))

  #### READ OBJECT
  final.dim = as.numeric(args[6])

  IVR_object.combined <- RunUMAP(IVR_object.combined, reduction = "pca", dims = 1:final.dim)
  IVR_object.combined <- FindNeighbors(IVR_object.combined, reduction = "pca", dims = 1:final.dim)
  IVR_object.combined <- FindClusters(IVR_object.combined, resolution = 0.5)

  save(IVR_object.combined, file = paste0(my_save_path, "final_object.RData"))
  my_labels <- c("libraries", "filtering_criteria", "dimensions")
  my_list <- list(my_filepaths, my_filter_csv, final.dim)
  save(my_list, file = paste0(my_save_path, "documented_criteria_for_final_object.RData"))
}


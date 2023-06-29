##########################
##### READ LIBRARIES #####
##########################

library(crayon, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(backports, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(withr, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(rstudioapi, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(cli, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(Seurat, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(Signac, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(GenomeInfoDb, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(Biobase, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(AnnotationDbi, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(AnnotationFilter, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(GenomicFeatures, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(ensembldb, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(biovizBase, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(EnsDb.Hsapiens.v86, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(hdf5r, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(matrixStats, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(MatrixGenerics, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(SummarizedExperiment, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(labeling, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(tidyverse, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
library(GenomicRanges, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")

library(future, lib.loc = "/Shared/IVR/apvoigt/programs/R-4.1.0/bin/PACKAGES/")
## enable parallel-izaton
options(future.globals.maxSize = 80 * 1024 ^ 3) # for 40 Gb RAM
plan("multiprocess", workers = 10)


#######################################
##### READ COMMAND LINE ARGUMENTS #####
#######################################

args <- commandArgs(trailingOnly = TRUE)

csv_args <- args[2]

#1. Specify filepaths to read in and specify corresponding library IDs.
#### FIRST ARGUMENT: CSV file specifying filepaths of cellranger output and library ID for each file to integrated
#### SECOND ARGUMENT: filepath where data (plots, itermediate objects) can be saved. RECOMMEND CREATING OUTPUT FOLDER WITH DATE.
#### THIRD ARGUMENT: CSV mgatk commands
#### 1. add_mgatk_data: TRUE or FALSE. If TRUE, mgatk data (including mito coverage) will be added as an assay
#### 2. output_path: Full path to directory that contains the mgatk output (different path for each sample).
#### 3. sample_name: sample name (this is used to read in the appropriate RDS file output from mgatk)
#### 4. allele_of_interest: mitochondrial allele of interest to add to meta data (eg MELAS - 3243)
#### 5. allele_1: uppercase character of allele 1 (eg: A)
#### 6. allele_2: uppercase character of allele 2 (eg: G)
#### FOURTH ARGUMENT: csv file of filtering criteria
#### 1. peak_region_fragments_lower
#### 2. peak_region_fragments_upper
#### 3. mtDNA_depth_lower
#### 4. mtDNA_depth_lower
#### FIFTH ARGUMENT: Step of analysis (either 1 or 2):
#### 1: read in data. Outputs: plots to determine filtering criteria.
#### 2. filter data. Outputs: merged intermediate RData object, elbow plot.
#### 3. finalize data object. Outputs: finalized RData object for manual analysis.

my_filepaths <- read.csv(as.character(args[2]))

my_save_path <- as.character(args[3])
mgatk_commands <- read.csv(as.character(args[4]))
filtering_criteria <- read.csv(as.character(args[5]))
my_step <- as.numeric(args[6])

if(my_step == 1){
  ######################
  ##### READ PEAKS #####
  ######################
  print("READ PEAKS")
  #This follows the merging vignette
  #https://github.com/timoast/signac/issues/133#issuecomment-615331021

  ### read in peak sets and convert to genomic ranges
  peaks_list <- vector(mode = "list", length = nrow(my_filepaths))
  for(i in 1:length(peaks_list)){
    print(i)
    my_peaks <- read.table(
      file = paste0(my_filepaths[i,1], "/outs/atac_peaks.bed"), # specify atac_peaks (vs gex)
      col.names = c("chr", "start", "end")
    )

    gr <- makeGRangesFromDataFrame(my_peaks)

    peaks_list[[i]] <- gr
  }

  ###################################
  ##### CREATE UNIFIED PEAK SET #####
  ###################################
  print("CREATE UNIFIED PEAK SET")
  combined.peaks = Signac::reduce(x = do.call(c, unlist(peaks_list, recursive=FALSE)))

  ############################################
  ##### FILTER BAD PEAKS BASED ON LENGTH #####
  ############################################
  print("FILTER PEAKS")
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
  combined.peaks


  ###################################
  ##### CREATE FRAGMENT OBJECTS #####
  ###################################
  print("CREATE FRAGMENT OBJECTS")

  fragments_list <- vector(mode = "list", length = nrow(my_filepaths))
  for(i in 1:length(fragments_list)){
    print(i)
    md <- read.table(
      file = paste0(my_filepaths[i,1], "/outs/per_barcode_metrics.csv"), # try to use the cellranger-arc output
      stringsAsFactors = FALSE,
      sep = ",",
      header = T,
      row.names = 1 # use GEX barcodes - these are in the .h5 matrix!
    )[-1,] # remove the first row

    ## I only want to create a fragment object for the CELLS identified by cellranger
    ## Therefore, I read in the filtered peaks from cellranger

    prelim_counts <- Read10X_h5(filename = paste0(my_filepaths[i,1], "/outs/filtered_feature_bc_matrix.h5"))
    counts <- prelim_counts$Peaks # add to specify atac data from the filtered bc matrix
    md <- md[which(rownames(md) %in% colnames(counts)), ]

    frags <- CreateFragmentObject(
      path = paste0(my_filepaths[i,1], "/outs/atac_fragments.tsv.gz"), # specify atac_fragments
      cells = rownames(md)
    )

    fragments_list[[i]] <- frags

  }


  ##########################
  ##### QUANTIFY PEAKS #####
  ##########################
  print("QUANTIFY PEAKS")

  objects_list <- vector(mode = "list", length = nrow(my_filepaths))
  for(i in 1:length(objects_list)){
    print(i)

    md <- read.table(
      file = paste0(my_filepaths[i,1], "/outs/per_barcode_metrics.csv"), # try to use the cellranger-arc output
      stringsAsFactors = FALSE,
      sep = ",",
      header = T,
      row.names = 1 # use GEX barcodes - these are in the .h5 matrix!
    )[-1,]

    sample.counts <- FeatureMatrix(
      fragments = fragments_list[[i]],
      features = combined.peaks,
      cells = rownames(md)
    )

    object_assay <- CreateChromatinAssay(sample.counts,
                                         fragments = fragments_list[[i]])
    seurat_obj <- CreateSeuratObject(object_assay,
                                     assay = "ATAC",
                                     meta.data = md,
                                     min.features = 200 ## Nate changed 09272022
    )

    ### add dataset information
    seurat_obj$dataset <- as.character(my_filepaths[i,2])

    objects_list[[i]] <- seurat_obj

  }

  save(objects_list, file = paste0(my_save_path, "/object_list.RData"))

  ########################
  ### MGATK GENOTYPING ###
  ########################
  print("MGATK GENOTYPING")

  extractme <- function(input_SE, position, letter){
    extracted <- assays(SE)[[paste0(letter, "_counts_fw")]][position, ] + assays(SE)[[paste0(letter, "_counts_rev")]][position, ]
    return(extracted)
  }
  e2 <- function(SE, position, letter1, letter2){
    df <- data.frame(
      m1 = extractme(SE, position, letter1),
      m2 = extractme(SE, position, letter2)
    )
    df$het <- round(df$m2 / (df$m2 + df$m1 + 0.0001), 3)
    df$cov <- df$m1 + df$m2
    colnames(df) <- c(paste0("m", as.character(position), "_", letter1), paste0("m", as.character(position), "_", letter2),
                      paste0("het", as.character(position), letter1, "_", letter2),
                      paste0("cov", as.character(position)))
    return(df)
  }

  if(mgatk_commands[["add_mgatk_data"]][1] == TRUE){

    for(i in 1:nrow(mgatk_commands)){
      # load mgatk output
      # this output is the mgatk.A/T/C/G.txt.gz files etc...
      mito.data <- ReadMGATK(dir = mgatk_commands[["output_path"]][i])

      # create an assay
      mito <- CreateAssayObject(counts = mito.data$counts)
      # Subset to cell present in the scATAC-seq assat
      mito <- subset(mito, cells = colnames(objects_list[[i]]))
      # add assay and metadata to the seurat object
      objects_list[[i]][["mito"]] <- mito
      objects_list[[i]] <- AddMetaData(objects_list[[i]],
                                       metadata = mito.data$depth,
                                       col.name = "mtDNA_depth")

      SE <- readRDS(paste0(mgatk_commands[["output_path"]][i], mgatk_commands[["sample_name"]][i], ".rds"))

      allele_coverage_per_cell <- e2(SE = SE, position = 3243, letter1 = "A", letter2 = "G") %>% rownames_to_column("barcode")

      ## ensure cells are in the same order as the metadata
      allele_coverage_md <- objects_list[[i]]@meta.data %>%
        rownames_to_column("barcode") %>%
        left_join(allele_coverage_per_cell)
      allele_coverage_md <- allele_coverage_md[colnames(allele_coverage_per_cell)]

      ## add into the meta data
      for(column in 2:ncol(allele_coverage_md)){ ## start at 2, skip barcode
        objects_list[[i]][[colnames(allele_coverage_md)[column]]] <- allele_coverage_md[[column]]
      }

    # ## create allele assay
    # print("CREATING ALLELE ASSAY")
    # variable.sites <- IdentifyVariants(objects_list[[i]], assay = "mito", refallele = mito.data$refallele)
    # high.conf <- subset(
    #   variable.sites, subset = n_cells_conf_detected >= 5 &
    #     strand_correlation >= 0.65 #&
    #     #vmr > 0.01
    # )

    # objects_list[[i]] <- AlleleFreq(
    #   object = objects_list[[i]],
    #   variants = high.conf$variant,
    #   assay = "mito"
    #  )

    }

  }

  ################################
  ##### CALCULATE QC METRICS #####
  ################################
  print("CALCULATE QC")

  for(i in 1:length(objects_list)){

    seurat_obj <- objects_list[[i]]

    # Augment QC metrics that were computed by cellranger-atac
    seurat_obj$pct_reads_in_peaks <- seurat_obj$atac_peak_region_fragments / seurat_obj$atac_fragments * 100

    # visualize QC metrics for each cell
    seurat_obj$pct_reads_in_peaks <- seurat_obj$atac_peak_region_fragments / seurat_obj$atac_fragments * 100

    # NKM added this
    seurat_obj$mitochondrial_ratio <- seurat_obj$atac_mitochondrial_reads/seurat_obj$atac_raw_reads * 100

    if(mgatk_commands[["add_mgatk_data"]][1] == TRUE){
       VlnPlot(
         object = seurat_obj,
         features = c('pct_reads_in_peaks',
                      'atac_peak_region_fragments',
                      'mitochondrial_ratio',
                      'mtDNA_depth'),
         ncol = 2
       ) + ggsave(paste0(my_save_path, as.character(my_filepaths[i,2]), ".pdf"), height = 7, width = 7)
    } else {
      VlnPlot(
         object = seurat_obj,
         features = c('pct_reads_in_peaks',
                      'atac_peak_region_fragments',
                      'mitochondrial_ratio'),
         ncol = 2
       ) + ggsave(paste0(my_save_path, as.character(my_filepaths[i,2]), ".pdf"), height = 7, width = 7)
   }
    objects_list[[i]] <- seurat_obj

  }
  save(objects_list, file = paste0(my_save_path, "/object_list_prefilter.RData"))


}

##################
##### FILTER #####
##################
if(my_step == 2){
  load(paste0(my_save_path, "/object_list_prefilter.RData"))

  print("FILTERING")

  peak_region_fragments_lower <- filtering_criteria[["peak_region_fragments_lower"]]
  peak_region_fragments_upper <- filtering_criteria[["peak_region_fragments_upper"]]
  mtDNA_depth_lower <- filtering_criteria[["mtDNA_depth_lower"]]
  mtDNA_depth_upper <- filtering_criteria[["mtDNA_depth_upper"]]

  for(i in 1:length(objects_list)){

    seurat_obj <- objects_list[[i]]

    if(mgatk_commands[["add_mgatk_data"]][1] == TRUE){

       seurat_obj_filtered <- subset(
         x = seurat_obj,
         subset = atac_peak_region_fragments > peak_region_fragments_lower &
           atac_peak_region_fragments < peak_region_fragments_upper &
           mtDNA_depth > mtDNA_depth_lower &
           mtDNA_depth < mtDNA_depth_upper
       )
    } else {
       seurat_obj_filtered <- subset(
         x = seurat_obj,
         subset = atac_peak_region_fragments > peak_region_fragments_lower &
           atac_peak_region_fragments < peak_region_fragments_upper
       )
   }

    objects_list[[i]] <- seurat_obj_filtered

  }

  #######################
  ##### COMPUTE LSI #####
  #######################
  print("COMPUTING LSI ON INDIVIDUAL OBJECTS")

  for(i in 1:length(objects_list)){
    print(i)
    seurat_obj <- objects_list[[i]]
    seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 10)
    seurat_obj <- RunTFIDF(seurat_obj)
    seurat_obj <- RunSVD(seurat_obj)
    objects_list[[i]] <- seurat_obj
  }

  save(objects_list, file = paste0(my_save_path, "/object_list_normalized.RData"))

  #################
  ##### MERGE #####
  #################
  print("MERGING")

  combined <- merge(
    x = objects_list[[1]],
    y = objects_list[c(2:length(objects_list))]
  )
  save(combined, file = paste0(my_save_path, "/merged.RData"))

  ##############################
  ##### COMPUTE MERGED LSI #####
  ##############################
  print("COMPUTING MERGED LSI")

  combined <- RunTFIDF(combined)
  combined <- FindTopFeatures(combined, min.cutoff = 10)
  combined <- RunSVD(combined)
  combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

  save(combined, file = paste0(my_save_path, "merged_final.RData"))

  #####################
  ##### INTEGRATE #####
  #####################
  print("FINDING INTEGRATION ANCHORS")

  integration.anchors <- FindIntegrationAnchors(
    object.list = objects_list,
    anchor.features = rownames(combined),
    reduction = "rlsi",
    dims = 2:30,
  )

  print(" INTEGRATING EMBEDDINGS")

  # integrate LSI embeddings
  integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = combined[["lsi"]],
    new.reduction.name = "integrated_lsi",
    k.weight = 50,
    dims.to.integrate = 1:30 ## TASK: make dims dynamic
  )

  # create a new UMAP using the integrated embeddings
  integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30) ## TASK: make dims dynamic

  ################
  ##### SAVE #####
  ################
  save(integrated, file = paste0(my_save_path, "/integrated_final.RData"))

  # lets save some things from the run data
  #current_date <- Sys.Date()
  seurat_version <- packageVersion("Seurat")
  signac_version <- packageVersion("Signac")

  save_list <- list(
    #date = current_date,
    seurat_version = seurat_version,
    signac_version = signac_version,
    input_files = my_filepaths,
    mgatk_commands = mgatk_commands,
    filtering_criteria = filtering_criteria
  )

  save(save_list, file = paste0(my_save_path, "/", "run-commands.RData"))


}


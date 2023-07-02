library(Seurat)
library(phateR)
library(tidyverse)
library(future)
load(".../integrated_data/nr2e3_nrl_all_timepoints_final_object.RData")

Idents(IVR_object.combined) <- "predicted.id"
DimPlot(IVR_object.combined)


#######################
### CLEAN META DATA ###
#######################

head(IVR_object.combined@meta.data)
md <- IVR_object.combined@meta.data[,c(1:10, 30, 31, 44)]
md <- md %>%
  mutate(line = ifelse(is.na(line) == FALSE, line,
                       ifelse(str_detect(orig.ident, "L75P"), "NRL_L75P", "NRL_WT")))

md <- md %>%
  mutate(timepoint = ifelse(is.na(timepoint) == FALSE, timepoint,
                       ifelse(str_detect(orig.ident, "d100"), "d100", 
                              ifelse(str_detect(orig.ident, "d101"), "d101", 
                                     ifelse(str_detect(orig.ident, "d103"), "d103", 
                                            ifelse(str_detect(orig.ident, "D260"), "d260", "d170"))))))


###############################
### CLASSIFY DIVERGENT RODS ###
###############################
table(md$predicted.id, md$timepoint)

md <- md %>%
  mutate(celltype = ifelse(line == "B342" & predicted.id == "Rod" & timepoint %in% c("D160", "D260"), "div_rod", predicted.id))


#####################
### CLASSIFY CODS ###
#####################
IVR_object.combined@meta.data <- md
Idents(IVR_object.combined) <- "line"
nrl <- subset(IVR_object.combined, idents = c("NRL_L75P", "NRL_WT"))
DefaultAssay(nrl) <- "RNA"

### cods have LOW expression of rod genes: ROM1, PDE6G, SAG, NRL NR2E3, CNGB1, GNAT1
### cods have LOW expression of cone genes: ARR3, CNGB3, GNAT2, GNB3, GNGT2, GUCA1A, PDE6C, PDE6H
### cods have HIGH expression of OPN1SW
rod_module <- list(c("ROM1", "PDE6G", "SAG", "NRL", "NR2E3", "CNGB1", "GNAT1"))
cone_module <- list(c("ARR3", "CNGB3", "GNAT2", "GNB3", "GNGT2", "GUCA1A", "PDE6C", "PDE6H"))

nrl <- AddModuleScore(
  object = nrl,
  features = rod_module,
  name = 'rod_features'
)

nrl <- AddModuleScore(
  object = nrl,
  features = cone_module,
  name = 'cone_features'
)

FeaturePlot(nrl, features = c("rod_features1", "cone_features1", "OPN1SW"), ncol = 3)

OPN1SW <- nrl@assays$RNA@data[which(rownames(nrl@assays$RNA@data) == "OPN1SW"),]
nrl[["OPN1SW"]] <- OPN1SW

cod_md <- nrl@meta.data %>%
  rownames_to_column("barcode") %>%
  as_tibble() %>%
  filter(rod_features1 < 2 &
         cone_features1 < 2 &
         OPN1SW > 0.5) %>%
  select(-c(nCount_RNA, nFeature_RNA, mito.genes, ribo.genes, 
    seurat_clusters, donor_annotation, percent.ribo.large, percent.ribo.small))

summary(as.factor(cod_md$line)) ## 95% of these cells belong to the NRL mutant line, 
                                ## which is very reassuring

cod_barcodes <- cod_md %>%
  filter(line == "NRL_L75P") %>%
  filter(predicted.id %in% c("Rod", "Cone", "T1", "T3", "PRC_Progenitors")) %>%
  pull(barcode)

md <- md %>% 
  rownames_to_column("barcode") %>%
  mutate(celltype = ifelse(barcode %in% cod_barcodes, "cod", celltype)) %>%
  column_to_rownames("barcode")

IVR_object.combined@meta.data <- md

###################################################
### VISUALIZE CODS AND DIVERGENT RODS WITH UMAP ###
###################################################
Idents(IVR_object.combined) <- "celltype"
DimPlot(IVR_object.combined, label = T) ## not super helpful


pr <- subset(IVR_object.combined, idents = c("Rod", "Cone", "Progenitors", "PRC_Progenitors", "T1", "T3", "cod", "div_rod"))
DimPlot(pr, label = T, cols = c(rep("grey", 5), "red", "grey", "blue")) ## <1000 'cods'


##################################################
### DIFF EXP BETWEEN COD AND DIV ROD ###
##################################################

DefaultAssay(IVR_object.combined) <- "RNA"
Idents(IVR_object.combined) <- "celltype"
dge <- FindMarkers(IVR_object.combined, ident.1 = "div_rod", ident.2 = "cod", min.pct = 0, logfc.threshold = 0)

write.csv(dge, ".../integrated_data/2023_01_22_div_rod_vs_cod_dge.csv")


##################################################
### DIFF EXP BETWEEN D170 COD AND D160 DIV ROD ###
##################################################
Idents(IVR_object.combined) <- "timepoint"
d160170 <- subset(IVR_object.combined, idents = c("D160", "d170"))


DefaultAssay(d160170) <- "RNA"
Idents(d160170) <- "celltype"
plan("multiprocess", workers = 4)
dge <- FindMarkers(d160170, ident.1 = "div_rod", ident.2 = "cod", min.pct = 0, logfc.threshold = 0)

write.csv(dge, ".../integrated_data/2023_01_28_d160_div_rod_vs_d160_cod_dge.csv")

####################################################
### DIFF EXP BETWEEN D160 ROD AND D170 ROD ###
####################################################
celltype_timepoint <- d160170@meta.data %>%
  as_tibble() %>%
  select(timepoint, celltype) %>%
  rowwise() %>%
  mutate(celltype_timepoint = paste0(celltype, "-",timepoint)) %>%
  pull(celltype_timepoint)

d160170[["celltype_timepoint"]] <- celltype_timepoint

DefaultAssay(d160170) <- "RNA"
Idents(d160170) <- "celltype_timepoint"
plan("multiprocess", workers = 6)
dge <- FindMarkers(d160170, ident.1 = "Rod-D160", ident.2 = "Rod-d170", min.pct = 0, logfc.threshold = 0)

write.csv(dge, ".../integrated_data/2023_01_28_d160_rod_vs_d170_rod_dge.csv")

##########################################################
### DIFF EXP BETWEEN D170 COD and D170 NRL CONTROL ROD ###
##########################################################

DefaultAssay(d160170) <- "RNA"
Idents(d160170) <- "celltype_timepoint"
plan("multiprocess", workers = 6)
dge <- FindMarkers(d160170, ident.1 = "cod-d170", ident.2 = "Rod-d170", min.pct = 0, logfc.threshold = 0)

write.csv(dge, ".../integrated_data/2023_01_28_d170_cod_vs_d170_rod_dge.csv")

#################################################################
### DIFF EXP BETWEEN DIVERGENT ROD AND D160 CONTROL NR2E3 ROD ###
#################################################################

DefaultAssay(d160170) <- "RNA"
Idents(d160170) <- "celltype_timepoint"
plan("multiprocess", workers = 6)
dge <- FindMarkers(d160170, ident.1 = "div_rod-D160", ident.2 = "Rod-D160", min.pct = 0, logfc.threshold = 0)

write.csv(dge, "~/LSS/IVR/apvoigt/experiments/external_sc_rna_seq/Kallman_NRL/integrated_data/2023_01_28_d160_div_rod_vs_d160_rod_dge.csv")

####################################################
### INTEGRATE DIFFERENTIAL EXPRESSION AND FILTER ###
####################################################

dge1 <- read_csv(".../integrated_data/2023_01_22_div_rod_vs_cod_dge.csv")
dge2 <- read_csv(".../integrated_data/2023_01_28_d160_div_rod_vs_d160_cod_dge.csv")
dge3 <- read_csv(".../integrated_data/2023_01_28_d160_rod_vs_d170_rod_dge.csv")
dge4 <- read_csv(".../integrated_data/2023_01_28_d170_cod_vs_d170_rod_dge.csv")
dge5 <- read_csv(".../integrated_data/2023_01_28_d160_div_rod_vs_d160_rod_dge.csv")

colnames(dge1)[1] <- "gene"
colnames(dge2)[1] <- "gene"
colnames(dge3)[1] <- "gene"
colnames(dge4)[1] <- "gene"
colnames(dge5)[1] <- "gene"

dge1 <- dge1 %>% select("gene", "avg_log2FC")
dge2 <- dge2 %>% select("gene", "avg_log2FC")
dge3 <- dge3 %>% select("gene", "avg_log2FC")
dge4 <- dge4 %>% select("gene", "avg_log2FC")
dge5 <- dge5 %>% select("gene", "avg_log2FC")

colnames(dge1) <- c("gene", "DivRod_Cod_FC") 
colnames(dge2) <- c("gene", "D160DivRod_D170Cod_FC")
colnames(dge3) <- c("gene", "D160Rod_D170Rod_FC") 
colnames(dge4) <- c("gene", "D170Cod_D170Rod_FC") 
colnames(dge5) <- c("gene", "D160DivRod_D160Rod_FC") 

dge_merged <- dge1 %>%
  left_join(dge2) %>%
  left_join(dge3) %>%
  left_join(dge4) %>%
  left_join(dge5)

write.csv(dge_merged, ".../integrated_data/2023_01_28_all_dge_merged.csv")

## We want to find genes that are...
## A: DIFFERENT between cod and divergent rod 
## B: SIMILAR between control rods 
## C: DIFFERENT between cods and (NRL control) rods 
## D: DIFFERENT between divergent rods and (NR2E3 control) rods 

dge_merged %>%
  filter(abs(D160DivRod_D170Cod_FC) > 0.75) %>% ## A
  filter(abs(D160Rod_D170Rod_FC) < 0.25) %>%    ## B
  filter(abs(D170Cod_D170Rod_FC) > 0.25) %>%    ## C
  filter(abs(D160DivRod_D160Rod_FC) > 0.25)     ## D

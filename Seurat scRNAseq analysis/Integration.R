
# Setup all libraries needed
library(Seurat)
library(SeuratDisk) #!!!!!!!
library(SeuratWrappers) #!!!!!!

library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(knitr)

library(reticulate) ## allows python functions to be called in R
library(Matrix) ## I need to read in a .mtx file

library(AnnotationDbi) ## to convert ensembl ID to symbol
library(org.Mm.eg.db)
#if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#BiocManager::install("org.Mm.eg.db")
library(viridis) ## beautiful and scientifically better heatmap colour schemes

setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1WuaGy1fdZTEwlHzeOEMCyIpn9pRIgui_/TFM Inigo Clemente")

################################################

# Read both datasets from scRNAseq alignment
p28.sham<-Read10X(
  data.dir = "Cellranger_HPC/P28_GFP/filtered_feature_bc_matrix/",
  gene.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
) 
p28.sham

p31.sham<-Read10X(
  data.dir = "Cellranger_HPC/P31_GFP/filtered_feature_bc_matrix/",
  gene.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
) 
p31.sham

################################################

# Change genes names from ENSEMBL ID to gene symbol
x<- select(x = org.Mm.eg.db, 
           keys = rownames(p28.sham), 
           column = "SYMBOL", 
           keytype = "ENSEMBL",
           multiVals = "first",
           asNA=F)
head(x)
x$FINAL <- ifelse(is.na(x$SYMBOL), x$ENSEMBL, x$SYMBOL) # issues with excess matches so collapse into a single col
for (name in rownames(p28.sham)) {
  rownames(p28.sham)[match(name, rownames(p28.sham))]<-x$FINAL[match(name, x$ENSEMBL)]
}
head(rownames(p28.sham))

x<- select(x = org.Mm.eg.db, 
           keys = rownames(p31.sham), 
           column = "SYMBOL", 
           keytype = "ENSEMBL",
           multiVals = "first",
           asNA=F)
head(x)
x$FINAL <- ifelse(is.na(x$SYMBOL), x$ENSEMBL, x$SYMBOL) # issues with excess matches so collapse into a single col
for (name in rownames(p31.sham)) {
  rownames(p31.sham)[match(name, rownames(p31.sham))]<-x$FINAL[match(name, x$ENSEMBL)]
}
head(rownames(p31.sham))

################################################

# Generate Seurat objects for both datasets
p28 <- CreateSeuratObject(
  p28.sham,
  project = "p28_filtered",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)
p28

p31 <- CreateSeuratObject(
  p31.sham,
  project = "p31_filtered",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)
p31

################################################

# Quality control for both datasets
p28[["percent.mt"]] <- PercentageFeatureSet(p28, pattern = "^Mt") ## if human, "MT"
p28[["percent.rbp"]] <- PercentageFeatureSet(p28, pattern = "^Rp[sl]") 
#VlnPlot(p28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4)
p28 <- subset(p28, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
# Filtering of nFeature_RNA allows get rid of empty droplets or duplets
# Filtering of percent.mt allows get rid of cells in bad cell cycle state (dying)

p31[["percent.mt"]] <- PercentageFeatureSet(p31, pattern = "^Mt") ## if human, "MT"
p31[["percent.rbp"]] <- PercentageFeatureSet(p31, pattern = "^Rp[sl]") 
#VlnPlot(p31, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4)
p31 <- subset(p31, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)

################################################

table(rownames(p28) %in% rownames(p31)) 

################################################

# Generate a list of the datasets to integrate
pbmc_list <- list()
pbmc_list[["p28"]] <- p28
pbmc_list[["p31"]] <- p31

set.seed(123)
for (i in 1:length(pbmc_list)) {
  pbmc_list[[i]] <- NormalizeData(pbmc_list[[i]], verbose = F)
  pbmc_list[[i]] <- FindVariableFeatures(pbmc_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}

################################################

# Setting the anchors for the integration
set.seed(123)
pbmc_anchors    <- FindIntegrationAnchors(object.list = pbmc_list, dims = 1:30)

################################################

# Integration of data based on the anchors
set.seed(123)
pbmc_seurat     <- IntegrateData(anchorset = pbmc_anchors, dims = 1:30)

################################################

# Representation of both datasets separated
DefaultAssay(pbmc_seurat) <- "RNA"

set.seed(123)
pbmc_seurat <- NormalizeData(pbmc_seurat, verbose = F)
set.seed(123)
pbmc_seurat <- FindVariableFeatures(pbmc_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
set.seed(123)
pbmc_seurat <- ScaleData(pbmc_seurat, verbose = F)
set.seed(123)
pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = F)
set.seed(123)
pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:30, verbose = F)

set.seed(123)
DimPlot(pbmc_seurat,reduction = "umap") + plot_annotation(title = "p28 vs p31 , before integration")

################################################

# Representation of both datasets integrated
DefaultAssay(pbmc_seurat) <- "integrated"

set.seed(123)
pbmc_seurat <- ScaleData(pbmc_seurat, verbose = F)
set.seed(123)
pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = F)
set.seed(123)
pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:30, verbose = F)

set.seed(123)
DimPlot(pbmc_seurat, reduction = "umap") + plot_annotation(title = "p28 vs p31, after integration (Seurat 3)")

set.seed(123)
DimPlot(pbmc_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend()
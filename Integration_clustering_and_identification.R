
################################################

# 0. Libraries setup

################################################

library(Seurat)
library(SeuratDisk) 
library(SeuratWrappers) 

library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(knitr)

library(reticulate)
library(Matrix)

library(AnnotationDbi)
library(org.Mm.eg.db)
#if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#BiocManager::install("org.Mm.eg.db")
library(viridis)

setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1WuaGy1fdZTEwlHzeOEMCyIpn9pRIgui_/TFM Inigo Clemente")

################################################

# 1. scRNAseq integration

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


# Change genes names from ENSEMBL ID to gene symbol
x<- select(x = org.Mm.eg.db, 
           keys = rownames(p28.sham), 
           column = "SYMBOL", 
           keytype = "ENSEMBL",
           multiVals = "first",
           asNA=F)
head(x)
x$FINAL <- ifelse(is.na(x$SYMBOL), x$ENSEMBL, x$SYMBOL)
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
x$FINAL <- ifelse(is.na(x$SYMBOL), x$ENSEMBL, x$SYMBOL)
for (name in rownames(p31.sham)) {
  rownames(p31.sham)[match(name, rownames(p31.sham))]<-x$FINAL[match(name, x$ENSEMBL)]
}
head(rownames(p31.sham))


# Generate Seurat objects for both datasets
p28 <- CreateSeuratObject(
  p28.sham,
  project = "pBIC",
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
  project = "pBIC10",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)
p31


# Quality control for both datasets
p28[["percent.mt"]] <- PercentageFeatureSet(p28, pattern = "^Mt")
p28[["percent.rbp"]] <- PercentageFeatureSet(p28, pattern = "^Rp[sl]") 
VlnPlot(p28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4)
p28 <- subset(p28, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
# Filtering of nFeature_RNA allows get rid of empty droplets or duplets
# Filtering of percent.mt allows get rid of cells in bad cell cycle state (dying)
VlnPlot(p28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4)

p31[["percent.mt"]] <- PercentageFeatureSet(p31, pattern = "^Mt")
p31[["percent.rbp"]] <- PercentageFeatureSet(p31, pattern = "^Rp[sl]") 
VlnPlot(p31, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4)
p31 <- subset(p31, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
VlnPlot(p31, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4)


table(rownames(p28) %in% rownames(p31)) 


# Generate a list of the datasets to integrate
pBIC_list <- list()
pBIC_list[["p28"]] <- p28
pBIC_list[["p31"]] <- p31

set.seed(123)
for (i in 1:length(pBIC_list)) {
  pBIC_list[[i]] <- NormalizeData(pBIC_list[[i]], verbose = F) # Normalized by log normalizarion
  pBIC_list[[i]] <- FindVariableFeatures(pBIC_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F) # Expected technical variation is corrected based on High Variable Genes (HVG)
}


# Setting the anchors for the integration
set.seed(123)
pBIC_anchors    <- FindIntegrationAnchors(object.list = pBIC_list, dims = 1:30) # Dims are stated up to 30 as the maximun, maximazing the differences between features and aiming to identify the best possible cell representation in both Seurat objects


# Integration of data based on the anchors
set.seed(123)
pBIC_seurat     <- IntegrateData(anchorset = pBIC_anchors, dims = 1:30) # Dataset integration is performed using the set of anchors found in the step before as a pre-computed AnchorSet


# Representation of both datasets separated
DefaultAssay(pBIC_seurat) <- "RNA"

set.seed(123)
pBIC_seurat <- NormalizeData(pBIC_seurat, verbose = F)
set.seed(123)
pBIC_seurat <- FindVariableFeatures(pBIC_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
set.seed(123)
pBIC_seurat <- ScaleData(pBIC_seurat, verbose = F)
set.seed(123)
pBIC_seurat <- RunPCA(pBIC_seurat, npcs = 30, verbose = F)
set.seed(123)
pBIC_seurat <- RunUMAP(pBIC_seurat, reduction = "pca", dims = 1:30, verbose = F)

DimPlot(pBIC_seurat,reduction = "umap") + 
  plot_annotation(title = "pBIC vs pBIC10 , before integration")

# Representation of both datasets integrated
DefaultAssay(pBIC_seurat) <- "integrated"

set.seed(123)
pBIC_seurat <- ScaleData(pBIC_seurat, verbose = F)
set.seed(123)
pBIC_seurat <- RunPCA(pBIC_seurat, npcs = 30, verbose = F)
set.seed(123)
pBIC_seurat <- RunUMAP(pBIC_seurat, reduction = "pca", dims = 1:30, verbose = F)

DimPlot(pBIC_seurat, reduction = "umap") + 
  plot_annotation(title = "pBIC vs pBIC10, after integration (Seurat 3)")
DimPlot(pBIC_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend()

################################################

# 2. Clustering

################################################

# Cluserting
set.seed(123)
pBIC_seurat <- FindNeighbors(pBIC_seurat, dims = 1:30, verbose = F) # K-nearest neighbors and then constructs the SNN graph by the Jaccard index for neigbourhood
set.seed(123)
pBIC_seurat <- FindClusters(pBIC_seurat, resolution = 0.5, verbose = F) # Clusters of cells are identified
# KEY POINT: Resolution of 0.5 is set to obtain a suitable number of clusters

DimPlot(pBIC_seurat,label = T)
DimPlot(pBIC_seurat, reduction = "umap", split.by = "orig.ident", label = T) + NoLegend()

# Clusters cell counts
count_table <- table(pBIC_seurat@meta.data$seurat_clusters, 
                     pBIC_seurat@meta.data$orig.ident)
count_table
populations <- data.frame(count_table)
populations_o <- populations %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  mutate(pBIC_pop = pBIC/sum(pBIC)*100, 
         pBIC10_pop = pBIC10/sum(pBIC10)*100)
names(populations_o)[names(populations_o) == 'Var1'] <- 'Cluster'
populations_o
write.csv(populations_o,"Cell_populations.csv") # Saves the dataframe in a csv for later analysis

################################################

# 3. Cluster identity

################################################

# Clusters identification
DefaultAssay(pBIC_seurat) <- "RNA"

#png(filename = "dot B cells feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Cd19","Ms4a1"),dot.scale = 15)
#dev.off()

#png(filename = "dot B cells  tumoral-non-tumoral feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Cd19","Ms4a1","GFP"),dot.scale = 15)
FeaturePlot(object = pBIC_seurat,features = c("Cd19","Ms4a1","GFP"),order=T, pt.size = 0.1)
#dev.off()

#png(filename = "dot T cells feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Cd3e"),dot.scale = 15)
#dev.off()

#png(filename = "dot T cells CD4-CD8 feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Cd3e","Cd4","Cd8a"),dot.scale = 15)
FeaturePlot(object = pBIC_seurat,features = c("Cd3e","Cd4","Cd8a"),order=T, pt.size = 0.1)
#dev.off()

#png(filename = "dot NK cells feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Klrk1","Ncr1","Klrd1","Prf1"),dot.scale = 15)
#dev.off()

#png(filename = "dot NK cells NK-NKT feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Klrk1","Ncr1","Klrd1","Prf1","Cd3e"),dot.scale = 15)
FeaturePlot(object = pBIC_seurat,features = c("Klrk1","Ncr1","Klrd1","Prf1","Cd3e"),order=T, pt.size = 0.1)
#dev.off()

#png(filename = "dot DC cells feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Itgam"),dot.scale = 15)
#dev.off()

#png(filename = "dot DC cells DC-GrMDSC-MoMDSC feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Itgam","Itgax","H2-Ab1","Ly6g","Ly6c1"),dot.scale = 15)
FeaturePlot(object = pBIC_seurat,features = c("Itgam","Itgax","H2-Ab1","Ly6g","Ly6c1"),order=T, pt.size = 0.1)
#dev.off()

#png(filename = "dot Mph cells feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = pBIC_seurat,features = c("Adgre1","Ly6c1","Ly6g"),dot.scale = 15)
FeaturePlot(object = pBIC_seurat,features = c("Adgre1","Ly6c1","Ly6g"),order=T, pt.size = 0.1)
#dev.off()

pBIC_seurat_copy <- pBIC_seurat
pBIC_seurat <- RenameIdents(object = pBIC_seurat, 
                            `0` = "T cells",`1` = "T cells",`2` = "DC",
                            `3` = "NKT",`4` = "B cells non-tumoral",
                            `5` = "GrMDSC",`6` = "B cells tumoral",
                            `7` = "T cells", `8` = "B cells tumoral",
                            `9` = "T cells",`10` = "B cells tumoral",
                            `11` = "Not classified",`12` = "T cells", 
                            `13` = "Not classified",`14` = "Mph", `15` = "DC",
                            `16` = "GrMDSC", `17` = "NK",`18` = "Mph", 
                            `19` = "Not classified", `20` = "MoMDSC") # Set the cell population identification to each cluster
DimPlot(pBIC_seurat, reduction = "umap", split.by = "orig.ident", label = T)

################################################

# 4. T cells subclustering and identification

################################################

DefaultAssay(pBIC_seurat_copy) <- "RNA"

T_Cells <- subset(pBIC_seurat_copy, idents= c("0","1","7","9","12"))
T_Cells.list <- SplitObject(T_Cells, split.by = "orig.ident")

set.seed(123)
for (i in 1:length(T_Cells.list)) {
  T_Cells.list[[i]] <- NormalizeData(T_Cells.list[[i]], verbose = F)
  T_Cells.list[[i]] <- FindVariableFeatures(T_Cells.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}

reference.list <- T_Cells.list
Integrated_Tcells <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
Integrated_Tcells <- IntegrateData(anchorset = Integrated_Tcells, dims = 1:30)
Integrated_Tcells <- ScaleData(Integrated_Tcells, verbose = FALSE)
Integrated_Tcells <- RunPCA(Integrated_Tcells, npcs = 30, verbose = FALSE)
Integrated_Tcells <- RunUMAP(Integrated_Tcells, reduction = "pca", dims = 1:30, verbose = F)
Integrated_Tcells <- FindNeighbors(Integrated_Tcells, dims = 1:30, verbose = F)
Integrated_Tcells <- FindClusters(Integrated_Tcells, resolution = 0.4, verbose =F)
DimPlot(Integrated_Tcells,label = T)
DimPlot(Integrated_Tcells, reduction = "umap", split.by = "orig.ident", label = T) + NoLegend()

DefaultAssay(Integrated) <- "RNA"
#png(filename = "Tcells feature T cells CD4-CD8 feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = Integrated_Tcells,features = c("Cd3e","Cd4","Cd8a"),dot.scale = 15)
FeaturePlot(object = Integrated_Tcells,features = c("Cd3e","Cd4","Cd8a"),order=T, pt.size = 0.1)
#dev.off()

#png(filename = "Tcells feature CD8 cells naive-central-memory feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = Integrated_Tcells,features = c("Cd3e","Cd8a","Cd44","Sell","Pdcd1"),dot.scale = 15)
FeaturePlot(object = Integrated_Tcells,features = c("Cd3e","Cd8a","Cd44","Sell","Pdcd1"),order=T, pt.size = 0.1)
#dev.off()

#png(filename = "Tcells feature CD4 cells cells naive-central-memory-Tfh-Treg feature.png", width = 1000, height = 1000, pointsize = 12,family = "sans", bg = "white")
DotPlot(object = Integrated_Tcells,features = c("Cd3e","Cd4","Cd44","Sell","Pdcd1","Cxcr5","Bcl6","Il2ra","Foxp3"),dot.scale = 15)
FeaturePlot(object = Integrated_Tcells,features = c("Cd3e","Cd4","Cd44","Sell","Pdcd1","Cxcr5","Bcl6","Il2ra","Foxp3"),order=T, pt.size = 0.1)
#dev.off()

Integrated_Tcells <- RenameIdents(object = Integrated_Tcells, 
                           `0` = "CD8 pre-effector memory",
                           `1` = "CD8 central memory",
                           `2` = "CD8 pre-effector memory",
                           `3` = "CD4 effector memory",
                           `4` = "CD8 effector memory",
                           `5` = "CD8 central memory",
                           `6` = "CD8 central memory",
                           `7` = "CD8 central memory", 
                           `8` = "CD8 central memory",
                           `9` = "CD8 effector memory",
                           `10` = "CD4 naive", 
                           `11` = "CD8 central memory",
                           `12` = "CD4 Treg")
DimPlot(Integrated_Tcells,label = T)
DimPlot(Integrated_Tcells, reduction = "umap", split.by = "orig.ident", label = T) + NoLegend()

count_table <- table(Integrated_Tcells@meta.data$seurat_clusters, Integrated_Tcells@meta.data$orig.ident)
count_table

populations <- data.frame(count_table)
populations_T <- populations %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  mutate(pBIC_pop = p28_filtered/sum(p28_filtered)*100, pBIC10_pop = p31_filtered/sum(p31_filtered)*100)
populations_T
write.csv(populations_T,"Tcell_populations.csv")

################################################

# 5. Automatic cluster identity with SingleR

################################################

library(SingleR)
library(celldex)
library(SingleCellExperiment)

mouse.ref <- celldex::MouseRNAseqData() # Out of the variety of reference genome, we keep MpuseRNAseqData

sce <- as.SingleCellExperiment(DietSeurat(pBIC_seurat))
sce

mouse.main <- SingleR(test = sce,assay.type.test = 1,ref = mouse.ref,labels = mouse.ref$label.main)
mouse.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouse.ref,labels = mouse.ref$label.fine)

table(mouse.fine$pruned.labels)

pBIC_seurat@meta.data$mouse.main   <- mouse.main$pruned.labels
pBIC_seurat@meta.data$mouse.fine   <- mouse.fine$pruned.labels

pBIC_seurat <- SetIdent(pBIC_seurat, value = "mouse.fine")

DimPlot(pBIC_seurat, reduction = "umap")
DimPlot(pBIC_seurat, reduction = "umap", split.by = "orig.ident")

################################################

# 6. Cell populations table

################################################

setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1WuaGy1fdZTEwlHzeOEMCyIpn9pRIgui_/TFM Inigo Clemente/Cell_populations")

populations <- read.csv("Cell_populations.csv")
populations_T <- read.csv("Tcell_populations.csv")

populations_r <- populations
populations_r[1,] <- populations[1,]+populations[2,]+populations[8,]+populations[10,]+populations[13,]
populations_r[2,] <- populations[3,]+populations[16,]
populations_r[3,] <- populations[5,]
populations_r[4,] <- populations[6,]+populations[17,]
populations_r[5,] <- populations[7,]+populations[9,]+populations[11,]
populations_r[6,] <- populations[12,]+populations[14,]+populations[20,]
populations_r[7,] <- populations[15,]+populations[19,]
populations_r[8,] <- populations[18,]
populations_r[9,] <- populations[4,]
populations_r[10,] <- populations[21,]
populations_r <- populations_r[-c(11:21), ] 
populations_r <- populations_r[, -c(1)] 
colnames(populations_r) <- c('Cell_population','pBIC_cell_count','pBIC10_cell_count','pBIC_%','pBIC10_%')
populations_r$Cell_population <- c('Tcells','DC','Bcells_nontumoral','GrMDSC','Bcells_tumoral','Not_classified','Mph','NK','NKT','MoMDSC')

cells <- c("Bcells_tumoral","Bcells_nontumoral","Tcells","DC","GrMDSC","MoMDSC","NK","NKT","Mph","Not_classified")
populations_r <- populations_r  %>%
  mutate(Cell_population =  factor(Cell_population, levels = cells)) %>%
  arrange(Cell_population) 
#write.csv(populations_r,"General_cell_populations.csv")

populations_t <- populations_T
populations_t[1,] <- populations_T[2,]+populations_T[5,]+populations_T[6,]+populations_T[7,]+populations_T[8,]+populations_T[9,]+populations_T[12,]
populations_t[2,] <- populations_T[1,]+populations_T[3,]
populations_t[3,] <- populations_T[10,]
populations_t[4,] <- populations_T[11,]
populations_t[5,] <- populations_T[4,]
populations_t[6,] <- populations_T[13,]
populations_t <- populations_t[-c(7:13), ] 
populations_t <- populations_t[, -c(1)] 
colnames(populations_t) <- c('Cell_population','pBIC_cell_count','pBIC10_cell_count','pBIC_%','pBIC10_%')
populations_t$Cell_population <- c('CD8 central memory','CD8 effector memory','CD8 pre-effector memory','CD8 naive','CD4 effector memory','CD4 Treg')
#write.csv(populations_t,"Clean_Tcell_populations.csv")

Cell_populations_final <- rbind(populations_t,populations_r)
Cell_populations_final_w <- Cell_populations_final[-c(9), ] 
Cell_populations_final_w <- Cell_populations_final_w %>% 
  mutate(`pBIC_%` = pBIC_cell_count/sum(pBIC_cell_count)*100, `pBIC10_%` = pBIC10_cell_count/sum(pBIC10_cell_count)*100)
Cell_populations_final_w[17,] <- c('Total',sum(Cell_populations_final_w$pBIC_cell_count),sum(Cell_populations_final_w$pBIC10_cell_count),sum(Cell_populations_final_w$`pBIC_%`),sum(Cell_populations_final_w$`pBIC10_%`))
Cell_populations_final_w <- Cell_populations_final_w[-c(16.1), ] 
#write.csv(Cell_populations_final_w,"Cell_populations_final.csv")

results_fisher <- Cell_populations_final_w[,c(4,5)]
results_fisher <- results_fisher[-c(16), ]
results_fisher$`pBIC_%` <- as.numeric(results_fisher$`pBIC_%`)
results_fisher$`pBIC10_%` <- as.numeric(results_fisher$`pBIC10_%`)
names_results_fisher <- Cell_populations_final_w$Cell_population[-c(16)]
rownames(results_fisher) <- names_results_fisher
label_mosaic <- c('CD8ce','CD8em','CD8pre-em','CD8naive','CD4em','CD4Treg',"Btumoral","Bnontumoral","DC","GrMDSC","MoMDSC","NK","NKT","Mph","NC")
rownames(results_fisher) <- label_mosaic
chisq.test(results_fisher)
mosaicplot(results_fisher, main='Differences pBIC and pBIC10', dir=c("v","h"), color=TRUE)


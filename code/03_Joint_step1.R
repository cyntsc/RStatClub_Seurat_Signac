#########################################################################
##
## RstatClub_3: Joint RNA and ATAC analysis
##
## Input: Seurat with GEX and ATAC data
## Output: Seurat with joint RNA and ATAC information
##         Stats and coverage plots to illustrate the rna and peak integration
##

## Authors. CSC 
## Code implemented from: 
##
########################################################################

# load libraries
library(Seurat)   
library(Signac)   
library(SeuratDisk) # required to load the pre-defined annotation
library(Matrix)
library(dplyr)
library(here)


here::here()


########################    Initials.  ########################  


dir_plots <- here("plots", "03_Join_v2")
if (!dir.exists(dir_plots)) { dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE) }


## Sample: Flash-Frozen Human Healthy Brain Tissue (3k)

## Load PRE QCed RDS Object
base_name <- 'FFB_Healthy'
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_QCed.rds'))
SeuratOBJ <- readRDS(rds_name)


## ================= Normalize and process the RNA assay first

## normalize the gene expression data using SCTransform, and reduce the dimensionality using PCA
DefaultAssay(SeuratOBJ) <- "RNA"
SeuratOBJ <- SCTransform(SeuratOBJ)
SeuratOBJ <- RunPCA(SeuratOBJ)

print(SeuratOBJ)

# # Alternative normalization with counts
# SeuratOBJ <- NormalizeData(SeuratOBJ) 
# SeuratOBJ <- FindVariableFeatures(SeuratOBJ, nfeatures = 500)
# SeuratOBJ <- ScaleData(SeuratOBJ)
# SeuratOBJ <- RunPCA(SeuratOBJ)

# head(SeuratOBJ[['pca']])
# PC_1          PC_2        PC_3         PC_4         PC_5
# CNTN5      -0.062734245 -0.0044122218  0.01042227  0.091991542  0.017495813
# ZNF385D    -0.069431650 -0.0006593901 -0.03138168  0.092411802  0.006744757

# head(Embeddings(SeuratOBJ, reduction = "pca")[, 1:5])

##  check  to decide which components of the SVD results to use
p1 <- ElbowPlot(SeuratOBJ, ndims = 30, reduction="pca")
p2 <- DepthCor(SeuratOBJ, n = 30, reduction="pca")
plt_rna <- p1 | p2

DefaultAssay(SeuratOBJ) <- "SCT"
SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")
SeuratOBJ@reductions



## ================= Normalize and process the ATAC assay

# process a scATAC-seq dataset, by performing latent semantic indexing (LSI)
DefaultAssay(SeuratOBJ) <- "ATAC"
SeuratOBJ <- FindTopFeatures(SeuratOBJ, min.cutoff = 5) # top feature selection
SeuratOBJ <- RunTFIDF(SeuratOBJ,  method = 1)           # normalization; method = 1 (log-transformation)
SeuratOBJ <- RunSVD(SeuratOBJ)                          # Linear dimension reduction

print(SeuratOBJ)

##  check  to decide which components of the SVD results to use
p1 <- ElbowPlot(SeuratOBJ, ndims = 30, reduction="lsi")
p2 <- DepthCor(SeuratOBJ, n = 30, reduction="lsi")
plt_atac <- p1 | p2

## plot the explained variance per component for pca and lsi

plt_elbow <- plt_rna / plt_atac
png(filename = here(dir_plots, 'elbow_rna_atac.png'))
plt_elbow
dev.off()

SeuratOBJ <- RunUMAP(SeuratOBJ, reduction = "lsi",
                     dims = 2:30, reduction.name = "umap_atac", reduction.key = "UMAPATAC_")
SeuratOBJ@reductions


rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_PCA_LSI_UMAP.rds'))
saveRDS(SeuratOBJ, rds_name)





library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()




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
#library(SeuratDisk) # required to load the pre-defined annotation
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(dplyr)
library(here)


here::here()


########################    Initials.  ########################  


dir_plots <- here("plots", "03_Join_v2")
if (!dir.exists(dir_plots)) { dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE) }


## Sample: Flash-Frozen Human Healthy Brain Tissue (3k)

## Load RDS Object with PCA and LSI reductions
base_name <- 'FFB_Healthy'
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_PCA_LSI_UMAP.rds'))
SeuratOBJ <- readRDS(rds_name)



## ======== Analysis on the RNA (SCT) assay

SeuratOBJ@reductions
DefaultAssay(SeuratOBJ) <- "SCT"

p1 <- DimPlot(SeuratOBJ, reduction = "umap_rna") & NoLegend() #& NoAxes()
head(VariableFeatures(SeuratOBJ))
p2 <- FeaturePlot(SeuratOBJ,
                  c("OBI1-AS1","RANBP3L","ADGRV1","LINC00499"), #cluster 1
                  #c('RNASE1','DBNDD2', 'PPP1R14A', 'CRYAB'),    #cluster 3
                  #c('KCNIP4', 'ROBO2', 'DPP10', 'SLC1A2'),      # variable features
                  reduction = "umap_rna") & NoAxes() & NoLegend()
p1 | p2




## we can do clustering and cluster marker identification


# SeuratOBJ <- FindNeighbors(SeuratOBJ, reduction = "pca",
#                         dims = 1:ncol(Embeddings(SeuratOBJ,"pca")))
# SeuratOBJ <- FindClusters(SeuratOBJ, resolution = 0.2)
# 
# unique(SeuratOBJ$seurat_clusters)
# 
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster1.markers <- FindMarkers(SeuratOBJ, ident.1 = 1, ident.2 = c(0, 2, 3, 4, 5, 6))
# cluster3.markers <- FindMarkers(SeuratOBJ, ident.1 = 3, ident.2 = c(0, 2, 1, 4, 5, 6))
# 
# top_markers <- cluster3.markers %>%
#   filter(avg_log2FC > 1 &
#            p_val_adj < 0.01) %>%
#   top_n(4, wt = p_val_adj)

# 'HMG20B', 'HMG20B', 'PPIF', 'SLC35A4'
# 'FDX2', 'PRKACA', 'LAPTM4B', 'TM6SF1'




## ======== Analysis on the ATAC assay

##  Feature selection
##    select peaks being detected in sufficient number of cells in the data

SeuratOBJ@reductions
DefaultAssay(SeuratOBJ) <- "ATAC"

p3 <- DimPlot(SeuratOBJ, reduction = "umap_atac") & NoLegend() #& NoAxes()
#head(VariableFeatures(SeuratOBJ))
p4 <- FeaturePlot(SeuratOBJ,
                  'chr1-30718015-30718946',
                  #c('chr3-93470143-93471053', 'chr1-30718015-30718946'),
                  reduction = "umap_atac") & NoAxes() & NoLegend()
p3 | p4



## =============  Bi-modal integrative analysis of the RNA-ATAC scMultiome data

## Weighted nearest neighbor analysis
##    It does within-modal and cross-modal prediction by averaging the data of one modality (not the raw data, but the dimension reduced embedding) of its kNNs of the same (within-modal) or the other (cross-modal) modality

SeuratOBJ <- FindMultiModalNeighbors(SeuratOBJ,
                                  reduction.list = list("pca", "lsi"),
                                  dims.list = list(1:ncol(Embeddings(SeuratOBJ,"pca")),
                                                   2:ncol(Embeddings(SeuratOBJ,"lsi"))),
                                  modality.weight.name = c("RNA.weight","ATAC.weight"),
                                  verbose = TRUE)

# NOTES:
#  In this example, we use the two dimension reduction representations (pca and lsi), when integration is done, it should be replaced with the integrated dimension reduction representations (css_rna and css_atac)

# knn.range = 200,  # number of approximate neighbors to compute
# smooth = FALSE,   # Smoothing modality score across each individual modality neighbors


# it generates a NN called "weighted.nn" and stores it at the neighbors slot 
SeuratOBJ@neighbors

# It also generates 2 neighbor graphs called "wknn" and "wsnn" which are both stored at the graphs slot
SeuratOBJ@graphs

## The "weighted.nn" can be then used as input to generate the UMAP embedding
SeuratOBJ <- RunUMAP(SeuratOBJ, nn.name = "weighted.nn", assay = "RNA")

## The "wsnn" graph can be used for clustering

SeuratOBJ <- FindClusters(SeuratOBJ, graph.name = "wsnn", resolution = 0.2)
unique(SeuratOBJ$seurat_clusters)
colnames(SeuratOBJ@meta.data)


DefaultAssay(SeuratOBJ) <- "SCT"

p1 <- UMAPPlot(SeuratOBJ, group.by = "wsnn_res.0.2", label=T) 
p2 <- FeaturePlot(SeuratOBJ,
                  c("OBI1-AS1","RANBP3L","ADGRV1","LINC00499"), #cluster 1
                  #c('RNASE1','DBNDD2', 'PPP1R14A', 'CRYAB'),    #cluster 3
                  #c('KCNIP4', 'ROBO2', 'DPP10', 'SLC1A2'),      # variable features
                  reduction = "umap") & NoAxes() & NoLegend()
p1 | p2 


rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_WNN.rds'))
saveRDS(SeuratOBJ, rds_name)



## we do a rough annotation of cells  

new.cluster.ids <- c("Clust1", "Clust2", "Clust3", "Clust4", "Clust5", "Clust6", "Clust7", "Clust8")
names(new.cluster.ids) <- levels(SeuratOBJ)
SeuratOBJ <- RenameIdents(SeuratOBJ, new.cluster.ids)
unique(Idents(SeuratOBJ))


p1 <- DimPlot(SeuratOBJ, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#p1 <- UMAPPlot(SeuratOBJ, group.by = "wsnn_res.0.2", label=T) & NoAxes()
p2 <- FeaturePlot(SeuratOBJ,
                  c("OBI1-AS1","RANBP3L","ADGRV1","LINC00499"), #cluster 1
                  #c('RNASE1','DBNDD2', 'PPP1R14A', 'CRYAB'),    #cluster 3
                  #c('KCNIP4', 'ROBO2', 'DPP10', 'SLC1A2'),      # variable features
                  order=T,
                  reduction = "umap") & NoAxes() & NoLegend()
p1 | p2



## ====== Linking peaks to genes

DefaultAssay(SeuratOBJ) <- "ATAC"

# Compute the GC content, region lengths, and dinucleotide base frequencies for regions in the assay and add to the feature metadata
SeuratOBJ <- RegionStats(SeuratOBJ, genome = BSgenome.Hsapiens.UCSC.hg38)

colnames(SeuratOBJ@meta.data)

## Find peaks that are correlated with the expression of nearby genes
# For each gene, this function computes the correlation coefficient between the gene expression and accessibility of each peak within a given distance from the gene TSS, and computes an expected correlation coefficient for each peak given the GC content, accessibility, and length of the peak. The expected coefficient values for the peak are then used to compute a z-score and p-value.

SeuratOBJ <- LinkPeaks(
  object = SeuratOBJ,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = c("OBI1-AS1","RANBP3L","ADGRV1","LINC00499"),
  distance = 5e+05, # Distance threshold for peaks to include in regression model
  min.cells = 10, # Minimum number of cells positive for the peak and gene needed to include in the results
  method = "pearson" # spearman
)


##  visualize these links using the CoveragePlot() function

unique(Idents(SeuratOBJ))
idents.plot <- c("Clust1", "Clust2", "Clust3")

p1 <- CoveragePlot(
  object = SeuratOBJ,
  region = "RANBP3L",
  features = "RANBP3L",
  expression.assay = "SCT", # RNA
  idents = idents.plot,
  extend.upstream = 500, # Number of bases to extend the region upstream
  extend.downstream = 10000 # Number of bases to extend the region downstream
)
p1

p2 <- CoveragePlot(
  object = SeuratOBJ,
  region = "RANBP3L",
  features = "RANBP3L",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 10000,
  extend.downstream = 10000
)
p2
p1 / p2

p3 <- CoveragePlot(
  object = SeuratOBJ,
  region = "RANBP3L",
  features = "RANBP3L",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 20000
)
p3
p1 / p2 / p3


p4 <- CoveragePlot(
  object = SeuratOBJ,
  region = "ADGRV1",
  features = "ADGRV1",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 20000
)
p4

p1 / p4


## Save RDS Object

rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_RNA_ATAC_join.rds'))
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat with ATAC clusters saved!')  



## ===== What is next ? TF binding motif enrichment analysis

# The previous analysis can help us to identify putative cis-regulatory elements that are critical for regulating cell type identity or cell state transitions. This is usually achieved via the binding of certain trans-regulators, e.g. TFs, to those open chromatin regions.

## To do that, we also need a database of TF binding motifs: TRANSFAC and JASPAR




library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()




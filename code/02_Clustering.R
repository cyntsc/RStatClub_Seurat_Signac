########################################################################
## Dimensionality reduction for scCellRanger-ARC datasets with Seurat & Signac
##
## Implementation done with public datasets: 
## Samples from : 
##        1) PBMC 10k (ARC 2.0.0) Cell Sorted: 11898 
##        2) PBMC 3k (ARC 2.0.0) Cell Sorted: 2711   
##        3) hippo 42_1 (ARC 2.0.0) Cell Sorted: 3k 
##        4) hippo 42_4 (ARC 2.0.0) Cell Sorted: 3k   

## Authors. CSC / HT
## Date. June 7th, 2023
## Last.Md June 12nd CSC
########################################################################

# load libraries
library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(here)
set.seed(1234)

# Paths in the JHPCE cluster
here::here()


# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/GEX_ATAC_preprocessing/"))) {
  dir.create(here("processed-data/GEX_ATAC_preprocessing/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/GEX_ATAC_preprocessing/"))) {
  dir.create(here("plots/GEX_ATAC_preprocessing/"))
}

# get QC violin plots
plot_violinQC <- function(seuratOBJ, sfeature, stitle) {
  p1 <- VlnPlot(object = seuratOBJ, features = sfeature, 
                group.by = 'orig.ident', pt.size = 0) & geom_boxplot() &
    theme(legend.position = 'none',
          axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
          axis.text.y = element_text(size=8), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) #&
  #labs(title = "", x = 'Samples', y ="")
  ggtitle(stitle)
  return(p1)
}

## plot reductions calculated: pca, umpa, CCA and Harmony
plot_clust <- function(sobj, f_name, reduct, ga2) {
  
  # integrate the samples and clusters
  p1 <- DimPlot(sobj, 
                reduction = reduct, group.by = c("orig.ident", ga2))
  png_file <- paste0(f_name,'_dimplot.png')
  png_name <- here('plots/02_merge_seurats', png_file)  
  ggsave(p1, filename = png_name, height = 5, width = 10)
  
  # visualize the two conditions side-by-side
  p1 <- DimPlot(sobj, 
                reduction = reduct, split.by = "orig.ident")
  png_file <- paste0(f_name,'_dimplot_splitted.png')
  png_name <- here('plots/02_merge_seurats', png_file)  
  ggsave(p1, filename = png_name, height = 5, width = 10)
  
}



########################    Initials.  ########################  

## Sample: Flash-Frozen Human Healthy Brain Tissue (3k)

## Load pre-existing Seurat object

sample_sample <- 'FFB_Healthy'

message('Loading sample: ',sample_sample)

rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_GEX_ATAC.rds'))
SeuratOBJ <- readRDS(rds_name)

message('Seurat with ATAC object loaded!')   

print(SeuratOBJ)




######## Dimensionality reduction

## Perform pre-processing and dimensional reduction on both assays independently
## Standard approaches for RNA and ATAC-seq data.
##    * PCA for GEX (k-means)
##    * LCA for ATAC (spherical clustering)           
##    * WNN for the integrated

## There are two approaches for running the RNA analysis
##    (1) Standard seurat workflow
##    (2) SCTransform 

## SCTransform model normalization is followed by PCA and UMAP dimensionality reduction
##       Function replaces NormalizeData(), ScaleData(), and FindVariableFeatures()


## (1) Standard seurat workflow

DefaultAssay(SeuratOBJ) <- "RNA"

SeuratOBJ <- NormalizeData(SeuratOBJ, 
                           normalization.method = "LogNormalize", 
                           scale.factor = 10000)
SeuratOBJ[["RNA"]]$data

SeuratOBJ <- FindVariableFeatures(SeuratOBJ, 
                                  selection.method = "vst", 
                                  nfeatures = 2000) 
# nfeatures define the top variable features to use
# only used when selection.method is set to 'dispersion' or 'vst'

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SeuratOBJ), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SeuratOBJ, raster=FALSE)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(SeuratOBJ)
SeuratOBJ <- ScaleData(SeuratOBJ, features = all.genes)


## To remove unwanted sources of variation
#SeuratOBJ <- ScaleData(SeuratOBJ, vars.to.regress = "percent.mt")




## Perform linear dimensional reduction

# We perform PCA on the scaled data

SeuratOBJ <- RunPCA(SeuratOBJ, features = VariableFeatures(object = SeuratOBJ))
# Outputs a list of genes with the most positive and negative loadings, representing modules of genes that exhibit either correlation (or anti-correlation) across single-cells in the dataset

SeuratOBJ[['pca']]
head(Embeddings(SeuratOBJ, reduction = "pca")[, 1:5])
head(Stdev(SeuratOBJ, reduction = "pca")[1:5])

##  To visualize both cells and features that define the PCA
DimHeatmap(SeuratOBJ, dims = 1, cells = 300, balanced = TRUE)
# cells and features are ordered according to their PCA scores.


VizDimLoadings(SeuratOBJ, dims = 1:2, reduction = "pca")

#DimPlot(SeuratOBJ)  + NoLegend()

## Determine the ‘dimensionality’ of the dataset
ElbowPlot(SeuratOBJ)


## Cluster the cells

SeuratOBJ <- FindNeighbors(SeuratOBJ, dims = 1:10)

SeuratOBJ <- FindClusters(SeuratOBJ, 
                          resolution = 0.8,
                          algorithm = 1,
                          cluster.name ='C.Lovain')
# Returns a Seurat object where the idents have been updated with new cluster info
# latest clustering results will be stored in object metadata under 'seurat_clusters'. 
# Note that 'seurat_clusters' will be overwritten everytime FindClusters is run

SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:10, reduction = "pca", reduction.name = "umap.lovain")

table(Idents(SeuratOBJ))

colnames(SeuratOBJ@meta.data)
head(SeuratOBJ, n=3)

DimPlot(SeuratOBJ, reduction = "umap.lovain")

SeuratOBJ <- FindClusters(SeuratOBJ, 
                          resolution = 2,
                          algorithm = 1,
                          cluster.name ='C.Lovain.r2')

SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:10, reduction = "pca", reduction.name = "umap.lovain.2")

## run the non-linear dimensional reduction (UMAP/tSNE)

SeuratOBJ@reductions

DimPlot(SeuratOBJ, reduction = "umap.lovain.2")




########### Next run ATAC analysis. ###########

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth

DefaultAssay(SeuratOBJ) <- "ATAC"

SeuratOBJ <- RunTFIDF(SeuratOBJ,
                      method = 1,
                      scale.factor = 10000)

SeuratOBJ <- FindTopFeatures(SeuratOBJ, 
                             min.cutoff = 'q5',
                             verbose = TRUE)
# Set 'q5':  include 95% most common features as the VariableFeatures.
# Set 10:  include features with >10 total counts in the set of VariableFeatures

SeuratOBJ <- RunSVD(SeuratOBJ)

SeuratOBJ <- RunUMAP(SeuratOBJ, 
                     reduction = 'lsi', 
                     dims = 2:50, 
                     reduction.name = "umap.atac", reduction.key = "atacUMAP_")


## Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

SeuratOBJ <- FindMultiModalNeighbors(SeuratOBJ, 
                                     reduction.list = list("pca", "lsi"),
                                     dims.list = list(1:50, 2:50))

SeuratOBJ <- RunUMAP(SeuratOBJ, 
                     nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")

SeuratOBJ <- FindClusters(SeuratOBJ, 
                          graph.name = "wsnn", 
                          algorithm = 3, verbose = FALSE)

SeuratOBJ@reductions

p1 <- DimPlot(SeuratOBJ, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(SeuratOBJ, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(SeuratOBJ, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pALL <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

table(Idents(SeuratOBJ))

## Base name to save plots
base_name <- levels(SeuratOBJ$`orig.ident`[1])

png_file <- paste0(base_name,'_CLuster_ALL.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(pALL, filename = png_name, height = 4, width = 8)


## perform sub-clustering on a specific clusterto find additional structure

SeuratOBJ <- FindSubCluster(SeuratOBJ, 
    cluster = 0,
    graph.name = 'wsnn', #Name of graph to use for the clustering algorithm
    subcluster.name = "sub.cluster", #name of sub cluster added in the meta.data
    resolution = 0.5,
    algorithm = 1) #1 = original Louvain algorithm

Idents(SeuratOBJ) <- "sub.cluster"

colnames(SeuratOBJ@meta.data)

table(SeuratOBJ@meta.data[['seurat_clusters']])

table(SeuratOBJ@meta.data[['sub.cluster']])


# add annotations
SeuratOBJ <- RenameIdents(SeuratOBJ, '0_0' = 'CD14-A', '0_1' ='CD14-B')
SeuratOBJ$celltype <- Idents(SeuratOBJ)


p1 <- DimPlot(SeuratOBJ, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(SeuratOBJ, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(SeuratOBJ, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pALL_sub <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


png_file <- paste0(base_name,'_CLuster_ALL_subclust.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(pALL_sub, filename = png_name, height = 4, width = 8)


# Save RDS Object
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_Clusters.rds'))
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat with ATAC clusters saved!')  




library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()


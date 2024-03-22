########################################################################
##
## RstatClub: Multiome Analysis with Seurat and Signac #2
##
## Input: Truncated H5, meta-data.cvs and fragments.tvs
## Output: rds Seurat objects with ATAC assay, one with QCs and other w/o QCs
##.        GEX and ATAC plots
##
## NOTES: 
## For a ~10k cells it is recommended ~40G free mem to process the TSS() ATAC score.  
## Without storing the base-resolution matrix of integration counts at each site you can use less memory
## For slurm env: $srun --pty --mem=40GB --x11 bash
##
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#packageVersion("Seurat")
library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(here)

here::here()

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/GEX_ATAC_preprocessing/"))) {
    dir.create(here("processed-data/GEX_ATAC_preprocessing/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/GEX_ATAC_preprocessing/"))) {
    dir.create(here("plots/GEX_ATAC_preprocessing/"))
}


########################    Initials.  ########################  

## Link to dataset: https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0 
## Flash-Frozen Human Healthy Brain Tissue (3k)

sample_sample <- 'FFB_Healthy'
#s_tissue <- 'human'

message('Processing sample: ',sample_sample, ' from ', s_tissue, ' brain tissue.')

base_path <- paste0('processed-data/cellranger-arc')    # Read H5 file
meta_path <- here(base_path, "human_brain_3k_per_barcode_metrics.csv")
bc_mtx_path <- here(base_path, "human_brain_3k_filtered_feature_bc_matrix.h5" )
atac_path <- here(base_path,"human_brain_3k_atac_fragments.tsv.gz")

    
# filtered bc mtx
mtx <- Read10X_h5(bc_mtx_path)             # Returns a sparse matrix with rows and columns labeled
rna_counts <- mtx$`Gene Expression`

metadata <- read.csv(file = meta_path, header = TRUE, row.names = 1)
# subset specific fields in the meta data df
meta_tmp = c('atac_peak_region_fragments','atac_fragments')
meta = metadata[meta_tmp]
print(head(meta, n = 3))

SeuratOBJ <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = sample_sample,
    meta.data = meta)

## Data verification
SeuratOBJ
str(SeuratOBJ)
head(SeuratOBJ, n=3)
Idents(SeuratOBJ)

rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(sample_sample,'.rds'))
saveRDS(SeuratOBJ, file = rds_name)

## Read pre-existing Seurat
## pre_existing_seurat <- '/Users/ccardin2/Documents/RStatClub2/processed-data/GEX_ATAC_preprocessing/'
#Seurat2 <- LoadSeuratRds(here('processed-data/GEX_ATAC_preprocessing', paste0(sample_sample,'.rds')))

## Calculate MITO levels
SeuratOBJ$log10GenesPerUMI <- log10(SeuratOBJ$nFeature_RNA) / log10(SeuratOBJ$nCount_RNA)

# Calculate MITO GENES from human
SeuratOBJ$percent.mt <- PercentageFeatureSet(SeuratOBJ, pattern = "^MT-")
SeuratOBJ$MTRatio <- SeuratOBJ$percent.mt / 100 
SeuratOBJ$percent.ribo <- PercentageFeatureSet(SeuratOBJ, pattern = "^RP[LS]")

head(SeuratOBJ, n=3)

message('^MITO and ^RIBO levels processed successfully')



########           Get Visualizations for the GEX           ######## 

# Build and plot:
#       UMI, Genes, MITO and RIBO violin plots
#       Number of cells per sample
#       UMI/transcripts per cell
#       Distribution of genes per cell (histogram)

# Call Vln plot from Seurat suite. (ggplot2 wrapped)    
sfeature <- c('nFeature_RNA', 'nCount_RNA', 'percent.mt')

# Violin plot with UMIs, Genes, ^MT and RIBO levels
p1_GEX <- VlnPlot(SeuratOBJ, features = sfeature, group.by = 'orig.ident')
png_file <- paste0(base_name,'_GEX_QCs.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(p1_GEX, filename = png_name, height = 4, width = 4)

## Plot Genes and UMIs by density per cell 

df_genes_per_cell <- as.data.frame(SeuratOBJ[[]])
df_genes_per_cell %>%
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5)) +
    geom_vline(xintercept = 300) +
    ylab("Log10(UMIs)") +
    xlab("Gene-counts") +
    ggtitle("Genes density by cell") 

df_genes_per_cell %>%
    ggplot(aes(x=orig.ident, y=(nFeature_RNA), fill=orig.ident)) +
    geom_boxplot(alpha = 0.7) +
    theme_classic() +
    theme(axis.text.x = element_text(vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ylab("Log10(nFeature_RNA)") +
    xlab("") +
    ggtitle("Genes distribution by cell")

# Correlation btw Genes/UMIs 
df_genes_per_cell %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt, group.by = 'orig.ident')) + # MTRatio
    #    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MTRatio)) + # MTRatio
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 200, linetype=2) +
    geom_hline(yintercept = 200, linetype=2) +
    facet_wrap(~orig.ident) +
    ylab("log10(nFeature_RNA)") +
    xlab("log10(UMIs)") +
    ggtitle('UMIs/Genes by MT levels')




########  Create ATAC assay and attach it to Seurat object     ##### 

message('Processing ATAC for sample ', levels(SeuratOBJ$`orig.ident`[1]))


## Create gene annotations for hg38 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# show(annotations) / # View(head(annotations,n=5)) /# names(genomeStyles('Homo_sapiens'))

## Format annotation 
seqlevelsStyle(annotations) <- "UCSC"
length(extractSeqlevelsByGroup(species = 'Homo_sapiens', style = 'UCSC', group = 'all')) 
genome(annotations) <- "hg38"
# names(mcols(annotations))

## Base name to save plots
base_name <- paste0(levels(SeuratOBJ$`orig.ident`[1]),'_QC')

## Pull peak fragment information and 
atac_counts <- mtx$Peaks

# Convert a genomic coordinate string to a GRanges object and attach to ATAC counts
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))      
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Create chromatin assay
chrom_assay <- CreateChromatinAssay(counts = atac_counts, 
                                    sep = c(":", "-"), 
                                    fragments = atac_path, 
                                    annotation = annotations)


SeuratOBJ[["ATAC"]] <- chrom_assay

## Annotations of the object are set
Annotation(SeuratOBJ[["ATAC"]]) <- annotations

head(SeuratOBJ, n = 3)

message('Chromatin assay attached to seurat object successfully')



## Measure ATAC quality matrics:  NFR and TSS

DefaultAssay(SeuratOBJ) <- "ATAC"

## Calculate the strength of the nucleosome signal per cell

SeuratOBJ <- NucleosomeSignal(SeuratOBJ)
SeuratOBJ$nucleosome_group <- ifelse(SeuratOBJ$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p1_NS <- FragmentHistogram(object = SeuratOBJ, group.by = 'nucleosome_group')

png_file <- paste0(base_name,'_NS.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(p1_NS, filename = png_name, height = 4, width = 4)

## Calculate the "Transcription Start Site (TSS)" enrichment score

tryCatch( {
    # Important NOTES: https://github.com/stuart-lab/signac/issues/374 
    # 1. Some times jumps an issue of ** memory allocation ** when "fast=FALSE", FALSE argument is mandatory to compute the TSS enrichment scores and visualize with TSSPlot(). You need more memory
    
    SeuratOBJ <- TSSEnrichment(SeuratOBJ, fast = FALSE) 
    # Group by cells with TSS enrichment scores in two groups.
    SeuratOBJ$high.tss <- ifelse(SeuratOBJ$TSS.enrichment > 2, 'High', 'Low')
    
    #colnames(SeuratOBJ@meta.data)
    p1_TSS <- TSSPlot(SeuratOBJ, group.by = 'high.tss') + NoLegend()
    png_file_TSS <- paste0(base_name,'_TSS.png')
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_TSS)
    ggsave(p1_TSS, filename = png_name, height = 4, width = 4)
        
    }
    , error = function(e) { print('An error occurred. Check the annotation or you probably have ATAC quality loss.') } )
    
    ## Add blacklist ratio and fraction of reads in peaks
    SeuratOBJ$blacklist_fraction <- FractionCountsInRegion(
        object = SeuratOBJ,
        assay = 'ATAC',
        regions = blacklist_hg38)
    
    ## Add blacklist ratio and fraction of reads in peaks
    SeuratOBJ$pct_reads_in_peaks <- SeuratOBJ$atac_peak_region_fragments / SeuratOBJ$atac_fragments * 100
    SeuratOBJ$blacklist_ratio <- SeuratOBJ$blacklist_fraction / SeuratOBJ$atac_peak_region_fragments
    #        Error in `x[[i, drop = TRUE]]`:
    #       ! 'blacklist_fraction' not found in this Seurat object
    
    # Plot Peaks in black ratio 
    p1_BlackR <-  VlnPlot(SeuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), group.by = "orig.ident", ncol = 2)  
    png_file_BlackR <- paste0(base_name,'_reads_in_peaks.png')
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_BlackR)
    ggsave(p1_BlackR, filename = png_name, height = 4, width = 5)
    
    ## Plot ATAC main feature scores
    p1_ATAC <- VlnPlot(SeuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), group.by = "orig.ident") 
    png_file_ATAC <- paste0(base_name,'_ATAC_QCs.png')
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_ATAC)
    ggsave(p1_ATAC, filename = png_name, height = 4, width = 7)
    
    message('ATAC QCs plots saved!')  
    

# Save RDS Object
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_GEX_ATAC.rds'))
#saveRDS(SeuratOBJ, file = rds_name)
message('Seurat with ATAC object saved!')   



########  Standard thresholds to remove low quality cells   ##### 

## Apply arbitrary scores as example

SeuratOBJ_QCed <- subset(x = SeuratOBJ,
    subset =  (nCount_RNA > 1000 & nCount_RNA < 10000) &
        (nCount_ATAC > 1000 & nCount_ATAC < 100000)  &
        (nucleosome_signal < 2) & TSS.enrichment > 1)

new_seurat <- paste0(base_name, '_QCed')
SeuratOBJ_QCed$orig.ident <- new_seurat

## Save Seurat QCed data

rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_QCed.rds'))
#saveRDS(SeuratOBJ_QCed, file = rds_name)


message('Seurat QCed saved!')   





############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# > library("sessioninfo")
# > print('Reproducibility information:')
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-03-21 14:49:08 EDT"
# > proc.time()
# user    system   elapsed 
# 1581.802    78.203 10703.538 
# > options(width = 120)
# > session_info()
# 
# backports                     1.4.1       2021-12-13 [2] CRAN (R 4.3.2)
# base64enc                     0.1-3       2015-07-28 [2] CRAN (R 4.3.2)
# beeswarm                      0.4.0       2021-06-01 [2] CRAN (R 4.3.2)
# Biobase                     * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache                 2.10.1      2023-10-26 [2] Bioconductor
# BiocGenerics                * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                      * 1.12.0      2023-10-24 [2] Bioconductor
# BiocManager                   1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocParallel                  1.36.0      2023-10-24 [2] Bioconductor
# biomaRt                       2.58.2      2024-01-30 [2] Bioconductor 3.18 (R 4.3.2)
# Biostrings                  * 2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
# biovizBase                    1.50.0      2023-10-24 [2] Bioconductor
# bit                           4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
# bit64                         4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
# bitops                        1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
# blob                          1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
# BSgenome                    * 1.70.1      2023-11-01 [2] Bioconductor
# BSgenome.Hsapiens.UCSC.hg38 * 1.4.5       2024-03-21 [1] Bioconductor
# cachem                        1.0.8       2023-05-01 [2] CRAN (R 4.3.2)
# callr                         3.7.3       2022-11-02 [2] CRAN (R 4.3.2)
# checkmate                     2.3.1       2023-12-04 [2] CRAN (R 4.3.2)
# cli                           3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
# cluster                       2.1.6       2023-12-01 [3] CRAN (R 4.3.2)
# codetools                     0.2-19      2023-02-01 [3] CRAN (R 4.3.2)
# colorspace                    2.1-0       2023-01-23 [2] CRAN (R 4.3.2)
# cowplot                       1.1.3       2024-01-22 [2] CRAN (R 4.3.2)
# crayon                        1.5.2       2022-09-29 [2] CRAN (R 4.3.2)
# curl                          5.2.0       2023-12-08 [2] CRAN (R 4.3.2)
# data.table                    1.15.0      2024-01-30 [2] CRAN (R 4.3.2)
# DBI                           1.2.1       2024-01-12 [2] CRAN (R 4.3.2)
# dbplyr                        2.4.0       2023-10-26 [2] CRAN (R 4.3.2)
# DelayedArray                  0.28.0      2023-10-24 [2] Bioconductor
# deldir                        2.0-2       2023-11-23 [2] CRAN (R 4.3.2)
# desc                          1.4.3       2023-12-10 [2] CRAN (R 4.3.2)
# devtools                      2.4.5       2022-10-11 [1] CRAN (R 4.3.2)
# dichromat                     2.0-0.1     2022-05-02 [2] CRAN (R 4.3.2)
# digest                        0.6.34      2024-01-11 [2] CRAN (R 4.3.2)
# dotCall64                     1.1-1       2023-11-28 [2] CRAN (R 4.3.2)
# dplyr                       * 1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
# ellipsis                      0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
# EnsDb.Hsapiens.v86          * 2.99.0      2024-03-21 [1] Bioconductor
# ensembldb                   * 2.26.0      2023-10-24 [2] Bioconductor
# evaluate                      0.23        2023-11-01 [2] CRAN (R 4.3.2)
# fansi                         1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
# farver                        2.1.1       2022-07-06 [2] CRAN (R 4.3.2)
# fastDummies                   1.7.3       2023-07-06 [2] CRAN (R 4.3.2)
# fastmap                       1.1.1       2023-02-24 [2] CRAN (R 4.3.2)
# fastmatch                     1.1-4       2023-08-18 [2] CRAN (R 4.3.2)
# filelock                      1.0.3       2023-12-11 [2] CRAN (R 4.3.2)
# fitdistrplus                  1.1-11      2023-04-25 [2] CRAN (R 4.3.2)
# forcats                     * 1.0.0       2023-01-29 [2] CRAN (R 4.3.2)
# foreign                       0.8-86      2023-11-28 [3] CRAN (R 4.3.2)
# Formula                       1.2-5       2023-02-24 [2] CRAN (R 4.3.2)
# fs                            1.6.3       2023-07-20 [2] CRAN (R 4.3.2)
# future                        1.33.1      2023-12-22 [2] CRAN (R 4.3.2)
# future.apply                  1.11.1      2023-12-21 [2] CRAN (R 4.3.2)
# generics                      0.1.3       2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb                * 1.38.5      2023-12-28 [2] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData              1.2.11      2024-02-09 [2] Bioconductor
# GenomicAlignments             1.38.2      2024-01-16 [2] Bioconductor 3.18 (R 4.3.2)
# GenomicFeatures             * 1.54.3      2024-01-31 [2] Bioconductor 3.18 (R 4.3.2)
# GenomicRanges               * 1.54.1      2023-10-29 [2] Bioconductor
# ggbeeswarm                    0.7.2       2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2                     * 3.4.4       2023-10-12 [2] CRAN (R 4.3.2)
# ggrastr                       1.0.2       2023-06-01 [2] CRAN (R 4.3.2)
# ggrepel                       0.9.5       2024-01-10 [2] CRAN (R 4.3.2)
# ggridges                      0.5.6       2024-01-23 [2] CRAN (R 4.3.2)
# globals                       0.16.2      2022-11-21 [2] CRAN (R 4.3.2)
# glue                          1.7.0       2024-01-09 [2] CRAN (R 4.3.2)
# goftest                       1.2-3       2021-10-07 [2] CRAN (R 4.3.2)
# gridExtra                     2.3         2017-09-09 [2] CRAN (R 4.3.2)
# gtable                        0.3.4       2023-08-21 [2] CRAN (R 4.3.2)
# hdf5r                         1.3.9       2024-01-14 [2] CRAN (R 4.3.2)
# here                        * 1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
# Hmisc                         5.1-1       2023-09-12 [2] CRAN (R 4.3.2)
# hms                           1.1.3       2023-03-21 [2] CRAN (R 4.3.2)
# htmlTable                     2.4.2       2023-10-29 [2] CRAN (R 4.3.2)
# htmltools                     0.5.7       2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets                   1.6.4       2023-12-06 [2] CRAN (R 4.3.2)
# httpuv                        1.6.14      2024-01-26 [2] CRAN (R 4.3.2)
# httr                          1.4.7       2023-08-15 [2] CRAN (R 4.3.2)
# ica                           1.0-3       2022-07-08 [2] CRAN (R 4.3.2)
# igraph                        2.0.1.9008  2024-02-09 [2] Github (igraph/rigraph@39158c6)
# IRanges                     * 2.36.0      2023-10-24 [2] Bioconductor
# irlba                         2.3.5.1     2022-10-03 [2] CRAN (R 4.3.2)
# jsonlite                      1.8.8       2023-12-04 [2] CRAN (R 4.3.2)
# KEGGREST                      1.42.0      2023-10-24 [2] Bioconductor
# KernSmooth                    2.23-22     2023-07-10 [3] CRAN (R 4.3.2)
# knitr                         1.45        2023-10-30 [2] CRAN (R 4.3.2)
# labeling                      0.4.3       2023-08-29 [2] CRAN (R 4.3.2)
# later                         1.3.2       2023-12-06 [2] CRAN (R 4.3.2)
# lattice                       0.22-5      2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval                      0.2.2       2019-03-15 [2] CRAN (R 4.3.2)
# leiden                        0.4.3.1     2023-11-17 [2] CRAN (R 4.3.2)
# lifecycle                     1.0.4       2023-11-07 [2] CRAN (R 4.3.2)
# listenv                       0.9.1       2024-01-29 [2] CRAN (R 4.3.2)
# lmtest                        0.9-40      2022-03-21 [2] CRAN (R 4.3.2)
# lubridate                   * 1.9.3       2023-09-27 [2] CRAN (R 4.3.2)
# magrittr                      2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
# MASS                          7.3-60.0.1  2024-01-13 [3] CRAN (R 4.3.2)
# Matrix                        1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics                1.14.0      2023-10-24 [2] Bioconductor
# matrixStats                   1.2.0       2023-12-11 [2] CRAN (R 4.3.2)
# memoise                       2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
# mgcv                          1.9-1       2023-12-21 [3] CRAN (R 4.3.2)
# mime                          0.12        2021-09-28 [2] CRAN (R 4.3.2)
# miniUI                        0.1.1.1     2018-05-18 [2] CRAN (R 4.3.2)
# munsell                       0.5.0       2018-06-12 [2] CRAN (R 4.3.2)
# nlme                          3.1-164     2023-11-27 [3] CRAN (R 4.3.2)
# nnet                          7.3-19      2023-05-03 [3] CRAN (R 4.3.2)
# parallelly                    1.36.0      2023-05-26 [2] CRAN (R 4.3.2)
# patchwork                     1.2.0       2024-01-08 [2] CRAN (R 4.3.2)
# pbapply                       1.7-2       2023-06-27 [2] CRAN (R 4.3.2)
# pillar                        1.9.0       2023-03-22 [2] CRAN (R 4.3.2)
# pkgbuild                      1.4.3       2023-12-10 [2] CRAN (R 4.3.2)
# pkgconfig                     2.0.3       2019-09-22 [2] CRAN (R 4.3.2)
# pkgload                       1.3.4       2024-01-16 [2] CRAN (R 4.3.2)
# plotly                        4.10.4      2024-01-13 [2] CRAN (R 4.3.2)
# plyr                          1.8.9       2023-10-02 [2] CRAN (R 4.3.2)
# png                           0.1-8       2022-11-29 [2] CRAN (R 4.3.2)
# polyclip                      1.10-6      2023-09-27 [2] CRAN (R 4.3.2)
# prettyunits                   1.2.0       2023-09-24 [2] CRAN (R 4.3.2)
# processx                      3.8.3       2023-12-10 [2] CRAN (R 4.3.2)
# profvis                       0.3.8       2023-05-02 [2] CRAN (R 4.3.2)
# progress                      1.2.3       2023-12-06 [2] CRAN (R 4.3.2)
# progressr                     0.14.0      2023-08-10 [2] CRAN (R 4.3.2)
# promises                      1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
# ProtGenerics                  1.34.0      2023-10-24 [2] Bioconductor
# ps                            1.7.6       2024-01-18 [2] CRAN (R 4.3.2)
# purrr                       * 1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
# R6                            2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
# ragg                          1.2.7       2023-12-11 [2] CRAN (R 4.3.2)
# RANN                          2.6.1       2019-01-08 [2] CRAN (R 4.3.2)
# rappdirs                      0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
# RColorBrewer                  1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                          1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
# RcppAnnoy                     0.0.22      2024-01-23 [2] CRAN (R 4.3.2)
# RcppHNSW                      0.6.0       2024-02-04 [2] CRAN (R 4.3.2)
# RcppRoll                      0.3.0       2018-06-05 [1] CRAN (R 4.3.2)
# RCurl                         1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
# readr                       * 2.1.5       2024-01-10 [2] CRAN (R 4.3.2)
# remotes                       2.4.2.1     2023-07-18 [2] CRAN (R 4.3.2)
# reshape2                      1.4.4       2020-04-09 [2] CRAN (R 4.3.2)
# restfulr                      0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
# reticulate                    1.35.0      2024-01-31 [2] CRAN (R 4.3.2)
# rjson                         0.2.21      2022-01-09 [2] CRAN (R 4.3.2)
# rlang                         1.1.3       2024-01-10 [2] CRAN (R 4.3.2)
# rmarkdown                     2.25        2023-09-18 [2] CRAN (R 4.3.2)
# ROCR                          1.0-11      2020-05-02 [2] CRAN (R 4.3.2)
# rpart                         4.1.23      2023-12-05 [3] CRAN (R 4.3.2)
# rprojroot                     2.0.4       2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools                     2.18.0      2023-10-24 [2] Bioconductor
# RSpectra                      0.16-1      2022-04-24 [2] CRAN (R 4.3.2)
# RSQLite                       2.3.5       2024-01-21 [2] CRAN (R 4.3.2)
# rstudioapi                    0.15.0      2023-07-07 [2] CRAN (R 4.3.2)
# rtracklayer                 * 1.62.0      2023-10-24 [2] Bioconductor
# Rtsne                         0.17        2023-12-07 [2] CRAN (R 4.3.2)
# S4Arrays                      1.2.0       2023-10-24 [2] Bioconductor
# S4Vectors                   * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# scales                        1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
# scattermore                   1.2         2023-06-12 [2] CRAN (R 4.3.2)
# sctransform                   0.4.1       2023-10-19 [2] CRAN (R 4.3.2)
# sessioninfo                 * 1.2.2       2021-12-06 [2] CRAN (R 4.3.2)
# Seurat                      * 5.0.1       2023-11-17 [2] CRAN (R 4.3.2)
# SeuratObject                * 5.0.1       2023-11-17 [2] CRAN (R 4.3.2)
# shiny                         1.8.0       2023-11-17 [2] CRAN (R 4.3.2)
# Signac                      * 1.12.9004   2024-03-21 [1] Github (stuart-lab/signac@19bda12)
# sp                          * 2.1-3       2024-01-30 [2] CRAN (R 4.3.2)
# spam                          2.10-0      2023-10-23 [2] CRAN (R 4.3.2)
# SparseArray                   1.2.3       2023-12-25 [2] Bioconductor 3.18 (R 4.3.2)
# spatstat.data                 3.0-4       2024-01-15 [2] CRAN (R 4.3.2)
# spatstat.explore              3.2-6       2024-02-01 [2] CRAN (R 4.3.2)
# spatstat.geom                 3.2-8       2024-01-26 [2] CRAN (R 4.3.2)
# spatstat.random               3.2-2       2023-11-29 [2] CRAN (R 4.3.2)
# spatstat.sparse               3.0-3       2023-10-24 [2] CRAN (R 4.3.2)
# spatstat.utils                3.0-4       2023-10-24 [2] CRAN (R 4.3.2)
# stringi                       1.8.3       2023-12-11 [2] CRAN (R 4.3.2)
# stringr                     * 1.5.1       2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment          1.32.0      2023-10-24 [2] Bioconductor
# survival                      3.5-7       2023-08-14 [3] CRAN (R 4.3.2)
# systemfonts                   1.0.5       2023-10-09 [2] CRAN (R 4.3.2)
# tensor                        1.5         2012-05-05 [2] CRAN (R 4.3.2)
# textshaping                   0.3.7       2023-10-09 [2] CRAN (R 4.3.2)
# tibble                      * 3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                       * 1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect                    1.2.0       2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse                   * 2.0.0       2023-02-22 [2] CRAN (R 4.3.2)
# timechange                    0.3.0       2024-01-18 [2] CRAN (R 4.3.2)
# tzdb                          0.4.0       2023-05-12 [2] CRAN (R 4.3.2)
# urlchecker                    1.0.1       2021-11-30 [2] CRAN (R 4.3.2)
# usethis                       2.2.2       2023-07-06 [2] CRAN (R 4.3.2)
# utf8                          1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
# uwot                          0.1.16      2023-06-29 [2] CRAN (R 4.3.2)
# VariantAnnotation             1.48.1      2023-11-15 [2] Bioconductor
# vctrs                         0.6.5       2023-12-01 [2] CRAN (R 4.3.2)
# vipor                         0.4.7       2023-12-18 [2] CRAN (R 4.3.2)
# viridisLite                   0.4.2       2023-05-02 [2] CRAN (R 4.3.2)
# withr                         3.0.0       2024-01-16 [2] CRAN (R 4.3.2)
# xfun                          0.42        2024-02-08 [2] CRAN (R 4.3.2)
# XML                           3.99-0.16.1 2024-01-22 [2] CRAN (R 4.3.2)
# xml2                          1.3.6       2023-12-04 [2] CRAN (R 4.3.2)
# xtable                        1.8-4       2019-04-21 [2] CRAN (R 4.3.2)
# XVector                     * 0.42.0      2023-10-24 [2] Bioconductor
# yaml                          2.3.8       2023-12-11 [2] CRAN (R 4.3.2)
# zlibbioc                      1.48.0      2023-10-24 [2] Bioconductor
# zoo                           1.8-12      2023-04-13 [2] CRAN (R 4.3.2)
# 
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library

 

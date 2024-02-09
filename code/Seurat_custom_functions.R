
########################################################################
##
## RStatClub
## FUNCTIONS TO HANDLE SEURAT OBJECTS 
##
########################################################################

get_seurat_obj <- function(seuratName, s_bc_mtx, s_tissue, s_meta, b_additional_feat = FALSE) {
    
    # This function creates a Seurat object with rna counts and meta-data attached
    #       @seuratName         string with the name to be assigned to the Seurat object
    #       @s_bc_mtx           string with the filtered barcode mtx file path
    #       @s_tissue           string used to distinguish btw human and mouse gene MITO phenotype.
    #       @s_meta             string with the meta_data file path 
    #       @b_additional_feat  Boolean to required calculate additional features; RIBO and largest genes (******** OPT *********)
    
    # filtered bc mtx
    mtx <- Read10X_h5(s_bc_mtx)             # Returns a sparse matrix with rows and columns labeled
    rna_counts <- mtx$`Gene Expression`
    print(head(rna_counts, n = 3))
    message('GEX counts loaded successfully')
    
    metadata <- read.csv(file = s_meta, header = TRUE, row.names = 1)
    # subset specific fields in the meta data df
    meta_tmp = c('atac_peak_region_fragments','atac_fragments')
    meta = metadata[meta_tmp]
    print(head(meta, n = 3))
    
    seur_obj <- CreateSeuratObject(
        counts = rna_counts,
        assay = "RNA",
        project = seuratName,
        meta.data = meta
    )
    
    message('Meta-data attached successfully')
    
    # Ensure seurat_obj is a Seurat object
    if (!("Seurat" %in% class(seur_obj))) { stop("Seurat object not found") } 
    #Layers(SeuratOBJ[["RNA"]])
    #if (!("counts" %in% Layers(seur_obj[["RNA"]]))) { stop("RNA counts layer not found!") } 
    #if (!("data" %in% Layers(seur_obj[["RNA"]]))) { stop("Data counts not found!") }
    
    # Calculate MITO levels
    seur_obj$log10GenesPerUMI <- log10(seur_obj$nFeature_RNA) / log10(seur_obj$nCount_RNA)
    
    # Select the correct Symbol for each genome. By default is assumed MT, belonging to humans.
    # For mouse genome we use the Symbol: mt-. 
    # https://www.ncbi.nlm.nih.gov/nuccore/34538597 
    
    # MITO GENES for human and mouse genomes. GRCh38 and mm10, respectively
    if (s_tissue=='human') {
        seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")
    } else { # it is mouse
        seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^Mt")
    }
    
    # RIBO GENES for human and mouse genomes. GRCh38 and mm10, respectively    
    if (s_tissue=='human') {
        seur_obj[["percent.ribo"]] <- PercentageFeatureSet(seur_obj, pattern = "^RP[LS]")
    } else { # it is mouse
        seur_obj[["percent.ribo"]] <- PercentageFeatureSet(seur_obj, pattern = "^Rp[ls]")
    }
    seur_obj[["MTRatio"]] <- seur_obj$percent.mt / 100 
    
    message('^MITO and ^RIBO levels processed successfully')
    # Calculate percentages of largest genes by single cell
    # Use the "b_additional_feat" property carefully because of grows the size object significantly.  
    # if (b_additional_feat==TRUE) {        
    #     seur_obj <- l_get_perc_largest_genes(seur_obj)
    #     message('Additional largest genes by cell features attached successfully')
    # }  
    
    message('Seurat completed successfully!')
    return(seur_obj) 
}


# Some descriptive stats for further analysis
get_basic_stats_GEX <- function(seuratOBJ, s_tissue) {
    ## Calculate statistics over the GEX seurat assay: nCount, nFeature, MITO, Ribo, #cells, etc.
    ## INPUT: 
    ##      @seuratOBJ: Seurat obj 
    ##      @s_tissue : string used to distinguish btw human and mouse gene MITO phenotype.    
    ## OUTPUT:
    ##      @tab_stats: table with the statistics about the GEX Seurat Assay
    ##      Print a file with the same data for further analysis
    ##
    
    tryCatch( {
        
        # Initialize metrics
        s_sample = toString(unique(seuratOBJ@meta.data$orig.ident))
        totalCells <- ncol(seuratOBJ)
        # Continuous sample quantile types 4: linear interpolation of the Empirical Cumulative Distribution Function 
        nCount = quantile(seuratOBJ$nCount_RNA, na.rm = TRUE, type = 4)       #  if na.rm = true, any NA and NaN's are removed
        nFeature = quantile(seuratOBJ$nFeature_RNA, na.rm = TRUE, type = 4)   
        #nMT = quantile(seuratOBJ$percent.mt, na.rm = TRUE, type = 4)
        #Calculate deciles and extract the 10th and 90th percentiles of the vector
        nMT = quantile(seuratOBJ$percent.mt, c(.01,.10,.25,.50,.75,.90,.99), na.rm = TRUE, type = 4)
        
        
        # Identify the number if mitochondrial genes by their names starting gene name pattern
        if (s_tissue=='human') {
            mt_genes = length(grep("^MT-",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
        } else { # it is mouse
            mt_genes = length(grep("^mt",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
        }
        
        
        # Identify the number of ribosomal genes, that usually tend to be very highly represented, 
        #           This can vary between cell types,check how prevalent they are in the data.
        if (s_tissue=='human') {
            ribo_genes = length(grep("^RP[LS]",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
        } else {    # it is mouse
            ribo_genes = length(grep("^Rp[ls]",rownames(seuratOBJ@assays$RNA@counts),value = TRUE))
        }
        
        # Build a table with the stats applied and the outputs gotten
        # Table with the cellranger-arc `gene expression statistics`
        #mtx_stats <- matrix(c(nCount, nFeature, nMT), ncol=15, byrow=TRUE)
        mtx_stats <- matrix(c(nCount, nFeature, nMT), ncol=17, byrow=TRUE)
        #class(mtx_stats)
        #dim(mtx_stats)
        colnames(mtx_stats) <- c('nCount0%','nCount25%','nCount50%','nCount75%','nCount100%',
                                 'nGenes0%','nGenes25%','nGenes50%','nGenes75%','nGenes100%',
                                 'nMito0%','nMito10p','nMito25%','nMito50%','nMito75%','nMito90p','nMito99%')
        
        # add additional columns: sample name, total cells, #mito genes and #ribo genes.
        
        #message('nCount, nFeature and nMT statistics done.')
        mtx_stats<-cbind(mtx_stats,mt_genes)    # number of mito genes
        mtx_stats<-cbind(mtx_stats,ribo_genes)  # number of ribo genes  
        sample <- s_sample
        mtx_stats<-cbind(mtx_stats,totalCells)
        mtx_stats<-cbind(mtx_stats,sample)
        # convert to table
        tab_stats <- as.table(mtx_stats)
        print(tab_stats)
        
        return(tab_stats) }
        , error = function(e) {print('An error ocurred. Verify you have an GEX object active') })
}

#### Method 1 ####
# GEX method 1 (Probabilities)
get_seurat_GEX_filteringM1p <- function(seurat_obj, lowRNA_prob = 0.01, highRNA_prob = 0.99, lowFeature_prob = 0.01, 
                                        highMT_prob = 0.90, save_combined=FALSE) {
    # Ensure seurat_obj is a Seurat object
    if(!("Seurat" %in% class(seurat_obj))) stop("seurat_obj should be a Seurat object")
    
    # Check if required slots exist
    required_slots <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
    if(any(!required_slots %in% names(seurat_obj@meta.data))) 
        stop("seurat_obj should have the slots nCount_RNA, nFeature_RNA, and percent.mt")
    
    # Get the name of seurat_obj
    seurat_name <- deparse(substitute(seurat_obj))
    
    # Calculate the thresholds
    countLOW = quantile(seurat_obj$nCount_RNA, probs=lowRNA_prob)
    countHIGH = quantile(seurat_obj$nCount_RNA, probs=highRNA_prob)
    featureLOW = quantile(seurat_obj$nFeature_RNA, probs=lowFeature_prob)
    mitHIGH = round(quantile(seurat_obj$percent.mt, probs=highMT_prob))
    
    # Print the thresholds
    message(paste("Seurat Object: ", seurat_name, 
                  "\nThresholds: ", 
                  "\ncountLOW_UMIs: ", countLOW,"%",
                  "\ncountHIGH_UMIs: ", countHIGH, "%",
                  "\nfeatureLOW_genes: ", featureLOW, "%",
                  "\nmitoHIGH: ", mitHIGH,"%"))
    
    # Subset the Seurat object based on the thresholds
    seuratOBJ.filtered.M1 <- subset(x = seurat_obj, subset = (nFeature_RNA >= featureLOW) &
                                        (nCount_RNA >= countLOW)  &
                                        (nCount_RNA < countHIGH) &
                                        (percent.mt < mitHIGH))
    
    if (save_combined) {
        # Check if processed_data directory exists, if not create it
        if (!dir.exists(here("processed-data"))) {
            dir.create(here("processed-data"))
        }
        
        # Define the file name based on combined parameter
        file_name <- ifelse(save_combined, paste0(seurat_name, ".combined.filtered.GEX.M1.Prob.RDS"), 
                            paste0(seurat_name, ".filtered.GEX.M1.Prob.RDS"))
        
        # Save the filtered Seurat object using here package
        write_rds(seuratOBJ.filtered.M1, here("processed-data", file_name), compress = ('gz'))
        message('Filtered object saved as ', here("processed-data", file_name))
    }      
    
    message('Process completed successfully')
    # Return the filtered Seurat object
    return(seuratOBJ.filtered.M1)
}

#### Method 2 ####
# GEX method 2 (Standard Deviation)
get_seurat_GEX_filteringM2sd <- function(seurat_obj, iSD=2, save_combined=FALSE){
    
    # Ensure seurat_obj is a Seurat object
    if(!("Seurat" %in% class(seurat_obj))) stop("seurat_obj should be a Seurat object")
    
    # Check if required slots exist
    required_slots <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
    if(any(!required_slots %in% names(seurat_obj@meta.data))) stop("seurat_obj should have the slots nCount_RNA, nFeature_RNA, and percent.mt")
    
    # Get the name of seurat_obj
    seurat_name <- deparse(substitute(seurat_obj))
    
    # Calculate count, feature, and mt thresholds
    count.max <- round(mean(seurat_obj$nCount_RNA) + iSD * sd(seurat_obj$nCount_RNA), digits = -2)
    count.min <- max(round(mean(seurat_obj$nCount_RNA) - iSD * sd(seurat_obj$nCount_RNA), digits = -2), 0)
    feat.max <- round(mean(seurat_obj$nFeature_RNA) + iSD * sd(seurat_obj$nFeature_RNA), digits = -2)
    feat.min <- max(round(mean(seurat_obj$nFeature_RNA) - iSD * sd(seurat_obj$nFeature_RNA), digits = -2), 0)
    mt.min <- round(mean(seurat_obj$percent.mt) - iSD * sd(seurat_obj$percent.mt))
    mt.max <- round(mean(seurat_obj$percent.mt) + iSD * sd(seurat_obj$percent.mt))
    
    # Print the thresholds
    message(paste("Seurat Object: ", seurat_name, 
                  "\nThresholds: ", 
                  "\ncount_min_UMIs: ", count.min,
                  "\ncount_max_UMIs: ", count.max, 
                  "\nfeature_min_genes: ", feat.min, 
                  "\nfeature_max_genes: ", feat.max,
                  "\nmt_min: ", mt.min,
                  "\nmt_max: ", mt.max))
    
    # Subset the Seurat object based on the thresholds
    seurat_obj.filtered <- subset(x = seurat_obj, subset = (
        ((nFeature_RNA > feat.min) & (nFeature_RNA < feat.max)) &
            ((nCount_RNA < count.max) & (nCount_RNA > count.min)) & 
            (percent.mt < mt.max))
    )
    
    if (save_combined) {
        # Check if processed_data directory exists, if not create it
        if (!dir.exists(here("processed-data"))) {
            dir.create(here("processed-data"))
        }
        
        # Define the file name based on combined parameter
        file_name <- ifelse(save_combined, paste0(seurat_name, ".combined.filtered.GEX.M2.SD.RData"), 
                            paste0(seurat_name, ".filtered.GEX.M2.SD.RData"))
        
        # Save the filtered Seurat object using here package
        save(seurat_obj.filtered, file = here("processed-data", file_name))
        message('Filtered object saved as ', here("processed-data", file_name))
    }
    
    message('Process completed successfully')
    # Return the filtered Seurat object
    return(seurat_obj.filtered)
}

###STEP 4: DECONVOLUTION



###MUSIC 
run_music <- function(spatial){
  Est.prop_sim = music_prop(bulk.mtx = spatial,  sc.sce = sc , clusters = 'celltypes',
                            samples = 'orig.ident')
  
  music_proportion <- data.matrix(Est.prop_sim$Est.prop.weighted)
  
  return(music_proportion)
}


music_baseline <- run_music(baseline1)
music_baseline2 <- run_music(baseline2)
music_baseline3 <- run_music(baseline3)

music_hard <- run_music(target_hard)
music_medium <- run_music(target_medium)



###CARD

run_CARD_deconvolution <- function(spatial_count) {
  # Assuming the following objects are preloaded and accessible:
  # seurat_object1_filtered: A Seurat object with necessary scRNA-seq data
  # locations: Spatial coordinates or identifiers for spatial data spots
  
  # Create a CARD object with provided spatial counts and other required data from the Seurat object
  CARD_obj_sim = createCARDObject(
    sc_count = counts(sc),
    sc_meta = colData(sc),
    spatial_count = spatial_count,
    spatial_location = layer_manual_MOB,
    ct.varname = "celltypes",
    ct.select = unique(colData(sc)$celltypes),
    sample.varname = NULL
  )
  
  
  
  # Perform deconvolution using the CARD method
  CARD_obj_sim = CARD_deconvolution(CARD_object = CARD_obj_sim)
  
  # Extract and return the proportion matrix from the CARD object
  card_proportion_sim = CARD_obj_sim@Proportion_CARD
  
  return(card_proportion_sim)
}

card_baseline1 <- run_CARD_deconvolution(baseline1)

card_baseline2 <- run_CARD_deconvolution(baseline2)

card_baseline3 <- run_CARD_deconvolution(baseline3)

card_hard <- run_CARD_deconvolution(target_hard)

card_medium <- run_CARD_deconvolution(target_medium)


###CIBERSORTx


calc_mean_expr_by_celltype <- function(seurat_obj) {
  # Ensure that the RNA assay is the default assay
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Split the Seurat object by cell type
  seurat_list <- SplitObject(seurat_obj, split.by = "celltypes")
  
  # Function to calculate mean expression for each gene in a subset
  calc_mean_expr <- function(subset) {
    data <- GetAssayData(subset, assay = "RNA", slot = "data")
    return(rowMeans(as.matrix(data)))
  }
  
  # Apply the function to each cell type subset
  mean_expr_list <- lapply(seurat_list, calc_mean_expr)
  
  # Combine into a single data matrix
  mean_expr_matrix <- do.call(cbind, mean_expr_list)
  
  # Name the columns as cell types
  colnames(mean_expr_matrix) <- names(mean_expr_list)
  
  # Extract gene names from any subset as they should be consistent across all
  gene_names <- rownames(GetAssayData(seurat_list[[1]], assay = "RNA", slot = "data"))
  
  # Assign gene names to the rows
  rownames(mean_expr_matrix) <- gene_names
  
  return(mean_expr_matrix)
}

sig_matrix <- calc_mean_expr_by_celltype(seurat_object1_markers)



cibersort_baseline1 <- IOBR::CIBERSORT(sig_matrix, as.data.frame(baseline1) , QN = FALSE, absolute = FALSE, perm = 10)
cibersort_baseline1<- cibersort_baseline1[, 1:4]


cibersort_baseline2 <- IOBR::CIBERSORT(sig_matrix, as.data.frame(baseline2) , QN = FALSE, absolute = FALSE, perm = 10)
cibersort_baseline2<- cibersort_baseline2[, 1:4]


cibersort_baseline3 <- IOBR::CIBERSORT(sig_matrix, as.data.frame(baseline3) , QN = FALSE, absolute = FALSE, perm = 10)
cibersort_baseline3<- cibersort_baseline3[, 1:4]



cibersort_medium <- IOBR::CIBERSORT(sig_matrix, as.data.frame(target_medium) , QN = FALSE, absolute = FALSE, perm = 10)
cibersort_medium <- cibersort_medium[, 1:4]


cibersort_hard <- IOBR::CIBERSORT(sig_matrix, as.data.frame(target_hard) , QN = FALSE, absolute = FALSE, perm = 10)
cibersort_hard <- cibersort_hard[, 1:4]



####RCTD####

library(spacexr)

run_rctd_analysis <- function(baseline) {
  # Ensure Seurat and required packages are loaded
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is not available. Please install it before running this function.")
  }
  
  # Correcting the cell type names directly in the Seurat object
  seurat_object1_markers@meta.data$celltypes <- gsub('M/TC', 'M.TC', seurat_object1_markers@meta.data$celltypes)
  
  # Extracting required data from Seurat object
  counts <- seurat_object1_markers[["RNA"]]$counts
  cluster <- as.factor(seurat_object1_markers$celltypes)
  names(cluster) <- colnames(seurat_object1_markers)
  nUMI <- seurat_object1_markers$nCount_RNA
  names(nUMI) <- colnames(seurat_object1_markers)
  
  # Creating a Reference object (assumes you have a function or constructor named `Reference`)
  reference <- Reference(counts, cluster, nUMI)
  
  # Prepare coordinates by removing certain columns
  coord <- layer_manual_MOB[, !names(layer_manual_MOB) %in% c("celltypes", "Spot", "ID", "Layer")]
  
  # Constructing a SpatialRNA object and running query (assumes SpatialRNA and create.RCTD functions are available)
  query <- SpatialRNA(coord, baseline, colSums(baseline))
  
  # Running RCTD (Reconstruction of Cell Type Dynamics)
  RCTD <- create.RCTD(query, reference, max_cores = 1, UMI_min_sigma = 200)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  # Extracting results and normalizing the weight matrix
  norm_weights <- t(apply(RCTD@results$weights, 1, function(x) x / sum(x)))
  
  # Convert the normalized weights to a data frame for easier handling
  rctd_baseline <- as.data.frame(norm_weights)
  
  # Return the final data frame
  return(rctd_baseline)
}


rctd_baseline1 <- run_rctd_analysis(baseline1)



rctd_baseline2 <- run_rctd_analysis(baseline3)



rctd_baseline3 <- run_rctd_analysis(baseline3)


rctd_medium <- run_rctd_analysis(target_medium)

rctd_hard <-  run_rctd_analysis(target_hard)




run_music <- function(spatial, sc ){
  
  
  Est.prop_sim = music_prop(bulk.mtx = spatial,  sc.sce = sc , clusters = 'celltypes',
                            samples = 'orig.ident')
  
  music_proportion <- data.matrix(Est.prop_sim$Est.prop.weighted)
  colnames(music_proportion) <- gsub("M/TC", "M.TC", colnames(music_proportion))
  return(music_proportion)
}

seurat_object1_200 <- subset(seurat_object1, features = deg200)
seurat_object1_400 <- subset(seurat_object1, features = deg400)
seurat_object1_600 <- subset(seurat_object1, features = deg600)
seurat_object1_norp <- subset(seurat_object1, features = deg_norp)
seurat_object1_sig_diff <- subset(seurat_object1, features = deg_sigdiff)
seurat_object1_cutoff <- subset(seurat_object1, features = deg_cutoff1.5)

sc <- as.SingleCellExperiment(seurat_object1)

dim(sc)

sc200 <-  as.SingleCellExperiment(seurat_object1_200)

sc400 <-  as.SingleCellExperiment(seurat_object1_400)

sc600 <-  as.SingleCellExperiment(seurat_object1_600)






scnorp <-  as.SingleCellExperiment(seurat_object1_norp)
scsigdiff <-  as.SingleCellExperiment(seurat_object1_sig_diff)
sc_cutoff <-  as.SingleCellExperiment(seurat_object1_cutoff)




music_200 <- run_music(as.matrix(spe@assays$spatial$counts), sc200)
music_400 <- run_music(as.matrix(spe@assays$spatial$counts), sc400)
music_600 <- run_music(as.matrix(spe@assays$spatial$counts), sc600)
music_norp <- run_music(as.matrix(spe@assays$spatial$counts), scnorp)
music_sigdiff <- run_music(as.matrix(spe@assays$spatial$counts), scsigdiff)

music_cutoff <- run_music(as.matrix(spe@assays$spatial$counts), sc_cutoff)


music_full <- run_music(as.matrix(spe@assays$spatial$counts), sc)




###CARD

run_CARD_deconvolution <- function(spatial_count, sc) {
  # Assuming the following objects are preloaded and accessible:
  # seurat_object1_filtered: A Seurat object with necessary scRNA-seq data
  # locations: Spatial coordinates or identifiers for spatial data spots
  
  locations <- strsplit(colnames(spatial_count), "x")
  locations <- matrix(unlist(locations), ncol = 2, byrow = TRUE)
  locations <- apply(locations, 2, as.numeric)
  colnames(locations) <- c("x","y")
  
  rownames(locations) <- colnames(spatial_count)
  locations <- as.data.frame(locations)
  
  locations$x <- round(locations$x)
  locations$y <- round(locations$y)
  
  
  
  # Create a CARD object with provided spatial counts and other required data from the Seurat object
  CARD_obj_sim = createCARDObject(
    sc_count = counts(sc),
    sc_meta = colData(sc),
    spatial_count = spatial_count,
    spatial_location = locations ,
    ct.varname = "celltypes",
    ct.select = unique(colData(sc)$celltypes),
    sample.varname = NULL
  )
  
  
  
  # Perform deconvolution using the CARD method
  CARD_obj_sim = CARD_deconvolution(CARD_object = CARD_obj_sim)
  
  # Extract and return the proportion matrix from the CARD object
  card_proportion_sim = CARD_obj_sim@Proportion_CARD
  colnames(card_proportion_sim) <- gsub("M/TC", "M.TC", colnames(card_proportion_sim))
  return(card_proportion_sim)
}


card_200 <-  run_CARD_deconvolution(as.matrix(spe@assays$spatial$counts), sc200)

card_400 <- run_CARD_deconvolution(as.matrix(spe@assays$spatial$counts), sc400)

card_600 <- run_CARD_deconvolution(as.matrix(spe@assays$spatial$counts), sc600)

card_norp <- run_CARD_deconvolution(as.matrix(spe@assays$spatial$counts), scnorp)

card_sigdiff <- run_CARD_deconvolution(as.matrix(spe@assays$spatial$counts), scsigdiff)


card_cutoff <- run_CARD_deconvolution(as.matrix(spe@assays$spatial$counts), sc_cutoff)

card_full <- run_CARD_deconvolution(as.matrix(spe@assays$spatial$counts), sc)
dim(sc)

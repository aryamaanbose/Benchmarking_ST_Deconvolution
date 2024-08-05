
load("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/Figure3A_layer_annote (1).RData")
layer_to_cell_type <- list(
  GCL = "GC",
  GL = "PGC",
  ONL = "OSNs",
  MCL = "M/TC"
)

layer_manual_MOB$Spot <- rownames(layer_manual_MOB)
layer_manual_MOB$celltypes <- sapply(layer_manual_MOB$Layer, function(x) layer_to_cell_type[[x]])
dim(spe)
dim(layer_manual_MOB)





# Match the row names of layer_manual_MOB with the column names of spe
common_ids <- intersect(colnames(spe), rownames(layer_manual_MOB))

# Ensure both datasets have the common IDs only
spe <- spe[, common_ids]
layer_manual_MOB <- layer_manual_MOB[common_ids, , drop = FALSE]

# Add metadata
spe <- AddMetaData(object = spe, metadata = layer_manual_MOB)
#Idents(spe) <- "Layer"
#SpatialDimPlot(spe, pt.size = 10)


dim(spe)

metadata <- spe@meta.data

intersect(rownames(spe), rownames(layer_manual_MOB))

dim(baseline3)


metadata$spatial_snn_res.0.5 <- as.character(metadata$spatial_snn_res.0.5)
metadata$spatial_snn_res.0.5[metadata$spatial_snn_res.0.5 == "0"] <- "PGC"
metadata$spatial_snn_res.0.5[metadata$spatial_snn_res.0.5 == "1"] <- "GC"
metadata$spatial_snn_res.0.5[metadata$spatial_snn_res.0.5 == "2"] <- "OSNs"
metadata$spatial_snn_res.0.5[metadata$spatial_snn_res.0.5 == "3"] <- "M/TC"


# Ensure all desired cell types are included
cell_types <- unique(metadata$spatial_snn_res.0.5)
samples <- unique(metadata$Spot)

# Create an empty data frame with samples as rows and cell types as columns
identity_matrix <- data.frame(matrix(0, nrow = length(samples), ncol = length(cell_types)))
rownames(identity_matrix) <- samples
colnames(identity_matrix) <- cell_types

# Fill the identity matrix
for (sample in samples) {
  for (cell_type in cell_types) {
    # Check if the sample contains the cell type
    if (any(metadata$Spot == sample & metadata$spatial_snn_res.0.5 == cell_type)) {
      identity_matrix[sample, cell_type] <- 1
    }
  }
}

# Print intermediate identity matrix to debug
print("Intermediate identity matrix with 1s filled:")
print(identity_matrix)

# Modify the identity matrix to set 1s to 0.70 and adjust the remaining values
for (i in 1:nrow(identity_matrix)) {
  num_celltypes <- sum(identity_matrix[i, ] == 1)
  if (num_celltypes > 0) {
    identity_matrix[i, identity_matrix[i, ] == 1] <- 0.30
    identity_matrix[i, identity_matrix[i, ] == 0] <- (1 - 0.30 * num_celltypes) / (length(cell_types) - num_celltypes)
  }
}

# Normalize rows to sum to 1, but only if the row has non-zero entries
for (i in 1:nrow(identity_matrix)) {
  row_sum <- sum(identity_matrix[i, ])
  if (row_sum > 0) {
    identity_matrix[i, ] <- identity_matrix[i, ] / row_sum
  }
}

# Print final identity matrix for verification
print("Final identity matrix after adjustments:")
print(identity_matrix)

rowSums(identity_matrix)
dim(identity_matrix)

# Reorder the columns

desired_column_order <- c("GC", "PGC", "OSNs", "M/TC")
final_identity_matrix <- identity_matrix[, desired_column_order]

datasets <- list(
  final_identity_matrix = "Expected30.csv"
)

rowSums(final_identity_matrix)

# Base directory within your project structure
base_dir <- here("Data", "Processed_BLADE")

# Iterate over the datasets and their file names
for (data_name in names(datasets)) {
  # Construct the file path
  file_path <- file.path(base_dir, datasets[[data_name]])
  
  # Ensure the directory exists
  if (!dir.exists(dirname(file_path))) {
    dir.create(dirname(file_path), recursive = TRUE)
  }
  
  # Write the dataset to a CSV file
  write.csv(get(data_name), file = file_path, row.names = TRUE)
}



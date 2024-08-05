###STEP3: Simulate ST data 

simulate_mixture <- function(proportions, sc, counts_matrix, n_cells = 50) {
  set.seed(123)
  
  # Clone the single-cell object to avoid modifying the original
  MOBSCSIM_sce <- sc
  counts(MOBSCSIM_sce) <- counts_matrix
  
  # Extract the unique cell types from the SingleCellExperiment object
  celltypes <- unique(colData(sc)$celltypes)
  
  # Apply the simulation logic to each row of the proportions matrix
  simulated_mixture <- apply(proportions, 1, function(row_proportions) {
    rowSums(sapply(celltypes, function(celltype) {
      # Find indices of cells that belong to the current cell type
      indices <- sample(which(colData(MOBSCSIM_sce)$celltypes == celltype), size = n_cells, replace = FALSE)
      
      # Calculate the mean counts for cells of the given cell type, scaled by the proportion
      rowMeans(counts_matrix[, indices]) * row_proportions[celltype]
    }))
  })
  
  # Round the counts to the nearest integer since raw read counts are integers
  simulated_mixture <- ceiling(simulated_mixture)
  
  return(simulated_mixture)
}



baseline1 <- simulate_mixture(p1, sc, MOBSC_newcount, n_cells = 50)
baseline2 <- simulate_mixture(p2, sc, MOBSC_newcount, n_cells = 50)
baseline3 <- simulate_mixture(p3, sc, MOBSC_newcount, n_cells = 50)
target_medium <-  simulate_mixture(p4, sc, MOBSC_newcount, n_cells = 50)
target_hard <-  simulate_mixture(p5, sc, MOBSC_newcount, n_cells = 50)




generate_identity_matrix <- function(count_data, sample_name, cell_types = c("GC", "PGC", "OSNs", "M/TC")) {
  
  # Function to create Seurat object and find clusters
  create_seurat_and_cluster <- function(count_data, sample_name) {
    seurat_obj <- CreateSeuratObject(counts = count_data)
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    seurat_obj <- FindNeighbors(seurat_obj)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
    
    # Plot PCA and save the plot
    pca_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "seurat_clusters") + 
      ggtitle(paste("PCA Plot for", sample_name))
    print(pca_plot)
    
    # Save PCA plot as an image
    ggsave(paste0(sample_name, "_PCA_plot.png"), plot = pca_plot)
    
    return(seurat_obj)
  }
  
  # Function to annotate clusters
  annotate_clusters <- function(seurat_obj) {
    metadata <- seurat_obj@meta.data
    metadata$RNA_snn_res.0.2 <- as.character(seurat_obj$seurat_clusters)
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "0"] <- "PGC"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "1"] <- "GC"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "2"] <- "OSNs"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "3"] <- "M/TC"
    return(metadata)
  }
  
  # Create Seurat object and find clusters
  seurat_obj <- create_seurat_and_cluster(count_data, sample_name)
  
  # Annotate clusters
  metadata <- annotate_clusters(seurat_obj)
  
  print(metadata)
  
  # Ensure all desired cell types are included
  samples <- rownames(metadata)
  
  # Create an empty data frame with samples as rows and cell types as columns
  identity_matrix <- data.frame(matrix(0, nrow = length(samples), ncol = length(cell_types)))
  rownames(identity_matrix) <- samples
  colnames(identity_matrix) <- cell_types
  
  # Fill the identity matrix
  for (sample in samples) {
    for (cell_type in cell_types) {
      if (any(rownames(metadata) == sample & metadata$RNA_snn_res.0.2 == cell_type)) {
        identity_matrix[sample, cell_type] <- 1
      }
    }
  }
  
  print("Intermediate identity matrix with 1s filled:")
  print(identity_matrix)
  
  # Modify the identity matrix
  for (i in 1:nrow(identity_matrix)) {
    num_celltypes <- sum(identity_matrix[i, ] == 1)
    if (num_celltypes > 0) {
      identity_matrix[i, identity_matrix[i, ] == 1] <- 0.70
      identity_matrix[i, identity_matrix[i, ] == 0] <- (1 - 0.70 * num_celltypes) / (length(cell_types) - num_celltypes)
    }
  }
  
  # Normalize rows to sum to 1
  for (i in 1:nrow(identity_matrix)) {
    row_sum <- sum(identity_matrix[i, ])
    if (row_sum > 0) {
      identity_matrix[i, ] <- identity_matrix[i, ] / row_sum
    }
  }
  
  print("Final identity matrix after adjustments:")
  print(identity_matrix)
  
  # Reorder the columns
  final_identity_matrix <- identity_matrix[, cell_types]
  
  print("Final identity matrix with desired column order:")
  print(final_identity_matrix)
  
  # Verify row sums and dimensions
  print(rowSums(final_identity_matrix))
  print(dim(final_identity_matrix))
  
  return(final_identity_matrix)
}




identity_matrix_baseline1 <- generate_identity_matrix(baseline1, "baseline1")
identity_matrix_baseline2 <- generate_identity_matrix(baseline2, "baseline1")
identity_matrix_baseline3 <- generate_identity_matrix(baseline3, "baseline1")
identity_matrix_medium <- generate_identity_matrix(target_medium, "baseline1")
identity_matrix_hard <- generate_identity_matrix(target_hard, "baseline1")



rowSums(identity_matrix_baseline3)






datasets <- list(
  identity_matrix_baseline1 = "identity_matrix_baseline1.csv",
  identity_matrix_baseline2 = "identity_matrix_baseline2.csv",
  identity_matrix_baseline3 = "identity_matrix_baseline3.csv",
  identity_matrix_medium = "identity_matrix_medium.csv",
  identity_matrix_hard = "identity_matrix_hard.csv",
  baseline1 = "baseline1.csv",
  baseline2 ="baseline2.csv",
  baseline3 = "baseline3.csv",
  target_medium = "target_medium.csv",
  target_hard = "target_hard.csv"
)

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

custom_colors <- c("1" = "lightblue", "3" = "orange", "2" = "red", "0" = "#9999FF")

custom_colors2 <- c("1" = "lightblue", "2" = "orange", "3" = "red", "0" = "#9999FF")



custom_colors3 <- c("0" = "lightblue", "2" = "orange", "3" = "red", "1" = "#9999FF")






# Function to create Seurat object, find clusters, generate identity matrix, and visualize clusters with t-SNE and PCA
generate_identity_matrix <- function(count_data, sample_name, cell_types = c("GC", "PGC", "OSNs", "M/TC"), perplexity = 20) {
  
  # Function to create Seurat object and find clusters
  create_seurat_and_cluster <- function(count_data, sample_name, perplexity) {
    seurat_obj <- CreateSeuratObject(counts = count_data)
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    seurat_obj <- FindNeighbors(seurat_obj)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
    seurat_obj <- RunTSNE(seurat_obj, perplexity = perplexity)
    
    # Plot t-SNE and save the plot
    tsne_plot <- DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") + 
      ggtitle(paste("t-SNE Plot for", sample_name))
    print(tsne_plot)
    
    # Save t-SNE plot as an image
    ggsave(paste0(sample_name, "_tSNE_plot.png"), plot = tsne_plot)
    
    # Plot PCA
    pca_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "seurat_clusters") + 
      ggtitle(paste("PCA Plot for", sample_name)) +
      scale_color_manual(values = custom_colors)
    
    return(list(seurat_obj = seurat_obj, pca_plot = pca_plot))
  }
  
  # Function to annotate clusters
  annotate_clusters <- function(seurat_obj) {
    metadata <- seurat_obj@meta.data
    metadata$RNA_snn_res.0.2 <- as.character(seurat_obj$seurat_clusters)
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "0" ] <- "PGC"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "2" ] <- "GC"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "3" ] <- "OSNs"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "1" ] <- "M/TC"
    return(metadata)
  }
  
  # Create Seurat object and find clusters
  result <- create_seurat_and_cluster(count_data, sample_name, perplexity)
  seurat_obj <- result$seurat_obj
  pca_plot <- result$pca_plot
  
  # Annotate clusters
  metadata <- annotate_clusters(seurat_obj)
  
  print(metadata)
  
  # Ensure all desired cell types are included
  samples <- rownames(metadata)
  
  # Create an empty data frame with samples as rows and cell types as columns
  identity_matrix <- data.frame(matrix(0, nrow = length(samples), ncol = length(cell_types)))
  rownames(identity_matrix) <- samples
  colnames(identity_matrix) <- cell_types
  
  # Fill the identity matrix
  for (sample in samples) {
    for (cell_type in cell_types) {
      if (any(rownames(metadata) == sample & metadata$RNA_snn_res.0.2 == cell_type)) {
        identity_matrix[sample, cell_type] <- 1
      }
    }
  }
  
  print("Intermediate identity matrix with 1s filled:")
  print(identity_matrix)
  
  # Modify the identity matrix
  for (i in 1:nrow(identity_matrix)) {
    num_celltypes <- sum(identity_matrix[i, ] == 1)
    if (num_celltypes > 0) {
      identity_matrix[i, identity_matrix[i, ] == 1] <- 0.70
      identity_matrix[i, identity_matrix[i, ] == 0] <- (1 - 0.70 * num_celltypes) / (length(cell_types) - num_celltypes)
    }
  }
  
  # Normalize rows to sum to 1
  for (i in 1:nrow(identity_matrix)) {
    row_sum <- sum(identity_matrix[i, ])
    if (row_sum > 0) {
      identity_matrix[i, ] <- identity_matrix[i, ] / row_sum
    }
  }
  
  print("Final identity matrix after adjustments:")
  print(identity_matrix)
  
  # Reorder the columns
  final_identity_matrix <- identity_matrix[, cell_types]
  
  print("Final identity matrix with desired column order:")
  print(final_identity_matrix)
  
  # Verify row sums and dimensions
  print(rowSums(final_identity_matrix))
  print(dim(final_identity_matrix))
  
  # Create the output list containing the identity matrix, the t-SNE plot, and the PCA plot
  output <- list(
    identity_matrix = final_identity_matrix,
    tsne_plot = DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") + 
      ggtitle(paste("t-SNE Plot for", sample_name)),
    pca_plot = pca_plot
  )
  
  return(output)
}

generate_identity_matrix2 <- function(count_data, sample_name, cell_types = c("GC", "PGC", "OSNs", "M/TC"), perplexity = 20) {
  
  # Function to create Seurat object and find clusters
  create_seurat_and_cluster <- function(count_data, sample_name, perplexity) {
    seurat_obj <- CreateSeuratObject(counts = count_data)
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    seurat_obj <- FindNeighbors(seurat_obj)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.2)
    seurat_obj <- RunTSNE(seurat_obj, perplexity = perplexity)
    
    # Plot t-SNE and save the plot
    tsne_plot <- DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") + 
      ggtitle(paste("t-SNE Plot for", sample_name))
    print(tsne_plot)
    
    # Save t-SNE plot as an image
    ggsave(paste0(sample_name, "_tSNE_plot.png"), plot = tsne_plot)
    
    # Plot PCA
    pca_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "seurat_clusters") + 
      ggtitle(paste("PCA Plot for", sample_name)) +
      scale_color_manual(values = custom_colors2)
    
    return(list(seurat_obj = seurat_obj, pca_plot = pca_plot))
  }
  
  # Function to annotate clusters
  annotate_clusters <- function(seurat_obj) {
    metadata <- seurat_obj@meta.data
    metadata$RNA_snn_res.0.2 <- as.character(seurat_obj$seurat_clusters)
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "0" ] <- "PGC"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "2"] <- "GC"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "3" ] <- "OSNs"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "1" ] <- "M/TC"
    return(metadata)
  }
  
  # Create Seurat object and find clusters
  result <- create_seurat_and_cluster(count_data, sample_name, perplexity)
  seurat_obj <- result$seurat_obj
  pca_plot <- result$pca_plot
  
  # Annotate clusters
  metadata <- annotate_clusters(seurat_obj)
  
  print(metadata)
  
  # Ensure all desired cell types are included
  samples <- rownames(metadata)
  
  # Create an empty data frame with samples as rows and cell types as columns
  identity_matrix <- data.frame(matrix(0, nrow = length(samples), ncol = length(cell_types)))
  rownames(identity_matrix) <- samples
  colnames(identity_matrix) <- cell_types
  
  # Fill the identity matrix
  for (sample in samples) {
    for (cell_type in cell_types) {
      if (any(rownames(metadata) == sample & metadata$RNA_snn_res.0.2 == cell_type)) {
        identity_matrix[sample, cell_type] <- 1
      }
    }
  }
  
  print("Intermediate identity matrix with 1s filled:")
  print(identity_matrix)
  
  # Modify the identity matrix
  for (i in 1:nrow(identity_matrix)) {
    num_celltypes <- sum(identity_matrix[i, ] == 1)
    if (num_celltypes > 0) {
      identity_matrix[i, identity_matrix[i, ] == 1] <- 0.70
      identity_matrix[i, identity_matrix[i, ] == 0] <- (1 - 0.70 * num_celltypes) / (length(cell_types) - num_celltypes)
    }
  }
  
  # Normalize rows to sum to 1
  for (i in 1:nrow(identity_matrix)) {
    row_sum <- sum(identity_matrix[i, ])
    if (row_sum > 0) {
      identity_matrix[i, ] <- identity_matrix[i, ] / row_sum
    }
  }
  
  print("Final identity matrix after adjustments:")
  print(identity_matrix)
  
  # Reorder the columns
  final_identity_matrix <- identity_matrix[, cell_types]
  
  print("Final identity matrix with desired column order:")
  print(final_identity_matrix)
  
  # Verify row sums and dimensions
  print(rowSums(final_identity_matrix))
  print(dim(final_identity_matrix))
  
  # Create the output list containing the identity matrix, the t-SNE plot, and the PCA plot
  output <- list(
    identity_matrix = final_identity_matrix,
    tsne_plot = DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") + 
      ggtitle(paste("t-SNE Plot for", sample_name)),
    pca_plot = pca_plot
  )
  
  return(output)
}


generate_identity_matrix3 <- function(count_data, sample_name, cell_types = c("GC", "PGC", "OSNs", "M/TC"), perplexity = 20) {
  
  # Function to create Seurat object and find clusters
  create_seurat_and_cluster <- function(count_data, sample_name, perplexity) {
    seurat_obj <- CreateSeuratObject(counts = count_data)
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    seurat_obj <- FindNeighbors(seurat_obj)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
    seurat_obj <- RunTSNE(seurat_obj, perplexity = perplexity)
    
    # Plot t-SNE and save the plot
    tsne_plot <- DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") + 
      ggtitle(paste("t-SNE Plot for", sample_name)) 
    print(tsne_plot)
    
    # Save t-SNE plot as an image
    ggsave(paste0(sample_name, "_tSNE_plot.png"), plot = tsne_plot)
    
    # Plot PCA
    pca_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "seurat_clusters") + 
      ggtitle(paste("PCA Plot for", sample_name)) +
      scale_color_manual(values = custom_colors3)
    
    return(list(seurat_obj = seurat_obj, pca_plot = pca_plot))
  }
  
  # Function to annotate clusters
  annotate_clusters <- function(seurat_obj) {
    metadata <- seurat_obj@meta.data
    metadata$RNA_snn_res.0.2 <- as.character(seurat_obj$seurat_clusters)
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "1" ] <- "PGC"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "0"] <- "GC"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "3" ] <- "OSNs"
    metadata$RNA_snn_res.0.2[metadata$RNA_snn_res.0.2 == "2" ] <- "M/TC"
    return(metadata)
  }
  
  # Create Seurat object and find clusters
  result <- create_seurat_and_cluster(count_data, sample_name, perplexity)
  seurat_obj <- result$seurat_obj
  pca_plot <- result$pca_plot
  
  # Annotate clusters
  metadata <- annotate_clusters(seurat_obj)
  
  print(metadata)
  
  # Ensure all desired cell types are included
  samples <- rownames(metadata)
  
  # Create an empty data frame with samples as rows and cell types as columns
  identity_matrix <- data.frame(matrix(0, nrow = length(samples), ncol = length(cell_types)))
  rownames(identity_matrix) <- samples
  colnames(identity_matrix) <- cell_types
  
  # Fill the identity matrix
  for (sample in samples) {
    for (cell_type in cell_types) {
      if (any(rownames(metadata) == sample & metadata$RNA_snn_res.0.2 == cell_type)) {
        identity_matrix[sample, cell_type] <- 1
      }
    }
  }
  
  print("Intermediate identity matrix with 1s filled:")
  print(identity_matrix)
  
  # Modify the identity matrix
  for (i in 1:nrow(identity_matrix)) {
    num_celltypes <- sum(identity_matrix[i, ] == 1)
    if (num_celltypes > 0) {
      identity_matrix[i, identity_matrix[i, ] == 1] <- 0.70
      identity_matrix[i, identity_matrix[i, ] == 0] <- (1 - 0.70 * num_celltypes) / (length(cell_types) - num_celltypes)
    }
  }
  
  # Normalize rows to sum to 1
  for (i in 1:nrow(identity_matrix)) {
    row_sum <- sum(identity_matrix[i, ])
    if (row_sum > 0) {
      identity_matrix[i, ] <- identity_matrix[i, ] / row_sum
    }
  }
  
  print("Final identity matrix after adjustments:")
  print(identity_matrix)
  
  # Reorder the columns
  final_identity_matrix <- identity_matrix[, cell_types]
  
  print("Final identity matrix with desired column order:")
  print(final_identity_matrix)
  
  # Verify row sums and dimensions
  print(rowSums(final_identity_matrix))
  print(dim(final_identity_matrix))
  
  # Create the output list containing the identity matrix, the t-SNE plot, and the PCA plot
  output <- list(
    identity_matrix = final_identity_matrix,
    tsne_plot = DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") + 
      ggtitle(paste("t-SNE Plot for", sample_name)),
    pca_plot = pca_plot
  )
  
  return(output)
}




# Example usage with your count data
result_baseline1 <- generate_identity_matrix(baseline1, "Scenario 1", perplexity = 20)
identity_matrix_baseline1 <- result_baseline1$identity_matrix
tsne_plot_baseline1 <- result_baseline1$tsne_plot
pca_plot_baseline1 <- result_baseline1$pca_plot

result_baseline2 <- generate_identity_matrix(baseline2, "Scenario 2", perplexity = 20)
identity_matrix_baseline2 <- result_baseline2$identity_matrix
tsne_plot_baseline2 <- result_baseline2$tsne_plot
pca_plot_baseline2 <- result_baseline2$pca_plot

result_baseline3 <- generate_identity_matrix(baseline3, "Scenario 3", perplexity = 20)
identity_matrix_baseline3 <- result_baseline3$identity_matrix
tsne_plot_baseline3 <- result_baseline3$tsne_plot
pca_plot_baseline3 <- result_baseline3$pca_plot

result_target_medium <- generate_identity_matrix3(target_medium, "Scenario 4", perplexity = 20)
identity_matrix_target_medium <- result_target_medium$identity_matrix
tsne_plot_target_medium <- result_target_medium$tsne_plot
pca_plot_target_medium <- result_target_medium$pca_plot

result_target_hard <- generate_identity_matrix3(target_hard, "Scenario 5", perplexity = 20)
identity_matrix_target_hard <- result_target_hard$identity_matrix
tsne_plot_target_hard <- result_target_hard$tsne_plot
pca_plot_target_hard <- result_target_hard$pca_plot

# Combine PCA plots
combined_pca_plot <- (pca_plot_baseline1 + pca_plot_baseline2 + pca_plot_baseline3 + pca_plot_target_medium + pca_plot_target_hard) +
  plot_layout(ncol = 2, guides = 'collect') &
  plot_annotation(title = "PCA Plots of Spatial Clusters in Synthetic Data Scenarios") &
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

# Label the plots
combined_pca_plot <- combined_pca_plot & plot_annotation(tag_levels = 'A')

# Print the combined PCA plot
print(combined_pca_plot)



coordinates <- do.call(rbind, strsplit(rownames(identity_matrix_medium), "x"))
coordinates <- as.data.frame(coordinates)
colnames(coordinates) <- c("x", "y")
coordinates$x <- as.numeric(as.character(coordinates$x))
coordinates$y <- as.numeric(as.character(coordinates$y))

rownames(coordinates) <- rownames(identity_matrix_medium)



CARD.visualize.pie(
  proportion = identity_matrix_target_hard,
  spatial_location = coordinates, 
  colors = my_colors,
  radius = 0.3)




datasets <- list(
  identity_matrix_baseline1 = "identity_matrix_baseline1.csv",
  identity_matrix_baseline2 = "identity_matrix_baseline2.csv",
  identity_matrix_baseline3 = "identity_matrix_baseline3.csv",
  identity_matrix_target_medium = "identity_matrix_medium.csv",
  identity_matrix_target_hard = "identity_matrix_hard.csv",
  baseline1 = "baseline1.csv"
)

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











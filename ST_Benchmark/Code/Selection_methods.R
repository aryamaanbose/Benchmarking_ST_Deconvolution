
locations <- strsplit(colnames(spe), "x")
locations <- matrix(unlist(locations), ncol = 2, byrow = TRUE)
locations <- apply(locations, 2, as.numeric)
colnames(locations) <- c("x","y")

rownames(locations) <- colnames(spe)
locations <- as.data.frame(locations)

locations$x <- round(locations$x)
locations$y <- round(locations$y)
locations$Location <- rownames(locations)





load("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/Figure3A_layer_annote (1).RData")
layer_to_cell_type <- list(
  GCL = "GC",
  GL = "PGC",
  ONL = "OSNs",
  MCL = "M.TC"
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

dim(seurat_object1)

dim(spe)

#####Comparison of Differential Sets##########
library(pheatmap)

Idents(seurat_object1) <- "celltypes"

seurat_object1 <- subset(seurat_object1, features = rownames(spe))



###Removal of genes ribosomal genes 


Rp_genes <- grep("^Rp", rownames(seurat_object1), value = TRUE)


genes_to_exclude <- c(Rp_genes)

# Invert the exclusion to get a list of genes to keep
genes_to_keep <- setdiff(rownames(seurat_object1), genes_to_exclude)


seurat_object1_norp <- subset(seurat_object1,features = genes_to_keep)



significantly_diff_genes 

# Invert the exclusion to get a list of genes to keep
genes_to_keep <- setdiff(rownames(seurat_object1), significantly_diff_genes)


seurat_object1_sig_diff <- subset(seurat_object1,features = genes_to_keep)






significantly_diff_genes_norp

# Invert the exclusion to get a list of genes to keep
genes_to_keep <- setdiff(rownames(seurat_object1_norp), significantly_diff_genes_norp)

seurat_object_sig_diffnorp <- subset(seurat_object1_norp, features = genes_to_keep)







##########
length(intersect(rownames(seurat_object1), rownames(spe)))

dim(seurat_object1)
dim(spe)


####Markers all 
markers_all <- FindAllMarkers(seurat_object1, only.pos = TRUE, min.pct = 0.1)
markers_all_norp <- FindAllMarkers(seurat_object1_norp, only.pos = TRUE, min.pct = 0.1)
markers_all_sigdiff <-  FindAllMarkers(seurat_object1_sig_diff, only.pos = TRUE, min.pct = 0.1)


markers_all_sigdiff_norp <- FindAllMarkers(seurat_object_sig_diffnorp, only.pos = TRUE, min.pct = 0.1)
markers_all_sigdiff_norp3SD <- FindAllMarkers(seurat_object_sig_diffnorp, only.pos = TRUE, min.pct = 0.1)



####Get top n markers 
get_top_markers <- function(markers_all, total_genes) {
  # Calculate the number of genes per cell type
  genes_per_type <- total_genes %/% 4
  
  # Filter and sort markers for each cell type
  sorted_markers_GC <- markers_all %>%
    filter(cluster == 'GC', avg_log2FC >= 0) %>%
    arrange(p_val_adj) %>%
    head(genes_per_type)
  
  sorted_markers_PGC <- markers_all %>%
    filter(cluster == 'PGC', avg_log2FC >= 0) %>%
    arrange(p_val_adj) %>%
    head(genes_per_type)
  
  sorted_markers_OSNS <- markers_all %>%
    filter(cluster == 'OSNs', avg_log2FC >= 0) %>%
    arrange(p_val_adj) %>%
    head(genes_per_type)
  
  sorted_markers_MTC <- markers_all %>%
    filter(cluster == 'M/TC', avg_log2FC >= 0) %>%
    arrange(p_val_adj) %>%
    head(genes_per_type)
  
  # Extract gene names
  genes_GC <- rownames(sorted_markers_GC)
  print(length(genes_GC))
  genes_PGC <- rownames(sorted_markers_PGC)
  print(length(genes_PGC))
  genes_OSNS <- rownames(sorted_markers_OSNS)
  print(length(genes_OSNS))
  genes_MTC <- rownames(sorted_markers_MTC)
  print(length(genes_MTC))
  
  # Combine all top markers into one list
  print(Reduce(intersect, list(genes_GC, genes_PGC, genes_OSNS, genes_MTC)))
  deg <- unique(c(genes_GC, genes_PGC, genes_OSNS, genes_MTC))
  
  # Return the result
  return(deg)
}
get_top_markersFC_cutoff <- function(markers_all, total_genes, cutoff = 1.5) {
  # Calculate the number of genes per cell type
  genes_per_type <- total_genes %/% 4
  
  # Filter and sort markers for each cell type
  sorted_markers_GC <- markers_all %>%
    filter(cluster == 'GC', avg_log2FC >= cutoff) %>%
    arrange(p_val_adj) %>%
    head(genes_per_type)
  
  sorted_markers_PGC <- markers_all %>%
    filter(cluster == 'PGC', avg_log2FC >= cutoff) %>%
    arrange(p_val_adj) %>%
    head(genes_per_type)
  
  sorted_markers_OSNS <- markers_all %>%
    filter(cluster == 'OSNs', avg_log2FC >= cutoff) %>%
    arrange(p_val_adj) %>%
    head(genes_per_type)
  
  sorted_markers_MTC <- markers_all %>%
    filter(cluster == 'M/TC', avg_log2FC >= cutoff) %>%
    arrange(p_val_adj) %>%
    head(genes_per_type)
  
  
  # Extract gene names
  genes_GC <- rownames(sorted_markers_GC)
  print(length(genes_GC))
  genes_PGC <- rownames(sorted_markers_PGC)
  print(length(genes_PGC))
  genes_OSNS <- rownames(sorted_markers_OSNS)
  print(length(genes_OSNS))
  genes_MTC <- rownames(sorted_markers_MTC)
  print(length(genes_MTC))
  
  # Combine all top markers into one list
  print(Reduce(intersect, list(genes_GC, genes_PGC, genes_OSNS, genes_MTC)))
  deg <- unique(c(genes_GC, genes_PGC, genes_OSNS, genes_MTC))
  
  
  markers <- markers_all[markers_all$cluster == 'GC', ]
  markersPGC <- markers_all[markers_all$cluster == 'PGC', ]
  markersOSNS <- markers_all[markers_all$cluster == 'OSNs', ]
  markersMTC <- markers_all[markers_all$cluster == 'M/TC', ]
  
  create_volcano_plot <- function(data, log2fc_cutoff, pval_cutoff, dataset_name, n_labels = 10) {
    
    # Add new columns to the data
    data <- data %>%
      mutate(logP = -log10(p_val_adj),
             Significant = ifelse(p_val_adj < pval_cutoff & abs(avg_log2FC) > log2fc_cutoff, "Yes", "No"),
             is_Rp = grepl("^Rp", gene)) # Identify Rp genes
    
    # Identify top significant genes
    top_genes <- data %>%
      filter(Significant == "Yes") %>%
      arrange(desc(logP), desc(abs(avg_log2FC))) %>%
      slice_head(n = n_labels) %>%
      mutate(rank = row_number()) # Add ranking numbers
    
    # Modify color mapping
    data <- data %>%
      mutate(Color = case_when(
        Significant == "Yes" & is_Rp == TRUE ~ "Rp",
        Significant == "Yes" & is_Rp == FALSE ~ "Significant",
        TRUE ~ "Non-Significant"
      ))
    
    # Create the volcano plot
    p <- ggplot(data, aes(x = avg_log2FC, y = logP, color = Color)) +
      geom_point(alpha = 0.5) +
      scale_color_manual(values = c("Non-Significant" = "grey", "Significant" = "red", "Rp" = "green")) +
      geom_text(data = top_genes, aes(label = paste(rank, gene, sep = ". ")), vjust = 2, size = 3, check_overlap = TRUE, color = "blue") + # Add ranking numbers to gene labels
      theme_minimal() +
      labs(title = paste("Volcano Plot for", dataset_name), x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
      geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), color = "blue", linetype = "dashed") +
      geom_hline(yintercept = -log10(pval_cutoff), color = "blue", linetype = "dashed")
    
    return(p)
  }
  
  gc_volcano_plot <- create_volcano_plot(markers, log2fc_cutoff = cutoff, pval_cutoff = 0.01, dataset_name = "GC", n_labels = 55)
  PGC_volcano_plot <- create_volcano_plot(markersPGC, log2fc_cutoff = cutoff, pval_cutoff = 0.01, dataset_name = "PGC", n_labels = 55)
  OSN_volcano_plot <- create_volcano_plot(markersOSNS, log2fc_cutoff = cutoff, pval_cutoff = 0.01, dataset_name = "OSNs", n_labels = 55)
  MTC_volcano_plot <- create_volcano_plot(markersMTC, log2fc_cutoff = cutoff, pval_cutoff = 0.01, dataset_name = "M/TC", n_labels = 55)
  
  plots <- gc_volcano_plot | PGC_volcano_plot |  OSN_volcano_plot | MTC_volcano_plot
  
  
  # Return the result
  return(deg)
}




deg200 <- get_top_markers(markers_all, 220)
deg400 <- get_top_markers(markers_all, 440)
deg600 <- get_top_markers(markers_all, 680)


deg_norp <-  get_top_markers(markers_all_norp, 220)
deg_cutoff1.5 <- get_top_markersFC_cutoff(markers_all, 220)
deg_sigdiff <- get_top_markers(markers_all_sigdiff, 220)


deg_cutoff1 <- get_top_markersFC_cutoff(markers_all, 220)
deg_cutoff0.5 <- get_top_markersFC_cutoff(markers_all, 220)
deg_sigdiff <- get_top_markers(markers_all_sigdiff, 220)
deg_sigdiff_norp <- get_top_markers(markers_all_sigdiff_norp, 220)
deg_sigdiff_norp_3SD <- get_top_markers(markers_all_sigdiff_norp3SD, 220)
autogene <- readLines(here("Spatial_Blade" ,"selected_genes.txt"))
autogene200 <- readLines(here("Spatial_Blade" ,"selected_genes200.txt"))
autogene400 <- readLines(here("Spatial_Blade" ,"selected_genes400.txt"))
autogene400_2 <- readLines(here("Spatial_Blade" ,"selected_genes400_2.txt"))
autogene200_2 <- readLines(here("Spatial_Blade" ,"selected_genes200_2.txt"))
autogene400_3 <- readLines(here("Spatial_Blade" ,"selected_genes400_3.txt"))

autogene600 <- readLines(here("Spatial_Blade" ,"selected_genes600.txt"))



deg_sigdiff_cutoff <- get_top_markers(markers_all_sigdiff,220)


deg_sigdiff_cutoff2 <- get_top_markersFC_cutoff(markers_all_sigdiff,220)






























process_and_save_geneset <- function(seurat_object, geneset, geneset_name, base_dir = here("Data", "Processed_BLADE")) {
  # Subset the Seurat object
  seurat_object_filtered <- subset(seurat_object, features = geneset)
  
  # Define helper functions
  calc_std_dev <- function(seurat_obj) {
    data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
    return(rowSds(as.matrix(data)))
  }
  
  calc_mean_expr <- function(seurat_obj) {
    data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
    return(rowMeans(as.matrix(data)))
  }
  
  calc_variance <- function(seurat_obj) {
    data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
    return(rowVars(as.matrix(data)))
  }
  
  calc_corrected_std <- function(mean_matrix, variance_matrix) {
    gene_list <- rownames(mean_matrix)
    ct_list <- colnames(mean_matrix)
    New_std <- matrix(nrow = length(gene_list), ncol = length(ct_list))
    
    for (i in seq_along(ct_list)) {
      trend <- fitTrendVar(mean_matrix[, i], variance_matrix[, i])$trend
      New_std[, i] <- sqrt(trend(mean_matrix[, i]))
    }
    
    rownames(New_std) <- gene_list
    colnames(New_std) <- ct_list
    return(New_std)
  }
  
  # Split the Seurat object by cell type
  seurat_list <- SplitObject(seurat_object_filtered, split.by = "celltypes")
  
  # Apply the function to each cell type subset
  std_dev_list <- lapply(seurat_list, calc_std_dev)
  std_dev_matrix <- do.call(cbind, std_dev_list)
  colnames(std_dev_matrix) <- names(std_dev_list)
  
  mean_expr_list <- lapply(seurat_list, calc_mean_expr)
  mean_expr_matrix <- do.call(cbind, mean_expr_list)
  colnames(mean_expr_matrix) <- names(mean_expr_list)
  
  variance_list <- lapply(seurat_list, calc_variance)
  variance_matrix <- do.call(cbind, variance_list)
  colnames(variance_matrix) <- names(variance_list)
  
  gene_names <- rownames(GetAssayData(seurat_object_filtered, assay = "RNA", slot = "data"))
  rownames(std_dev_matrix) <- gene_names
  rownames(mean_expr_matrix) <- gene_names
  rownames(variance_matrix) <- gene_names
  
  New_std <- calc_corrected_std(mean_expr_matrix, variance_matrix)
  
  # Save the datasets
  datasets <- list(
    mean_expr_matrix = paste0("deg_mean", geneset_name, ".csv"),
    New_std = paste0("deg_sd", geneset_name, ".csv")
  )
  
  base_dir_full <- here(base_dir)
  
  for (data_name in names(datasets)) {
    file_path <- file.path(base_dir_full, datasets[[data_name]])
    
    if (!dir.exists(dirname(file_path))) {
      dir.create(dirname(file_path), recursive = TRUE)
    }
    
    write.csv(get(data_name), file = file_path, row.names = TRUE)
    message(paste("Saved", datasets[[data_name]]))
  }
  
  return(list(mean_expr_matrix = mean_expr_matrix, New_std = New_std))
}


results200 <- process_and_save_geneset(seurat_object1, deg200, "200")
results400 <- process_and_save_geneset(seurat_object1, deg400, "400")
results600 <- process_and_save_geneset(seurat_object1, deg600, "600")
results_norp <- process_and_save_geneset(seurat_object1_norp, deg_norp, "_norp")
autogene200 <- process_and_save_geneset(seurat_object1, autogene200_2, "_ag200")
autogene200_2 <- process_and_save_geneset(seurat_object1, autogene200_2, "_ag200_2") 
autogene400 <- process_and_save_geneset(seurat_object1, autogene400_2, "_ag400")
autogene200 <- process_and_save_geneset(seurat_object1, autogene400_2, "_ag400_2")
autogene400_3 <- process_and_save_geneset(seurat_object1, autogene400_3, "_ag400_3")
results_cutoff <-  process_and_save_geneset(seurat_object1, deg_cutoff1.5, "_cutoff1.5")
results_cutoff1 <-  process_and_save_geneset(seurat_object1, deg_cutoff1, "_cutoff1")
results_cutoff0.5 <- process_and_save_geneset(seurat_object1, deg_cutoff0.5, "_cutoff0.5")
results_sigdiff <-  process_and_save_geneset(seurat_object1_sig_diff, deg_sigdiff, "_sigdiff200")
results_sigdiff_norp <-  process_and_save_geneset(seurat_object_sig_diffnorp, deg_sigdiff_norp, "_sigdiff_norp200")
results_sigdiff_norp3SD <-  process_and_save_geneset(seurat_object_sig_diffnorp, deg_sigdiff_norp_3SD, "_sigdiff_norp_3SD_200")


results_sigdiff_cutoff <-  process_and_save_geneset(seurat_object1_sig_diff, deg_sigdiff_cutoff, "_sigdiff_cutoff")


results_sigdiff_cutoff2 <-  process_and_save_geneset(seurat_object1_sig_diff, deg_sigdiff_cutoff2, "_sigdiff_cutoff2")



autogene600 <- process_and_save_geneset(seurat_object1, autogene600, "_ag600")














compute_and_plot_correlations <- function(expression_matrix, geneset_name) {
  # Compute the correlation matrix
  corr_matrix <- cor(expression_matrix)
  
  # Melt the matrix for ggplot
  melted_corr_matrix <- melt(corr_matrix)
  
  # Create a ggplot object for the correlation matrix
  plot <- ggplot(melted_corr_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low= "blue", high  = "red", midpoint = 0.5,
                         limits = c(0, 1), oob = scales::squish) +
    labs(title = paste("Correlation Matrix for", geneset_name)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank())
  
  return(plot)
}

# List of gene sets with their expression matrices
gene_sets <- list(
  "200" = results200$mean_expr_matrix,
  "ag200" = autogene200$mean_expr_matrix
)

# Compute and collect plots
plots <- lapply(names(gene_sets), function(geneset_name) {
  compute_and_plot_correlations(gene_sets[[geneset_name]], geneset_name)
})

# Combine all plots into one figure
combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)




dim(spe)

####Important for real data analyisis please run
calculate_ari_and_purity <- function(proportions_sim, true_labels_df) {
  # Helper function to calculate Node Purity
  calculateNodePurity <- function(true_labels, predicted_labels) {
    combined <- table(true_labels, predicted_labels)
    sum_max_per_cluster <- sum(apply(combined, 2, max))  # Max per column (predicted label)
    total_predictions <- sum(combined)
    purity <- sum_max_per_cluster / total_predictions
    return(purity)
  }
  
  cell_type_names <- as.factor(colnames(proportions_sim))
  
  # Determine dominant cell types for each spot
  dominant_idx <- apply(proportions_sim, 1, which.max)
  dominant_types <- cell_type_names[dominant_idx]
  names(dominant_types) <- rownames(proportions_sim)
  
  # Create a dataframe for the results
  df_sim <- data.frame(Spot = names(dominant_types), DominantType = dominant_types, stringsAsFactors = FALSE)
  
  # Ensure the dataframe is ordered correctly
  df_sim <- df_sim[order(df_sim$Spot),]
  
  # Reorder df_sim to match the order of spots in true_labels_df
  df_sim_ordered <- df_sim[match(true_labels_df$Spot, df_sim$Spot),]
  
  # Verify alignment
  if (!all(true_labels_df$Spot == df_sim_ordered$Spot)) {
    stop("Spot alignment issue in the dataset")
  }
  
  # Extract predicted and true labels
  predicted_labels <- df_sim_ordered$DominantType
  true_labels <- true_labels_df$celltypes
  
  # Calculate ARI and Node Purity
  ARI <- mclust::adjustedRandIndex(true_labels, predicted_labels)
  Node_Purity <- calculateNodePurity(true_labels, predicted_labels)
  
  return(list(ARI = ARI, Node_Purity = Node_Purity))
}


process_dataset <- function(dataset_name, true_labels_df, base_dir = here("Spatial_Blade", "Results_new")) {
  # Load the dataset
  proportion_matrix <- read.csv(here(base_dir, dataset_name))
  print(proportion_matrix)
  print(dim(proportion_matrix))
  # Ensure rownames match the column names of true_labels_df
  rownames(proportion_matrix) <- colnames(spe)
  
  # Calculate ARI and Node Purity
  results <- calculate_ari_and_purity(proportion_matrix, true_labels_df)
  
  # Print and return the results
  print(paste("ARI for", dataset_name, ":", results$ARI))
  print(paste("Node Purity for", dataset_name, ":", results$Node_Purity))
  
  
  
  correlation <- calculate_correlations(spe, marker_genes, proportion_matrix)
  
  results$GC <- correlation$GC
  results$PGC <- correlation$PGC
  results$MTC <- correlation$M.TC
  results$OSN <- correlation$OSNs
  
  
  plot_dominant_ct <- function(proportions_sim, layer_manual_MOB, point_size = 3) {
    # Common plot theme
    common_theme <- theme_minimal() + 
      theme(legend.position = "bottom")
    
    # Common labels
    common_labs <- labs(color = "Cell Type / Layer")
    
    # Calculate the dominant indices and types
    dominant_idx <- apply(proportions_sim, 1, which.max)
    cell_type_names <- colnames(proportions_sim)
    dominant_types <- cell_type_names[dominant_idx]
    names(dominant_types) <- rownames(proportions_sim)
    
    # Create the data frame for plotting
    df_dominant <- data.frame(
      Spot = names(dominant_types),
      DominantType = dominant_types,
      stringsAsFactors = FALSE
    )
    
    # Add the x and y coordinates directly from layer_manual_MOB
    if (!all(c("x", "y") %in% colnames(layer_manual_MOB))) {
      stop("The input data frame 'layer_manual_MOB' must contain 'x' and 'y' columns.")
    }
    
    coordinates <- layer_manual_MOB[, c("x", "y")]
    coordinates$Spot <- rownames(layer_manual_MOB)
    
    # Merge coordinates with dominant cell types
    df_dominant <- merge(df_dominant, coordinates, by = "Spot", all.x = TRUE)
    
    # Plotting with ggplot2
    ggplot(df_dominant, aes(x = x, y = y, color = DominantType)) +
      geom_point(size = point_size) + 
      common_theme +
      labs(title = "Dominant Cell Types") +
      common_labs +
      scale_color_manual(values = c("lightblue",
                                    "orange",
                                    "red",
                                    "#9999FF"))  # Adjust colors as needed
  }
  
  results$plot <- plot_dominant_ct(proportion_matrix, layer_manual_MOB, point_size = 3)
  
  my_colors <- c( "lightblue", "orange", "red", "#9999FF")
  
  
  results$prop_matrix <- proportion_matrix
  
  return(results)
}




###for R based methods
process_dataset2 <- function(dataset,dataset_name, true_labels_df) {
  # Load the dataset
  proportion_matrix <- dataset 
  rownames(proportion_matrix) <- colnames(spe)
  # Calculate ARI and Node Purity
  results <- calculate_ari_and_purity(proportion_matrix, true_labels_df)
  
  # Print and return the results
  print(paste("ARI for", dataset_name, ":", results$ARI))
  print(paste("Node Purity for", dataset_name, ":", results$Node_Purity))
  
  
  
  correlation <- calculate_correlations(spe, marker_genes, proportion_matrix)
  
  results$GC <- correlation$GC
  results$PGC <- correlation$PGC
  results$MTC <- correlation$M.TC
  results$OSN <- correlation$OSNs
  
  
  plot_dominant_ct <- function(proportions_sim, layer_manual_MOB, point_size = 5) {
    # Common plot theme
    common_theme <- theme_minimal() + 
      theme(legend.position = "bottom")
    
    # Common labels
    common_labs <- labs(color = "Cell Type / Layer")
    
    # Calculate the dominant indices and types
    dominant_idx <- apply(proportions_sim, 1, which.max)
    cell_type_names <- colnames(proportions_sim)
    dominant_types <- cell_type_names[dominant_idx]
    names(dominant_types) <- rownames(proportions_sim)
    
    # Create the data frame for plotting
    df_dominant <- data.frame(
      Spot = names(dominant_types),
      DominantType = dominant_types,
      stringsAsFactors = FALSE
    )
    
    # Add the x and y coordinates directly from layer_manual_MOB
    if (!all(c("x", "y") %in% colnames(layer_manual_MOB))) {
      stop("The input data frame 'layer_manual_MOB' must contain 'x' and 'y' columns.")
    }
    
    coordinates <- layer_manual_MOB[, c("x", "y")]
    coordinates$Spot <- rownames(layer_manual_MOB)
    
    # Merge coordinates with dominant cell types
    df_dominant <- merge(df_dominant, coordinates, by = "Spot", all.x = TRUE)
    
    # Plotting with ggplot2
    ggplot(df_dominant, aes(x = x, y = y, color = DominantType)) +
      geom_point(size = point_size) + 
      common_theme +
      labs(title = "Dominant Cell Types") +
      common_labs +
      scale_color_manual(values = c("lightblue",
                                    "orange",
                                    "red",
                                    "#9999FF"))  # Adjust colors as needed
  }
  
  results$plot <- plot_dominant_ct(proportion_matrix, layer_manual_MOB, point_size = 3)
  
  my_colors <- c( "lightblue", "orange", "red", "#9999FF")
  
  prop <- CARD.visualize.pie(
    proportion = proportion_matrix,
    spatial_location = locations, 
    colors = my_colors, 
    radius = 0.4)
  
  results$prop <- prop
  
  
  
  results$prop_matrix <- proportion_matrix
  
  return(results)
}




results_BLADE1 <- process_dataset("deg200.csv", layer_manual_MOB)
results_BLADE2 <- process_dataset("deg400.csv", layer_manual_MOB)
results_BLADE3 <- process_dataset("deg600.csv", layer_manual_MOB)



sp1 <-  process_dataset("", layer_manual_MOB)



results_card200 <- process_dataset2(card_200,"CARD 200", layer_manual_MOB)
results_card400 <- process_dataset2(card_400, "CARD 400",layer_manual_MOB)
results_card600 <- process_dataset2(card_600, "CARD 600",layer_manual_MOB)



results_music200 <- process_dataset2(music_200,"MUSIC 200", layer_manual_MOB)
results_music400 <- process_dataset2(music_400, "MUSIC 400",layer_manual_MOB)
results_music600 <- process_dataset2(music_600,"MUSIC 600", layer_manual_MOB)



results_BLADE_norp <- process_dataset("deg_norp.csv", layer_manual_MOB)
results_CARDnorp  <-  process_dataset2(card_norp,"CARD No RP", layer_manual_MOB)
results_MUSICnorp <- process_dataset2(music_norp,"Music No RP", layer_manual_MOB)


results_BLADE_sigdiff <- process_dataset("deg_sigdiff200.csv", layer_manual_MOB)
results_CARD_sigdiff  <-  process_dataset2(card_sigdiff,"CARD No sigdiff", layer_manual_MOB)
results_MUSIC_sigdiff <- process_dataset2(music_sigdiff,"Music No sigdiff", layer_manual_MOB)



results_CARD_cutoff <- process_dataset2(card_cutoff,"L2FC cutoff 1.5 CARD",  layer_manual_MOB)
results_BLADE_cutoff <- process_dataset("deg_cutoff1.5.csv", layer_manual_MOB)
results_MUSIC_cutoff <- process_dataset2(music_cutoff,"L2FC cutoff 1.5 MUSIC", layer_manual_MOB)


bladesp200 <-  process_dataset("BLADE_sp_200.csv", layer_manual_MOB)
bladesp400 <-  process_dataset("BLADE_sp_400.csv", layer_manual_MOB)
bladesp600 <-  process_dataset("BLADE_sp_600.csv", layer_manual_MOB)




bladespnorp <-  process_dataset("BLADE_sp_norp.csv", layer_manual_MOB)
bladespsigdiff <-  process_dataset("BLADE_sp_sigdiff.csv", layer_manual_MOB)
bladespfc <-  process_dataset("BLADE_sp_fc.csv", layer_manual_MOB)








results_BLADE_norp$plot +
  results_BLADE_cutoff$plot













results_BLADE4 <- process_dataset("deg_ag200.csv", layer_manual_MOB)
results_BLADE5 <- process_dataset("deg_ag400.csv", layer_manual_MOB)
results_BLADE600 <- process_dataset("deg_ag600.csv", layer_manual_MOB)









results_card200 <- process_dataset2(card_200,"CARD 200", layer_manual_MOB)
results_card400 <- process_dataset2(card_400, "CARD 400",layer_manual_MOB)
results_card600 <- process_dataset2(card_600, "CARD 600",layer_manual_MOB)


results_music200 <- process_dataset2(music_200,"MUSIC 200", layer_manual_MOB)
results_music400 <- process_dataset2(music_400, "MUSIC 400",layer_manual_MOB)
results_music600 <- process_dataset2(music_600,"MUSIC 600", layer_manual_MOB)



results_BLADE_norp <- process_dataset("deg_norp.csv", layer_manual_MOB)
results_CARDnorp  <-  process_dataset2(card_norp,"CARD No RP", layer_manual_MOB)
results_MUSICnorp <- process_dataset2(music_norp,"Music No RP", layer_manual_MOB)






results_90_manual <- process_dataset("deg200_sp2.csv", layer_manual_MOB)
results_90_SP <- process_dataset("deg_sp3.csv", layer_manual_MOB)
results_50_SP <- process_dataset("deg_sp4.csv", layer_manual_MOB)
results_30_SP <- process_dataset("deg_sp5.csv", layer_manual_MOB)











results_BLADE_sigdiff <- process_dataset("deg_sigdiff200.csv", layer_manual_MOB)
results_CARD_sigdiff  <-  process_dataset2(card_sigdiff,"CARD No sigdiff", layer_manual_MOB)
results_MUSIC_sigdiff <- process_dataset2(music_sigdiff,"Music No sigdiff", layer_manual_MOB)



results_CARD_cutoff <- process_dataset2(card_cutoff,"L2FC cutoff 1.5 CARD",  layer_manual_MOB)
results_BLADE_cutoff <- process_dataset("deg_cutoff1.5.csv", layer_manual_MOB)
results_MUSIC_cutoff <- process_dataset2(music_cutoff,"L2FC cutoff 1.5 MUSIC", layer_manual_MOB)






results_card_full <- process_dataset2(card_full, "CARD Full",layer_manual_MOB)
results_music_full <-  process_dataset2(music_full, "MUSIC Full",layer_manual_MOB)





results_card_full$plot
results_




results_BLADE7 <- process_dataset("deg_ag200_2.csv", layer_manual_MOB)
results_BLADE8 <- process_dataset("deg_ag400_2.csv", layer_manual_MOB)
results_BLADE10 <- process_dataset("deg_cutoff1.5.csv", layer_manual_MOB)



results_BLADE11 <- process_dataset("deg_ag400_3.csv", layer_manual_MOB)
results_BLADE12 <- process_dataset("deg_cutoff1.csv", layer_manual_MOB)
results_BLADE13 <- process_dataset("deg_cutoff0.5.csv", layer_manual_MOB)
results_BLADE9 <- process_dataset("deg200_sp2.csv", layer_manual_MOB)
results_BLADE14 <- process_dataset("deg_sp3.csv", layer_manual_MOB)






results_card200 <- process_dataset2(card_200,"CARD 200", layer_manual_MOB)
results_cardnorp <- process


results_card400 <- process_dataset2(card_400, "CARD 400",layer_manual_MOB)
results_card600 <- process_dataset2(card_600, "CARD 600",layer_manual_MOB)


results_music200 <- process_dataset2(music_200,"MUSIC 200", layer_manual_MOB)
results_music400 <- process_dataset2(music_400, "MUSIC 400",layer_manual_MOB)
results_music600 <- process_dataset2(music_600,"MUSIC 600", layer_manual_MOB)

results_card_full <- process_dataset2(card_full, "CARD Full",layer_manual_MOB)
results_music_full <-  process_dataset2(music_full, "MUSIC Full",layer_manual_MOB)

results_card_full$plot
results_music_full$plot


results_CARD_cutoff <- process_dataset2(card_cutoff,"L2FC cutoff 1.5 CARD",  layer_manual_MOB)
results_BLADE_cutoff <- process_dataset("deg_cutoff1.5.csv", layer_manual_MOB)
results_MUSIC_cutoff <- process_dataset2(music_cutoff,"L2FC cutoff 1.5 MUSIC", layer_manual_MOB)






results_music200 <- process_dataset2(music_200, layer_manual_MOB)
results_card400 <- process_dataset2(card_400, layer_manual_MOB)
results_music400 <- process_dataset2(music_400, layer_manual_MOB)
results_BLADE15 <- process_dataset("deg_sp4.csv", layer_manual_MOB)
results_BLADE16 <- process_dataset("deg_sigdiff200.csv", layer_manual_MOB)
results_BLADE17 <- process_dataset("deg_sigdiff_norp200.csv", layer_manual_MOB)
results_BLADE18 <- process_dataset("deg_sigdiff_norp_3SD_200.csv", layer_manual_MOB)
results_BLADE19 <- process_dataset("deg_sp5.csv", layer_manual_MOB)
results_BLADE20<- process_dataset("deg_group.csv", layer_manual_MOB)
results_BLADE21<- process_dataset("deg_group_norp.csv", layer_manual_MOB)






results_BLADE9 <- process_dataset("deg200_sp2.csv", layer_manual_MOB)
results_BLADE14 <- process_dataset("deg_sp3.csv", layer_manual_MOB)
results_BLADE15 <- process_dataset("deg_sp4.csv", layer_manual_MOB)
results_BLADE19 <- process_dataset("deg_sp5.csv", layer_manual_MOB)




results_BLADE6 <- process_dataset("deg_norp.csv", layer_manual_MOB)
results_BLADE16 <- process_dataset("deg_sigdiff200.csv", layer_manual_MOB)
results_BLADE17 <- process_dataset("deg_sigdiff_norp200.csv", layer_manual_MOB)
results_BLADE18 <- process_dataset("deg_sigdiff_norp_3SD_200.csv", layer_manual_MOB)





results_BLADE20<- process_dataset("deg_group.csv", layer_manual_MOB)
results_BLADE21<- process_dataset("deg_group_norp.csv", layer_manual_MOB)



results_BLADE22<- process_dataset("deg_sig_diff_cutoff.csv", layer_manual_MOB)



results_BLADE23<- process_dataset("deg_sig_diff_cutoff2.csv", layer_manual_MOB)
























performance_metrics <- data.frame(
  Method = c(
    "deg200", "deg400", "deg600", "deg_ag200", "deg_ag400", 
    "deg_norp", "deg_ag200_2", "deg_ag400_2", "deg200_sp2", "deg_cutoff1.5", 
    "deg_ag400_3", "deg_sp3", "CARD 200", "CARD 400", "MUSIC 200", "MUSIC 400", 
    "deg_sp4", "deg_sigdiff200", "deg_sigdiff_norp200", "deg_sigdiff_norp_3SD_200"
  ),
  ARI = c(
    results_BLADE1$ARI, results_BLADE2$ARI, results_BLADE3$ARI, results_BLADE4$ARI, results_BLADE5$ARI,
    results_BLADE6$ARI, results_BLADE7$ARI, results_BLADE8$ARI, results_BLADE9$ARI, results_BLADE10$ARI,
    results_BLADE11$ARI, results_BLADE14$ARI, results_CARD200$ARI, results_card400$ARI, results_music200$ARI, results_music400$ARI,
    results_BLADE15$ARI, results_BLADE16$ARI, results_BLADE17$ARI, results_BLADE18$ARI
  ),
  Node_Purity = c(
    results_BLADE1$Node_Purity, results_BLADE2$Node_Purity, results_BLADE3$Node_Purity, results_BLADE4$Node_Purity, results_BLADE5$Node_Purity,
    results_BLADE6$Node_Purity, results_BLADE7$Node_Purity, results_BLADE8$Node_Purity, results_BLADE9$Node_Purity, results_BLADE10$Node_Purity,
    results_BLADE11$Node_Purity, results_BLADE14$Node_Purity, results_CARD200$Node_Purity, results_card400$Node_Purity, results_music200$Node_Purity, results_music400$Node_Purity,
    results_BLADE15$Node_Purity, results_BLADE16$Node_Purity, results_BLADE17$Node_Purity, results_BLADE18$Node_Purity
  ),
  GC = c(
    results_BLADE1$GC, results_BLADE2$GC, results_BLADE3$GC, results_BLADE4$GC, results_BLADE5$GC,
    results_BLADE6$GC, results_BLADE7$GC, results_BLADE8$GC, results_BLADE9$GC, results_BLADE10$GC,
    results_BLADE11$GC, results_BLADE14$GC, results_CARD200$GC, results_card400$GC, results_music200$GC, results_music400$GC,
    results_BLADE15$GC, results_BLADE16$GC, results_BLADE17$GC, results_BLADE18$GC
  ),
  OSN = c(
    results_BLADE1$OSN, results_BLADE2$OSN, results_BLADE3$OSN, results_BLADE4$OSN, results_BLADE5$OSN,
    results_BLADE6$OSN, results_BLADE7$OSN, results_BLADE8$OSN, results_BLADE9$OSN, results_BLADE10$OSN,
    results_BLADE11$OSN, results_BLADE14$OSN, results_CARD200$OSN, results_card400$OSN, results_music200$OSN, results_music400$OSN,
    results_BLADE15$OSN, results_BLADE16$OSN, results_BLADE17$OSN, results_BLADE18$OSN
  ),
  PGC = c(
    results_BLADE1$PGC, results_BLADE2$PGC, results_BLADE3$PGC, results_BLADE4$PGC, results_BLADE5$PGC,
    results_BLADE6$PGC, results_BLADE7$PGC, results_BLADE8$PGC, results_BLADE9$PGC, results_BLADE10$PGC,
    results_BLADE11$PGC, results_BLADE14$PGC, results_CARD200$PGC, results_card400$PGC, results_music200$PGC, results_music400$PGC,
    results_BLADE15$PGC, results_BLADE16$PGC, results_BLADE17$PGC, results_BLADE18$PGC
  ),
  MTC = c(
    results_BLADE1$MTC, results_BLADE2$MTC, results_BLADE3$MTC, results_BLADE4$MTC, results_BLADE5$MTC,
    results_BLADE6$MTC, results_BLADE7$MTC, results_BLADE8$MTC, results_BLADE9$MTC, results_BLADE10$MTC,
    results_BLADE11$MTC, results_BLADE14$MTC, results_CARD200$MTC, results_card400$MTC, results_music200$MTC, results_music400$MTC,
    results_BLADE15$MTC, results_BLADE16$MTC, results_BLADE17$MTC, results_BLADE18$MTC
  ),
  Total = c(
    results_BLADE1$total, results_BLADE2$total, results_BLADE3$total, results_BLADE4$total, results_BLADE5$total,
    results_BLADE6$total, results_BLADE7$total, results_BLADE8$total, results_BLADE9$total, results_BLADE10$total,
    results_BLADE11$total, results_BLADE14$total, results_CARD200$total, results_card400$total, results_music200$total, results_music400$total,
    results_BLADE15$total, results_BLADE16$total, results_BLADE17$total, results_BLADE18$total
  )
)




performance_metrics_long <- melt(performance_metrics, id.vars = "Method")

# Define a color palette
colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(performance_metrics_long, aes(x = Method, y = variable, fill = value)) +
  geom_tile(color = "white") +  # Use white to outline each tile
  scale_fill_gradientn(colors = colors, name = "Metric Value") +  # Gradient fill based on values
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +  # Add text labels to each tile
  theme_minimal() +  # Minimal theme
  labs(title = "Performance Metrics Heatmap", x = "Method", y = "Metric") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Slant text labels for x-axis
    axis.text.y = element_text(size = 10, face = "bold"),  # Bold y-axis text
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Bold and centered title
  )





print(heatmap_plot)












performance_metrics <- data.frame(
  Method = c("deg200", "deg400", "deg600", "deg_ag200", "deg_ag400", "deg_norp", "deg_ag200_2", "deg_ag400_2", "deg200_sp", "deg_cutoff", "deg_ag400_3", "Spatial Clusters 90%", "CARD 200 DE", "CARD 400 DE", "MUSIC 200 DE", "MUSIC 400 DE", "Spatial Clusters 50%", "No Significantly Different Genes BLADE 200", "No Significantly Different Genes BLADE 200 NO RP", "No Significantly Different Genes BLADE 200 NO RP +3SD"),
  ARI = c(results_BLADE1$ARI, results_BLADE2$ARI, results_BLADE3$ARI, results_BLADE4$ARI, results_BLADE5$ARI, results_BLADE6$ARI, results_BLADE7$ARI, results_BLADE8$ARI, results_BLADE9$ARI, results_BLADE10$ARI, results_BLADE11$ARI, results_BLADE14$ARI, results_CARD200$ARI, results_card400$ARI, results_music200$ARI, results_music400$ARI, results_BLADE15$ARI, results_BLADE16$ARI, results_BLADE17$ARI,results_BLADE18$ARI ),
  Purity = c(results_BLADE1$Node_Purity, results_BLADE2$Node_Purity, results_BLADE3$Node_Purity, results_BLADE4$Node_Purity, results_BLADE5$Node_Purity,results_BLADE6$Node_Purity, results_BLADE7$Node_Purity, results_BLADE8$Node_Purity, results_BLADE9$Node_Purity, results_BLADE10$Node_Purity, results_BLADE11$Node_Purity, results_BLADE14$Node_Purity, results_CARD200$Node_Purity, results_card400$Node_Purity, results_music200$Node_Purity, results_music400$Node_Purity, results_BLADE15$Node_Purity, results_BLADE16$Node_Purity, results_BLADE17$Node_Purity, results_BLADE18$Node_Purity)
)

# Reshape the data frame to long format
performance_metrics_long <- melt(performance_metrics, id.vars = "Method")

colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)


# Create the heatmap using ggplot2
heatmap_plot_combine <- ggplot(performance_metrics_long, aes(x = Method, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = colors, limits = c(0, 1), name = "Value") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 5) +
  theme_minimal() +
  labs(title = "Performance Metrics Heatmap", x = "Gene Selection Method", y = "Metric", size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display the heatmap
print(heatmap_plot_combine)



create_funky_heatmap <- function(data) {
  # Arrange the data based on ARI
  data <- data %>%
    arrange(desc(ARI))
  
  # Define the column metadata
  column_info <- tribble(
    ~id,     ~group,         ~name,                      ~geom,        ~palette,    ~options,
    "Method", "",            "Method",                   "text",       NA,          list(hjust = 0, width = 6),
    "ARI",   "Cluster Accuracy",      "Adjusted Rand Index",      "bar",        "palette1",  list(width = 4, legend = FALSE, text = TRUE),
    "Purity", "Cluster Accuracy",     "Node Purity",              "bar",        "palette1",  list(width = 4, legend = FALSE, text = TRUE),
    "Total",  "Correllation",       "Total",                    "funkyrect",  "palette2",  list(text = TRUE),
    "GC",    "Correllation",        "GC",                       "funkyrect",  "palette2",  list(text = TRUE),
    "OSN",   "Correllation",        "OSN",                      "funkyrect",  "palette2",  list(text = TRUE),
    "PGC",   "Correllation",        "PGC",                      "funkyrect",  "palette2",  list(text = TRUE),
    "MTC",   "Correllation",        "MTC",                      "funkyrect",  "palette2",  list(text = TRUE)
  )
  
  # Create the funky heatmap
  funky_heatmap(data, column_info = column_info, expand = list(xmax = 4))
}





performance_metrics <- data.frame(
  Method = c("DEG 200", "DEG 400", "DEG 600", "AG 200", "AG 400"),
  ARI = c(results_BLADE1$ARI, results_BLADE2$ARI, results_BLADE3$ARI, results_BLADE4$ARI, results_BLADE5$ARI),
  Purity = c(results_BLADE1$Node_Purity, results_BLADE2$Node_Purity, results_BLADE3$Node_Purity, results_BLADE4$Node_Purity, results_BLADE5$Node_Purity),
  Total = c(
    results_BLADE1$total, results_BLADE2$total, results_BLADE3$total, results_BLADE4$total, results_BLADE5$total
  ),
  GC = c(
    results_BLADE1$GC, results_BLADE2$GC, results_BLADE3$GC, results_BLADE4$GC, results_BLADE5$GC
  ),
  OSN = c(
    results_BLADE1$OSN, results_BLADE2$OSN, results_BLADE3$OSN, results_BLADE4$OSN, results_BLADE5$OSN
  ),
  PGC = c(
    results_BLADE1$PGC, results_BLADE2$PGC, results_BLADE3$PGC, results_BLADE4$PGC, results_BLADE5$PGC
  ),
  MTC = c(
    results_BLADE1$MTC, results_BLADE2$MTC, results_BLADE3$MTC, results_BLADE4$MTC, results_BLADE5$MTC
  ))



create_funky_heatmap(performance_metrics)








# Reshape the data frame to long format
performance_metrics_long <- melt(performance_metrics, id.vars = "Method")

# Create the heatmap using ggplot2
heatmap_plot1 <- ggplot(performance_metrics_long, aes(x = Method, y = variable, fill = value)) +
  geom_tile(color = "black")+
  geom_text(aes(label = round(value, 2)), color = "black", size = 5) +
  theme_minimal() +
  labs(title = "Performance Metrics Heatmap", x = "Gene Selection Method", y = "Metric", size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display the heatmap
print(heatmap_plot1)






performance_metrics <- data.frame(
  Method = c("DEG 0.5 Cutoff", "DEG 1 Cutoff", "DEG 1.5 Cutoff", "AG 200 Distance", "AG 400 Distance", "AG 400 Distance 2"),
  ARI = c(results_BLADE13$ARI, results_BLADE12$ARI, results_BLADE10$ARI, results_BLADE7$ARI, results_BLADE8$ARI, results_BLADE11$ARI),
  Purity = c(results_BLADE13$Node_Purity, results_BLADE12$Node_Purity, results_BLADE10$Node_Purity, results_BLADE7$Node_Purity, results_BLADE8$Node_Purity,results_BLADE11$Node_Purity )
)




# Reshape the data frame to long format
performance_metrics_long <- melt(performance_metrics, id.vars = "Method")

# Create the heatmap using ggplot2
heatmap_plot2 <- ggplot(performance_metrics_long, aes(x = Method, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = colors, limits = c(0, 1), name = "Value") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 5) +
  theme_minimal() +
  labs(title = "Performance Metrics Heatmap", x = "Gene Selection Method", y = "Metric", size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display the heatmap
print(heatmap_plot2)







performance_metrics <- data.frame(
  Method = c("CARD 200", "CARD 400", "Music 200", "Music 400"),
  ARI = c(results_CARD200$ARI, results_card400$ARI, results_music200$ARI, results_music400$ARI),
  Purity = c(results_CARD200$Node_purity, results_card400$Node_Purity, results_music200$Node_Purity, results_music400$Node_Purity)
)




# Reshape the data frame to long format
performance_metrics_long <- melt(performance_metrics, id.vars = "Method")

# Create the heatmap using ggplot2
heatmap_plot3 <- ggplot(performance_metrics_long, aes(x = Method, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = colors, limits = c(0, 1), name = "Value") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 5) +
  theme_minimal() +
  labs(title = "Performance Metrics Heatmap", x = "Gene Selection Method", y = "Metric", size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display the heatmap
print(heatmap_plot3)



###Create heatmap for 
































deg <- deg[deg %in% rownames(spe)]
autogene <- autogene[autogene %in% rownames(spe)]



overlap_genes <- intersect(deg, autogene)



# Calculate the percentage of overlap
percentage_overlap_deg <- length(overlap_genes) / length(deg) * 100
percentage_overlap_autogene <- length(overlap_genes) / length(autogene) * 100





create_heatmap <- function(seurat_obj, genes, group_by, title) {
  # Subset data for the genes of interest
  gene_data <- FetchData(seurat_obj, vars = genes)
  
  # Add cluster annotations
  gene_data$cluster <- seurat_obj@meta.data[[group_by]]
  
  # Calculate average expression per cluster
  avg_exp <- gene_data %>% group_by(cluster) %>% summarise_all(mean, na.rm = TRUE)
  
  # Convert to matrix for heatmap
  avg_exp_matrix <- as.matrix(avg_exp[,-1])
  rownames(avg_exp_matrix) <- avg_exp$cluster
  
  # Generate heatmap
  p <- pheatmap(avg_exp_matrix, 
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                main = title, 
                scale = "row", 
                silent = TRUE,
                show_colnames = FALSE)
  return(p)
}

# Create heatmaps for the two gene sets
heatmap_deg <- create_heatmap(spe, deg, "Layer", "Differentially expressed Genes ")
heatmap_autogene <- create_heatmap(spe, autogene, "Layer", "AutoGene Genes")

# Display heatmaps side by side
grid.arrange(heatmap_deg[[4]], heatmap_autogene[[4]], ncol = 2)



ggplot(data = layer_manual_MOB, aes(x = x, y = y, color = celltypes)) +
  geom_point(aes(colour = celltypes), size = 5) +
  common_theme +
  labs(title = "Mouse olfactory bulb annotated layers", color = "Layer") +
  common_labs +
  scale_color_manual(values = c("GC" = "lightblue",    # lighter red
                                "M.TC" = "orange",  # lighter green
                                "OSNs" = "red",  # lighter blue
                                "PGC" = "#9999FF"))


























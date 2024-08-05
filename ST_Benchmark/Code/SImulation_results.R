library(scales) 
library(spdep)



load("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/Figure3A_layer_annote (1).RData")
layer_to_cell_type <- list(
  GCL = "GC",
  GL = "PGC",
  ONL = "OSNs",
  MCL = "M.TC"
)

layer_manual_MOB$Spot <- rownames(layer_manual_MOB)
layer_manual_MOB$celltypes <- sapply(layer_manual_MOB$Layer, function(x) layer_to_cell_type[[x]])










align_data_frames <- function(reference, result) {
  numeric_columns <- function(df) {
    df[sapply(df, is.numeric)]
  }
  
  if (!is.data.frame(reference) || !is.data.frame(result)) {
    stop("Either reference or result is not a data frame.")
  }
  
  reference <- numeric_columns(reference)
  result <- numeric_columns(result)
  
  common_columns <- intersect(colnames(reference), colnames(result))
  if (length(common_columns) == 0) {
    stop("No common numeric columns found.")
  }
  
  reference <- reference[, common_columns, drop = FALSE]
  result <- result[, common_columns, drop = FALSE]
  
  common_rows <- intersect(rownames(reference), rownames(result))
  if (length(common_rows) == 0) {
    stop("No common rows found.")
  }
  
  reference <- reference[common_rows, , drop = FALSE]
  result <- result[common_rows, , drop = FALSE]
  
  return(list(reference = reference, result = result))
}

calculate_ari_and_purity <- function(proportions_sim, proportions_ref) {
  # Helper function to calculate dominant and secondary indices and the count of non-zero entries
  results <- list() 
  
  calculateNodePurity <- function(true_labels, predicted_labels) {
    combined <- table(true_labels, predicted_labels)
    sum_max_per_cluster <- sum(apply(combined, 2, max))  # Max per column (predicted label)
    total_predictions <- sum(combined)
    purity <- sum_max_per_cluster / total_predictions
    return(purity)
  }
  
  
  get_types <- function(proportions) {
    dominant_idx <- apply(proportions, 1, function(x) order(x, decreasing = TRUE)[1])
    secondary_idx <- apply(proportions, 1, function(x) order(x, decreasing = TRUE)[2])
    cell_type_names <- colnames(proportions)
    
    list(
      Dominant = cell_type_names[dominant_idx],
      Secondary = cell_type_names[secondary_idx],
      Count = rowSums(proportions > 0)  # Count non-zero entries for each spot
    )
  }
  types_sim <- get_types(proportions_sim)
  types_ref <- get_types(proportions_ref)
  
  # Prepare data frames ensuring consistent structure and data types
  df_sim <- data.frame(
    Spot = rownames(proportions_sim),
    DominantType = types_sim$Dominant,
    SecondaryType = types_sim$Secondary,
    Count = types_sim$Count,
    Set = rep("Simulated", length(rownames(proportions_sim))),
    stringsAsFactors = FALSE
  )
  
  df_ref <- data.frame(
    Spot = rownames(proportions_ref),
    DominantType = types_ref$Dominant,
    SecondaryType = types_ref$Secondary,
    Count = types_ref$Count,
    Set = rep("Reference", length(rownames(proportions_ref))),
    stringsAsFactors = FALSE
  )
  
  # Merge data frames based on 'Spot' to ensure matching entries
  df_merged <- merge(df_sim, df_ref, by = "Spot", suffixes = c(".sim", ".ref"))
  
  # Calculate ARI and Node Purity for matched entries
  ARI_dominant <- adjustedRandIndex(df_merged$DominantType.ref, df_merged$DominantType.sim)
  ARI_secondary <- adjustedRandIndex(df_merged$SecondaryType.ref, df_merged$SecondaryType.sim)
  
  purity_dominant <- calculateNodePurity(df_merged$DominantType.ref, df_merged$DominantType.sim)
  purity_secondary <- calculateNodePurity(df_merged$SecondaryType.ref, df_merged$SecondaryType.sim)
  
  results$ARI <- ARI_dominant
  results$ARI2 <- ARI_secondary
  
  results$Node_Purity <- purity_dominant
  results$Node_Purity2 <- purity_secondary
  # Return a list of ARI and Purity values
  
  
  
  
  
  return(results)
}


calculateNodePurity <- function(true_labels, predicted_labels) {
  combined <- table(true_labels, predicted_labels)
  sum_max_per_cluster <- sum(apply(combined, 2, max))  # Max per column (predicted label)
  total_predictions <- sum(combined)
  purity <- sum_max_per_cluster / total_predictions
  return(purity)
}
calculate_morans_i_proportions <- function(layer_manual_MOB, proportions) {
  # Ensure x and y columns are present
  if (!("x" %in% colnames(layer_manual_MOB)) || !("y" %in% colnames(layer_manual_MOB))) {
    stop("The dataset must contain 'x' and 'y' columns.")
  }
  
  # Extract coordinates
  coords <- layer_manual_MOB[, c("x", "y")]
  
  # Create neighbor object using knn
  nb <- knn2nb(knearneigh(coords, k = 4), sym = TRUE)
  
  # Convert neighbors to weights
  weights <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Create a dataframe to store Moran's I results
  moran_results <- data.frame(CellType = character(), MoranI = numeric(), p.value = numeric())
  
  # Loop through each proportion column
  for (i in 1:ncol(proportions)) {
    # Perform Moran's I test
    moran_result <- moran.test(proportions[, i], weights, zero.policy = TRUE)
    
    # Add results to the dataframe
    moran_results <- rbind(moran_results, data.frame(
      CellType = colnames(proportions)[i],
      MoranI = moran_result$estimate["Moran I statistic"],
      p.value = moran_result$p.value
    ))
  }
  return(moran_results)
}
plot_error_comparisons <- function(reference, result_list) {
  reference <- reference[sapply(reference, is.numeric)]  # Filter numeric columns
  
  if (is.null(names(result_list))) {
    stop("result_list must have names for each method")
  }
  
  error_metrics <- list()
  
  for (i in seq_along(result_list)) {
    method_name <- names(result_list)[i]
    result <- result_list[[i]][, sapply(result_list[[i]], is.numeric)]  # Filter numeric columns
    
    if (!all(dim(reference) == dim(result))) {
      stop(paste("Dimensions of reference and result do not match for method:", method_name))
    }
    
    # Calculate row-wise RMSE
    rowwise_rmse <- sqrt(rowMeans((reference - result)^2))
    overall_rmse <- mean(rowwise_rmse)
    
    # Calculate RMSE per cell type (column-wise)
    colwise_rmse <- sqrt(colMeans((reference - result)^2))
    
    # Store RMSE data
    error_metrics[[method_name]] <- list(
      overall_rmse = overall_rmse,
      rowwise_rmse = rowwise_rmse,
      colwise_rmse = colwise_rmse
    )
  }
  
  return(error_metrics)
}


run_analysis <- function(ground_truth, result_matrix, layer_manual_MOB) {
  aligned_data <- align_data_frames(ground_truth, result_matrix)
  reference <- aligned_data$reference
  result <- aligned_data$result
  
  # Calculate ARI and Purity
  ari_and_purity_results <- calculate_ari_and_purity(result, reference)
  
  # Calculate Moran's I for reference and result
  moran_results_reference <- calculate_morans_i_proportions(layer_manual_MOB, reference)
  moran_results_result <- calculate_morans_i_proportions(layer_manual_MOB, result)
  
  # Merge Moran's I results by cell type and calculate differences
  moran_differences <- merge(moran_results_reference, moran_results_result, by = "CellType", suffixes = c("_Ref", "_Res"))
  moran_differences$MoranI_Difference <- moran_differences$MoranI_Res - moran_differences$MoranI_Ref
  
  # Prepare the final structure for differences
  moran_diff_list <- setNames(as.list(moran_differences$MoranI_Difference), paste("MoransI", moran_differences$CellType, sep = "_"))
  
  error_results <- plot_error_comparisons(reference, list(result = result))
  
  # Plot dominant cell types
  plot <- plot_dominant_ct(result, layer_manual_MOB)
  
  return(list(
    ARI_and_Purity = ari_and_purity_results,
    Moran_Differences = moran_diff_list,
    error_results = error_results
  ))
}


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
  p <- ggplot(df_dominant, aes(x = x, y = y, color = DominantType)) +
    geom_point(size = point_size) + 
    common_theme +
    common_labs +
    scale_color_manual(values = c("lightblue",
                                  "orange",
                                  "red",
                                  "#9999FF")) 
  print(p)
}




####read results 
read_result <- function(name){
  
  proportion_matrix <- read.csv(here("Spatial_Blade", "Results_new", name))
  
  # Ensure rownames match the column names of true_labels_df
  rownames(proportion_matrix) <- rownames(layer_manual_MOB)
  
  return(proportion_matrix)
  
}

#####
create_combined_heatmap <- function(df, title_label) {
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(patchwork)
  
  # Define the groups
  group1_metrics <- df %>% filter(Metric %in% c("ARI", "Node_Purity"))
  group2_metrics <- df %>% filter(Metric == "Overall_RMSE")
  group3_metrics <- df %>% filter(Metric %in% c("MoransI_GC", "MoransI_PGC", "MoransI_OSNs", "MoransI_MTC"))
  
  # Melt the data frames
  group1_melted <- melt(group1_metrics, id.vars = c("Condition", "Metric"))
  group2_melted <- melt(group2_metrics, id.vars = c("Condition", "Metric"))
  group3_melted <- melt(group3_metrics, id.vars = c("Condition", "Metric"))
  
  # Create heatmaps for each group
  heatmap1 <- ggplot(group1_melted, aes(x = Condition, y = Metric, fill = value)) +
    geom_tile(color = "black") +  # Added border around tiles
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3, check_overlap = TRUE) +  # Added text labels
    scale_fill_gradientn(colors = c("white", "red"), 
                         values = rescale(c(0, 0.5, 1)),
                         limits = c(0, 1),
                         space = "Lab",   
                         guide = FALSE) +  # Removed color scale
    labs(title = NULL) +  # Removed title
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),  # Customize y-axis text
          axis.title.y = element_blank(),  # Removed y-axis title
          axis.text.x = element_blank(),  # Removed x-axis text
          axis.title.x = element_blank())  # Removed x-axis title
  
  
  heatmap2 <- ggplot(group2_melted, aes(x = Condition, y = Metric, fill = value)) +
    geom_tile(color = "black") +  # Added border around tiles
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3, check_overlap = TRUE) +  # Added text labels
    scale_fill_gradientn(colors = c("red","white"),
                         space = "Lab",   
                         guide = FALSE) +  # Removed color scale
    labs(title = NULL) +  # Removed title
    theme_minimal() +
    theme(axis.text.y = element_blank(),  # Removed y-axis text
          axis.title.y = element_blank(),  # Removed y-axis title
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  # Customize x-axis text
  
  heatmap3 <- ggplot(group3_melted, aes(x = Condition, y = Metric, fill = value)) +
    geom_tile(color = "black") +  # Added border around tiles
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3, check_overlap = TRUE) +  # Added text labels
    scale_fill_gradientn(colors = c("white","red"), 
                         space = "Lab",   
                         guide = FALSE) +  # Removed color scale
    labs(title = NULL) +  # Removed title
    theme_minimal() +
    theme(axis.text.y = element_blank(),  # Removed y-axis text
          axis.title.y = element_blank(),  # Removed y-axis title
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  # Customize x-axis text
  
  # Add annotation for metric set label
  heatmap1 <- heatmap1 + annotate("text", x = 0.5, y = 1.1, label = title_label, size = 6, hjust = 0.5)
  
  # Combine the heatmaps using patchwork
  combined_heatmap <- heatmap1 / heatmap2 / heatmap3  # Stack the heatmaps
  
  return(combined_heatmap)
}


blade_baseline1_ungroup <- read_result("BLADE_baseline1_independent.csv")
blade_baseline1 <- read_result("BLADE_baseline1_group.csv")
blade_baseline2 <- read_result("BLADE_baseline2_group.csv")
blade_baseline3 <- read_result("BLADE_baseline3_group.csv")
blade_hard <- read_result("BLADE_hard.csv")
blade_medium <- read_result("BLADE_medium.csv")




bladesp_baseline1 <- read_result("BLADE_baseline1sp.csv")
bladesp_baseline2 <- read_result("BLADE_baseline2sp.csv")
bladesp_baseline3 <- read_result("BLADE_baseline3sp.csv")
bladesp_hard <- read_result("BLADE_baseline_hardsp.csv")
bladesp_medium <- read_result("BLADE_baseline_mediumsp.csv")







run_analysis





card_baseline1 <- run_analysis(data.frame(p1), data.frame(card_baseline1), layer_manual_MOB)
card_baseline2 <- run_analysis(data.frame(p2), data.frame(card_baseline2), layer_manual_MOB)
card_baseline3 <- run_analysis(data.frame(p3), data.frame(card_baseline3), layer_manual_MOB)
card_medium <- run_analysis(data.frame(p4), data.frame(card_medium), layer_manual_MOB)
card_hard <- run_analysis(data.frame(p5), data.frame(card_hard), layer_manual_MOB)




music_baseline1 <- run_analysis(data.frame(p1), data.frame(music_baseline), layer_manual_MOB)
music_baseline2 <- run_analysis(data.frame(p2), data.frame(music_baseline2), layer_manual_MOB)
music_baseline3 <- run_analysis(data.frame(p3), data.frame(music_baseline3), layer_manual_MOB)
music_medium <- run_analysis(data.frame(p4), data.frame(music_medium), layer_manual_MOB)
music_hard <- run_analysis(data.frame(p5), data.frame(music_hard), layer_manual_MOB)






cibersort_baseline1 <- run_analysis(data.frame(p1), data.frame(cibersort_baseline1), layer_manual_MOB)
cibersort_baseline2 <- run_analysis(data.frame(p2), data.frame(cibersort_baseline2), layer_manual_MOB)
cibersort_baseline3 <- run_analysis(data.frame(p3), data.frame(cibersort_baseline3), layer_manual_MOB)
cibersort_medium <- run_analysis(data.frame(p4), data.frame(cibersort_medium), layer_manual_MOB)
cibersort_hard <- run_analysis(data.frame(p5), data.frame(cibersort_hard), layer_manual_MOB)






rctd_baseline1 <- run_analysis(data.frame(p1), data.frame(rctd_baseline1), layer_manual_MOB)
rctd_baseline2 <- run_analysis(data.frame(p2), data.frame(rctd_baseline2), layer_manual_MOB)
rctd_baseline3 <- run_analysis(data.frame(p3), data.frame(rctd_baseline3), layer_manual_MOB)
rctd_medium <- run_analysis(data.frame(p4), data.frame(rctd_medium), layer_manual_MOB)
rctd_hard <- run_analysis(data.frame(p5), data.frame(rctd_hard), layer_manual_MOB)




blade_baseline1 <- run_analysis(data.frame(p1), data.frame(blade_baseline1), layer_manual_MOB)

blade_baseline1_ungroup <- run_analysis(data.frame(p1), data.frame(blade_baseline1_ungroup), layer_manual_MOB)



blade_baseline2 <- run_analysis(data.frame(p2), data.frame(blade_baseline2), layer_manual_MOB)
blade_baseline3 <- run_analysis(data.frame(p3), data.frame(blade_baseline3), layer_manual_MOB)
blade_medium <- run_analysis(data.frame(p4), data.frame(blade_medium), layer_manual_MOB)
blade_hard <- run_analysis(data.frame(p5), data.frame(blade_hard), layer_manual_MOB)



bladesp_baseline1 <- run_analysis(data.frame(p1), data.frame(bladesp_baseline1), layer_manual_MOB)
bladesp_baseline2 <- run_analysis(data.frame(p2), data.frame(bladesp_baseline2), layer_manual_MOB)
bladesp_baseline3 <- run_analysis(data.frame(p3), data.frame(bladesp_baseline3), layer_manual_MOB)
bladesp_medium <- run_analysis(data.frame(p4), data.frame(bladesp_medium), layer_manual_MOB)
bladesp_hard <- run_analysis(data.frame(p5), data.frame(bladesp_hard), layer_manual_MOB)


blade_baseline2$error_results$result$rowwise_rmse


blade_metrics <- data.frame(
  Condition = rep(c('Baseline1', 'Baseline2', 'Baseline3', 'Medium', 'Hard'), times = 7),
  Metric = rep(c('ARI', 'Node_Purity', 'MoransI_GC', 'MoransI_MTC', 'MoransI_OSNs', 'MoransI_PGC', 'Overall_RMSE'), each = 5),
  Value = c(
    blade_baseline1$ARI_and_Purity$ARI, blade_baseline2$ARI_and_Purity$ARI, blade_baseline3$ARI_and_Purity$ARI, blade_medium$ARI_and_Purity$ARI, blade_hard$ARI_and_Purity$ARI,
    blade_baseline1$ARI_and_Purity$Node_Purity, blade_baseline2$ARI_and_Purity$Node_Purity, blade_baseline3$ARI_and_Purity$Node_Purity, blade_medium$ARI_and_Purity$Node_Purity, blade_hard$ARI_and_Purity$Node_Purity,
    blade_baseline1$Moran_Differences$MoransI_GC, blade_baseline2$Moran_Differences$MoransI_GC, blade_baseline3$Moran_Differences$MoransI_GC, blade_medium$Moran_Differences$MoransI_GC, blade_hard$Moran_Differences$MoransI_GC,
    blade_baseline1$Moran_Differences$MoransI_M.TC, blade_baseline2$Moran_Differences$MoransI_M.TC, blade_baseline3$Moran_Differences$MoransI_M.TC, blade_medium$Moran_Differences$MoransI_M.TC, blade_hard$Moran_Differences$MoransI_M.TC,
    blade_baseline1$Moran_Differences$MoransI_OSNs, blade_baseline2$Moran_Differences$MoransI_OSNs, blade_baseline3$Moran_Differences$MoransI_OSNs, blade_medium$Moran_Differences$MoransI_OSNs, blade_hard$Moran_Differences$MoransI_OSNs,
    blade_baseline1$Moran_Differences$MoransI_PGC, blade_baseline2$Moran_Differences$MoransI_PGC, blade_baseline3$Moran_Differences$MoransI_PGC, blade_medium$Moran_Differences$MoransI_PGC, blade_hard$Moran_Differences$MoransI_PGC,
    blade_baseline1$error_results$result$overall_rmse, blade_baseline2$error_results$result$overall_rmse, blade_baseline3$error_results$result$overall_rmse, blade_medium$error_results$result$overall_rmse
    , blade_hard$error_results$result$overall_rmse
  )
)






bladesp_metrics <- data.frame(
  Condition = rep(c('Baseline1', 'Baseline2', 'Baseline3', 'Medium', 'Hard'), times = 7),
  Metric = rep(c('ARI', 'Node_Purity', 'MoransI_GC', 'MoransI_MTC', 'MoransI_OSNs', 'MoransI_PGC', 'Overall_RMSE'), each = 5),
  Value = c(
    bladesp_baseline1$ARI_and_Purity$ARI, bladesp_baseline2$ARI_and_Purity$ARI, bladesp_baseline3$ARI_and_Purity$ARI, bladesp_medium$ARI_and_Purity$ARI, bladesp_hard$ARI_and_Purity$ARI,
    bladesp_baseline1$ARI_and_Purity$Node_Purity, bladesp_baseline2$ARI_and_Purity$Node_Purity, bladesp_baseline3$ARI_and_Purity$Node_Purity, bladesp_medium$ARI_and_Purity$Node_Purity, bladesp_hard$ARI_and_Purity$Node_Purity,
    bladesp_baseline1$Moran_Differences$MoransI_GC, bladesp_baseline2$Moran_Differences$MoransI_GC, bladesp_baseline3$Moran_Differences$MoransI_GC, bladesp_medium$Moran_Differences$MoransI_GC, bladesp_hard$Moran_Differences$MoransI_GC,
    bladesp_baseline1$Moran_Differences$MoransI_M.TC, bladesp_baseline2$Moran_Differences$MoransI_M.TC, bladesp_baseline3$Moran_Differences$MoransI_M.TC, bladesp_medium$Moran_Differences$MoransI_M.TC, bladesp_hard$Moran_Differences$MoransI_M.TC,
    bladesp_baseline1$Moran_Differences$MoransI_OSNs, bladesp_baseline2$Moran_Differences$MoransI_OSNs, bladesp_baseline3$Moran_Differences$MoransI_OSNs, bladesp_medium$Moran_Differences$MoransI_OSNs, bladesp_hard$Moran_Differences$MoransI_OSNs,
    bladesp_baseline1$Moran_Differences$MoransI_PGC, bladesp_baseline2$Moran_Differences$MoransI_PGC, bladesp_baseline3$Moran_Differences$MoransI_PGC, bladesp_medium$Moran_Differences$MoransI_PGC, bladesp_hard$Moran_Differences$MoransI_PGC,
    bladesp_baseline1$error_results$result$overall_rmse, bladesp_baseline2$error_results$result$overall_rmse, bladesp_baseline3$error_results$result$overall_rmse, bladesp_medium$error_results$result$overall_rmse
    , bladesp_hard$error_results$result$overall_rmse
  )
)



card_metrics <- data.frame(
  Condition = rep(c('Baseline1', 'Baseline2', 'Baseline3', 'Medium', 'Hard'), times = 7),
  Metric = rep(c('ARI', 'Node_Purity', 'MoransI_GC', 'MoransI_MTC', 'MoransI_OSNs', 'MoransI_PGC', 'Overall_RMSE'), each = 5),
  Value = c(
    card_baseline1$ARI_and_Purity$ARI, card_baseline2$ARI_and_Purity$ARI, card_baseline3$ARI_and_Purity$ARI, card_medium$ARI_and_Purity$ARI, card_hard$ARI_and_Purity$ARI,
    card_baseline1$ARI_and_Purity$Node_Purity, card_baseline2$ARI_and_Purity$Node_Purity, card_baseline3$ARI_and_Purity$Node_Purity, card_medium$ARI_and_Purity$Node_Purity, card_hard$ARI_and_Purity$Node_Purity,
    card_baseline1$Moran_Differences$MoransI_GC, card_baseline2$Moran_Differences$MoransI_GC, card_baseline3$Moran_Differences$MoransI_GC, card_medium$Moran_Differences$MoransI_GC, card_hard$Moran_Differences$MoransI_GC,
    card_baseline1$Moran_Differences$MoransI_M.TC, card_baseline2$Moran_Differences$MoransI_M.TC, card_baseline3$Moran_Differences$MoransI_M.TC, card_medium$Moran_Differences$MoransI_M.TC, card_hard$Moran_Differences$MoransI_M.TC,
    card_baseline1$Moran_Differences$MoransI_OSNs, card_baseline2$Moran_Differences$MoransI_OSNs, card_baseline3$Moran_Differences$MoransI_OSNs, card_medium$Moran_Differences$MoransI_OSNs, card_hard$Moran_Differences$MoransI_OSNs,
    card_baseline1$Moran_Differences$MoransI_PGC, card_baseline2$Moran_Differences$MoransI_PGC, card_baseline3$Moran_Differences$MoransI_PGC, card_medium$Moran_Differences$MoransI_PGC, card_hard$Moran_Differences$MoransI_PGC,
    card_baseline1$error_results$result$overall_rmse, card_baseline2$error_results$result$overall_rmse, card_baseline3$error_results$result$overall_rmse, card_medium$error_results$result$overall_rmse
    , card_hard$error_results$result$overall_rmse
  )
)

music_metrics <- data.frame(
  Condition = rep(c('Baseline1', 'Baseline2', 'Baseline3', 'Medium', 'Hard'), times = 7),
  Metric = rep(c('ARI', 'Node_Purity', 'MoransI_GC', 'MoransI_MTC', 'MoransI_OSNs', 'MoransI_PGC', 'Overall_RMSE'), each = 5),
  Value = c(
    music_baseline1$ARI_and_Purity$ARI, music_baseline2$ARI_and_Purity$ARI, music_baseline3$ARI_and_Purity$ARI, music_medium$ARI_and_Purity$ARI, music_hard$ARI_and_Purity$ARI,
    music_baseline1$ARI_and_Purity$Node_Purity, music_baseline2$ARI_and_Purity$Node_Purity, music_baseline3$ARI_and_Purity$Node_Purity, music_medium$ARI_and_Purity$Node_Purity, music_hard$ARI_and_Purity$Node_Purity,
    music_baseline1$Moran_Differences$MoransI_GC, music_baseline2$Moran_Differences$MoransI_GC, music_baseline3$Moran_Differences$MoransI_GC, music_medium$Moran_Differences$MoransI_GC, music_hard$Moran_Differences$MoransI_GC,
    music_baseline1$Moran_Differences$MoransI_M.TC, music_baseline2$Moran_Differences$MoransI_M.TC, music_baseline3$Moran_Differences$MoransI_M.TC, music_medium$Moran_Differences$MoransI_M.TC, music_hard$Moran_Differences$MoransI_M.TC,
    music_baseline1$Moran_Differences$MoransI_OSNs, music_baseline2$Moran_Differences$MoransI_OSNs, music_baseline3$Moran_Differences$MoransI_OSNs, music_medium$Moran_Differences$MoransI_OSNs, music_hard$Moran_Differences$MoransI_OSNs,
    music_baseline1$Moran_Differences$MoransI_PGC, music_baseline2$Moran_Differences$MoransI_PGC, music_baseline3$Moran_Differences$MoransI_PGC, music_medium$Moran_Differences$MoransI_PGC, music_hard$Moran_Differences$MoransI_PGC,
    music_baseline1$error_results$result$overall_rmse, music_baseline2$error_results$result$overall_rmse, music_baseline3$error_results$result$overall_rmse, music_medium$error_results$result$overall_rmse
    , music_hard$error_results$result$overall_rmse
  )
)


cibersort_metrics <- data.frame(
  Condition = rep(c('Baseline1', 'Baseline2', 'Baseline3', 'Medium', 'Hard'), times = 7),
  Metric = rep(c('ARI', 'Node_Purity', 'MoransI_GC', 'MoransI_MTC', 'MoransI_OSNs', 'MoransI_PGC', 'Overall_RMSE'), each = 5),
  Value = c(
    cibersort_baseline1$ARI_and_Purity$ARI, cibersort_baseline2$ARI_and_Purity$ARI, cibersort_baseline3$ARI_and_Purity$ARI, cibersort_medium$ARI_and_Purity$ARI, cibersort_hard$ARI_and_Purity$ARI,
    cibersort_baseline1$ARI_and_Purity$Node_Purity, cibersort_baseline2$ARI_and_Purity$Node_Purity, cibersort_baseline3$ARI_and_Purity$Node_Purity, cibersort_medium$ARI_and_Purity$Node_Purity, cibersort_hard$ARI_and_Purity$Node_Purity,
    cibersort_baseline1$Moran_Differences$MoransI_GC, cibersort_baseline2$Moran_Differences$MoransI_GC, cibersort_baseline3$Moran_Differences$MoransI_GC, cibersort_medium$Moran_Differences$MoransI_GC, cibersort_hard$Moran_Differences$MoransI_GC,
    cibersort_baseline1$Moran_Differences$MoransI_M.TC, cibersort_baseline2$Moran_Differences$MoransI_M.TC, cibersort_baseline3$Moran_Differences$MoransI_M.TC, cibersort_medium$Moran_Differences$MoransI_M.TC, cibersort_hard$Moran_Differences$MoransI_M.TC,
    cibersort_baseline1$Moran_Differences$MoransI_OSNs, cibersort_baseline2$Moran_Differences$MoransI_OSNs, cibersort_baseline3$Moran_Differences$MoransI_OSNs, cibersort_medium$Moran_Differences$MoransI_OSNs, cibersort_hard$Moran_Differences$MoransI_OSNs,
    cibersort_baseline1$Moran_Differences$MoransI_PGC, cibersort_baseline2$Moran_Differences$MoransI_PGC, cibersort_baseline3$Moran_Differences$MoransI_PGC, cibersort_medium$Moran_Differences$MoransI_PGC, cibersort_hard$Moran_Differences$MoransI_PGC,
    cibersort_baseline1$error_results$result$overall_rmse, cibersort_baseline2$error_results$result$overall_rmse, cibersort_baseline3$error_results$result$overall_rmse, cibersort_medium$error_results$result$overall_rmse
    , cibersort_hard$error_results$result$overall_rmse
  )
)


rctd_metrics <- data.frame(
  Condition = rep(c('Baseline1', 'Baseline2', 'Baseline3', 'Medium', 'Hard'), times = 7),
  Metric = rep(c('ARI', 'Node_Purity', 'MoransI_GC', 'MoransI_MTC', 'MoransI_OSNs', 'MoransI_PGC', 'Overall_RMSE'), each = 5),
  Value = c(
    rctd_baseline1$ARI_and_Purity$ARI, rctd_baseline2$ARI_and_Purity$ARI, rctd_baseline3$ARI_and_Purity$ARI, rctd_medium$ARI_and_Purity$ARI, rctd_hard$ARI_and_Purity$ARI,
    rctd_baseline1$ARI_and_Purity$Node_Purity, rctd_baseline2$ARI_and_Purity$Node_Purity, rctd_baseline3$ARI_and_Purity$Node_Purity, rctd_medium$ARI_and_Purity$Node_Purity, rctd_hard$ARI_and_Purity$Node_Purity,
    rctd_baseline1$Moran_Differences$MoransI_GC, rctd_baseline2$Moran_Differences$MoransI_GC, rctd_baseline3$Moran_Differences$MoransI_GC, rctd_medium$Moran_Differences$MoransI_GC, rctd_hard$Moran_Differences$MoransI_GC,
    rctd_baseline1$Moran_Differences$MoransI_M.TC, rctd_baseline2$Moran_Differences$MoransI_M.TC, rctd_baseline3$Moran_Differences$MoransI_M.TC, rctd_medium$Moran_Differences$MoransI_M.TC, rctd_hard$Moran_Differences$MoransI_M.TC,
    rctd_baseline1$Moran_Differences$MoransI_OSNs, rctd_baseline2$Moran_Differences$MoransI_OSNs, rctd_baseline3$Moran_Differences$MoransI_OSNs, rctd_medium$Moran_Differences$MoransI_OSNs, rctd_hard$Moran_Differences$MoransI_OSNs,
    rctd_baseline1$Moran_Differences$MoransI_PGC, rctd_baseline2$Moran_Differences$MoransI_PGC, rctd_baseline3$Moran_Differences$MoransI_PGC, rctd_medium$Moran_Differences$MoransI_PGC, rctd_hard$Moran_Differences$MoransI_PGC,
    rctd_baseline1$error_results$result$overall_rmse, rctd_baseline2$error_results$result$overall_rmse, rctd_baseline3$error_results$result$overall_rmse, rctd_medium$error_results$result$overall_rmse
    , rctd_hard$error_results$result$overall_rmse
  )
)




blade_baseline1$Moran_Differences$MoransI_M.TC
create_combined_heatmap(card_metrics, "Card Metrics") |
  create_combined_heatmap(music_metrics, "Music Metrics") |
  create_combined_heatmap(rctd_metrics, "RCTD Metrics") |
  create_combined_heatmap(cibersort_metrics, "Cibersort Metrics") |
  create_combined_heatmap(blade_metrics, "Blade Metrics")



create_combined_heatmap <- function(df, show_x_labels = FALSE, show_scale = FALSE, show_y_labels = TRUE, name) {
  # Define the groups
  group1_metrics <- df %>% filter(Metric %in% c("ARI", "Node_Purity"))
  group2_metrics <- df %>% filter(Metric == "Overall_RMSE")
  group3_metrics <- df %>% filter(Metric %in% c("MoransI_GC", "MoransI_PGC", "MoransI_OSNs", "MoransI_M.TC"))
  
  # Melt the data frames
  group1_melted <- melt(group1_metrics, id.vars = c("Condition", "Metric"))
  group2_melted <- melt(group2_metrics, id.vars = c("Condition", "Metric"))
  group3_melted <- melt(group3_metrics, id.vars = c("Condition", "Metric"))
  
  # Ensure the 'value' column is numeric
  group1_melted$value <- as.numeric(group1_melted$value)
  group2_melted$value <- as.numeric(group2_melted$value)
  group3_melted$value <- as.numeric(group3_melted$value)
  
  # Create heatmaps for each group
  heatmap1 <- ggplot(group1_melted, aes(x = Condition, y = Metric, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3, check_overlap = TRUE) +
    scale_fill_gradientn(colors = c("white", "red"), 
                         values = rescale(c(0, 0.5, 1)),
                         limits = c(0, 1),
                         space = "Lab",   
                         guide = if (show_scale) "colourbar" else "none") + 
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = if (show_y_labels) element_text() else element_blank()) + ggtitle(name)
  
  heatmap3 <- ggplot(group2_melted, aes(x = Condition, y = Metric, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3, check_overlap = TRUE) +
    scale_fill_gradientn(colors = c("red", "white"),
                         space = "Lab",   
                         guide = if (show_scale) "colourbar" else "none") +
    theme_minimal() +
    theme(axis.text.x = if (show_x_labels) element_text(angle = 45, hjust = 1) else element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = if (show_y_labels) element_text() else element_blank())
  
  heatmap2 <- ggplot(group3_melted, aes(x = Condition, y = Metric, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3, check_overlap = TRUE) +
    scale_fill_gradientn(colors = c("white", "red"), 
                         space = "Lab",   
                         guide = if (show_scale) "colourbar" else "none") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = if (show_y_labels) element_text() else element_blank())
  
  # Combine the heatmaps using patchwork
  combined_heatmap <- heatmap1  / heatmap2 / heatmap3  # Stack the heatmaps
  
  return(combined_heatmap)
}

# Combine all methods and apply the custom settings

create_combined_heatmap(card_metrics, show_x_labels = TRUE, show_scale = FALSE, show_y_labels = TRUE, "CARD") |
  create_combined_heatmap(music_metrics, show_x_labels = TRUE, show_scale = FALSE, show_y_labels = FALSE, "MUSIC") |
  create_combined_heatmap(rctd_metrics, show_x_labels = TRUE, show_scale = FALSE, show_y_labels = FALSE, "RCTD") |
  create_combined_heatmap(cibersort_metrics, show_x_labels = TRUE, show_scale = FALSE, show_y_labels = FALSE, "CIBSERSORTx") |
  create_combined_heatmap(blade_metrics, show_x_labels = TRUE, show_scale = TRUE, show_y_labels = FALSE, "BLADE") 


blade_comparison_metrics <- data.frame(
  Condition = rep(c('Independent', 'Group'), times = 7),
  Metric = rep(c('ARI', 'Node_Purity', 'MoransI_GC', 'MoransI_MTC', 'MoransI_OSNs', 'MoransI_PGC', 'Overall_RMSE'), each = 2),
  Value = c(
    blade_baseline1_ungroup$ARI_and_Purity$ARI, blade_baseline1$ARI_and_Purity$ARI,
    blade_baseline1_ungroup$ARI_and_Purity$Node_Purity, blade_baseline1$ARI_and_Purity$Node_Purity,
    blade_baseline1_ungroup$Moran_Differences$MoransI_GC, blade_baseline1$Moran_Differences$MoransI_GC,
    blade_baseline1_ungroup$Moran_Differences$MoransI_M.TC, blade_baseline1$Moran_Differences$MoransI_M.TC,
    blade_baseline1_ungroup$Moran_Differences$MoransI_OSNs, blade_baseline1$Moran_Differences$MoransI_OSNs,
    blade_baseline1_ungroup$Moran_Differences$MoransI_PGC, blade_baseline1$Moran_Differences$MoransI_PGC,
    blade_baseline1_ungroup$error_results$result$overall_rmse, blade_baseline1$error_results$result$overall_rmse
  )
)


create_combined_heatmap(blade_comparison_metrics, show_x_labels = TRUE, show_scale = TRUE, show_y_labels = TRUE, "BLADE ALPHA GROUP vs INDEPENDENT") 

print(blade_baseline1_ungroup$error_results$result$overall_rmse)
print(blade_baseline1$error_results$result$overall_rmse)










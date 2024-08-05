######Spatial Data Simulator#####

###Aryamaan Bose


library(vegan)






###Data extracted from layer_manual_MOB
##Its a manual annotation of tissue which is He stained 
#the same tissue is used for the spatial transcriptomics experiment 


####Check which layers are more correllated 








###SPOTS 
##The spots are the rownames of the spatial experiment (spe) and the manual annotation layer_manual_MOB
##It is stored as Ycoordinate x Ycoordinate
org1 <- read.csv("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/org1.csv")
colnames(org1)[which(colnames(org1) == "cluster")] <- "celltypes"
org1$Spots <- paste0("spot_", seq_len(nrow(org1)))


##
load("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/Figure3A_layer_annote (1).RData")

layer_to_cell_type <- list(
  GCL = "GC",
  GL = "PGC",
  ONL = "OSNs",
  MCL = "M/TC"
)

layer_manual_MOB$Spot <- rownames(layer_manual_MOB)
layer_manual_MOB$celltypes <- sapply(layer_manual_MOB$Layer, function(x) layer_to_cell_type[[x]])

######
compute_proportions6 <- function(layer_manual_MOB, pn = 0.0, sc = 0.0, specified_dominant_proportions = NULL, neighborhood_radius, add_labels = FALSE, adjust_zeroes = FALSE) {
  set.seed(13)
  
  # Extract the spot identifiers directly from the row names of layer_manual_MOB
  spots <- rownames(layer_manual_MOB)
  n_spots <- length(spots)
  
  # Extract unique cell types from the 'celltypes' column of layer_manual_MOB
  unique_cell_types <- unique(layer_manual_MOB$celltypes)
  
  # Initialize a matrix to hold the proportions
  proportions <- matrix(0, nrow = n_spots, ncol = length(unique_cell_types))
  colnames(proportions) <- unique_cell_types
  rownames(proportions) <- spots
  
  # Assign dominant proportions (either specified or using a default value)
  if (is.null(specified_dominant_proportions)) {
    avg_proportions <- rep(0.75, length(unique_cell_types))
    names(avg_proportions) <- unique_cell_types
  } else {
    avg_proportions <- specified_dominant_proportions
    names(avg_proportions) <- unique_cell_types
  }
  
  # Extract coordinates for all spots to compute a distance matrix
  coords <- layer_manual_MOB[, c("x", "y")]
  dist_matrix <- as.matrix(dist(coords))
  
  # Initialize label columns if requested
  if (add_labels) {
    labels <- data.frame(Dominant = character(n_spots), Secondary = character(n_spots), Tertiary = character(n_spots), stringsAsFactors = FALSE)
  }
  
  # Compute the preliminary dominant proportions
  for (i in 1:n_spots) {
    spot <- spots[i]
    dominant_cell_type <- layer_manual_MOB$celltypes[which(rownames(layer_manual_MOB) == spot)]
    
    # Change the dominant cell type based on `pn` probability
    if (runif(1) < pn) {
      dominant_cell_type <- sample(setdiff(unique_cell_types, dominant_cell_type), 1)
    }
    
    # Apply the dominant proportion directly
    dominant_proportion <- avg_proportions[[dominant_cell_type]] + runif(1, -sc, sc)  # Add noise if `sc` isn't zero
    dominant_proportion <- max(min(dominant_proportion, 0.9), 0)  # Ensure it's within bounds
    proportions[i, dominant_cell_type] <- dominant_proportion
    
    # Identify neighboring spots within a given radius
    neighbors <- which(dist_matrix[i, ] <= neighborhood_radius & dist_matrix[i, ] > 0)
    remaining_percentage <- 1 - dominant_proportion
    
    if (length(neighbors) > 0) {
      # Dynamically allocate remaining proportion to neighbors based on their count
      equal_split <- remaining_percentage / length(neighbors)
      neighbor_cell_types <- layer_manual_MOB$celltypes[neighbors]
      
      for (type in unique(neighbor_cell_types)) {
        proportions[i, type] <- proportions[i, type] + equal_split * sum(neighbor_cell_types == type)
      }
    }
    
    # Assign labels if required
    if (add_labels) {
      labels$Dominant[i] <- dominant_cell_type
      labels$Secondary[i] <- NA
      labels$Tertiary[i] <- NA
      if (length(neighbors) > 0) {
        sorted_types <- names(sort(table(neighbor_cell_types), decreasing = TRUE))
        if (length(sorted_types) >= 1) labels$Secondary[i] <- sorted_types[1]
        if (length(sorted_types) >= 2) labels$Tertiary[i] <- sorted_types[2]
      }
    }
  }
  
  # Normalize each row to sum to 1
  proportions <- proportions / rowSums(proportions)
  
  # Adjust rows with zeroes if requested
  if (adjust_zeroes) {
    proportions[proportions == 0] <- 0.01
    proportions <- t(apply(proportions, 1, function(row) row / sum(row)))
  }
  
  if (add_labels) {
    return(cbind(as.data.frame(proportions), labels))
  } else {
    return(as.data.frame(proportions))
  }
}




###Finding optimal radii based on n nearest nieghbours 
calculate_optimal_radius <- function(coords, k=270) {
  # Calculate the distance matrix
  dist_matrix <- as.matrix(dist(coords))
  
  # Initialize a vector to store optimal radii
  optimal_radii <- numeric(nrow(coords))
  
  # Loop over each spot to compute the average distance to the k nearest neighbors
  for (i in 1:nrow(coords)) {
    # Exclude the current spot itself
    distances <- sort(dist_matrix[i, -i])
    
    # Get the average distance to the k-nearest neighbors
    optimal_radii[i] <- mean(distances[1:k])
  }
  
  return(optimal_radii)
}

coords <- layer_manual_MOB[, c("x", "y")]

table(layer_manual_MOB$celltypes)
optimal_radii <- calculate_optimal_radius(coords)

# To analyze diversity and uniformity, you can use these radii values
summary(optimal_radii)



######Checking shannon diversity for checking cell type diversity per spot when neighbourhood radius is increased 

calculate_diversity_per_radius <- function(layer_manual_MOB, radii) {
  # Extract coordinates
  coords <- layer_manual_MOB[, c("x", "y")]
  dist_matrix <- as.matrix(dist(coords))
  
  # List to store diversity values for each radius
  diversity_values <- data.frame(radius = numeric(), shannon_diversity = numeric())
  
  # Loop over each radius
  for (radius in radii) {
    # Initialize an empty vector to store diversity indices per spot
    spot_diversity <- numeric(nrow(coords))
    
    # Calculate diversity for each spot based on the radius
    for (i in 1:nrow(coords)) {
      # Find neighboring spots within the radius
      neighbors <- which(dist_matrix[i, ] <= radius & dist_matrix[i, ] > 0)
      
      if (length(neighbors) > 0) {
        # Extract cell types for neighbors
        neighbor_types <- layer_manual_MOB$celltypes[neighbors]
        type_counts <- table(neighbor_types)
        
        # Calculate Shannon Diversity Index for this neighborhood
        spot_diversity[i] <- diversity(type_counts, index = "shannon")
      } else {
        # No neighbors means zero diversity
        spot_diversity[i] = 0
      }
    }
    
    # Append the average diversity value for this radius
    diversity_values <- rbind(
      diversity_values,
      data.frame(radius = radius, shannon_diversity = mean(spot_diversity, na.rm = TRUE))
    )
  }
  
  return(diversity_values)
}

# Example usage:
radii <- c(1, 2, 3, 4, 5, 6,7,8,9,10,12,13)  # Adjust as needed
diversity_results <- calculate_diversity_per_radius(layer_manual_MOB, radii)
# Plot the results
ggplot(diversity_results, aes(x = radius, y = shannon_diversity)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Shannon Diversity Index vs. Neighborhood Radius",
    x = "Neighborhood Radius",
    y = "Shannon Diversity Index"
  ) +
  theme_minimal()


layer_manual_MOB




##Specify Dominant celltype prportion prportion 
##Initial Proportions if using noise parameter sc 
##If not specified will simulate between 30 to 40% for each cell type 
##Can be changed in the compute_proportions() if there are more celltypes
specified_props <- c(0.70, 0.70,  0.70, 0.70)
names(specified_props) <- c("GC", "M/TC", "OSNs", "PGC")


specified_props2 <- c(0.50, 0.50,  0.50, 0.50)
names(specified_props2) <- c("GC", "M/TC", "OSNs", "PGC")


specified_props3 <- c(0.30, 0.30,  0.30, 0.30)
names(specified_props3) <- c("GC", "M/TC", "OSNs", "PGC")


p1 <- compute_proportions6(layer_manual_MOB, pn = 0.00, sc = 0.00, specified_dominant_proportions = specified_props, neighborhood_radius = 2, add_labels = FALSE)


p2 <- compute_proportions6(layer_manual_MOB, pn = 0.0, sc = 0.00, specified_dominant_proportions = specified_props, neighborhood_radius = 3.8, add_labels = FALSE)


p3 <- compute_proportions6(layer_manual_MOB, pn = 0.0, sc = 0.00, specified_dominant_proportions = specified_props, neighborhood_radius = 7, add_labels = FALSE, adjust_zeroes = TRUE)



p4 <- compute_proportions6(layer_manual_MOB, pn = 0.0, sc = 0.00, specified_dominant_proportions = specified_props2, neighborhood_radius = 7, add_labels = FALSE, adjust_zeroes = TRUE)



p5 <- compute_proportions6(layer_manual_MOB, pn = 0.0, sc = 0.00, specified_dominant_proportions = specified_props3, neighborhood_radius = 7, add_labels = FALSE, adjust_zeroes = TRUE)




####VISUALIZE THE PROPORTIONS 
##Specify colours for the cell types 
my_colors <- c( "lightblue", "orange", "red", "#ed78ff")

##Visulize using CARDs built in function 

p1_plot <- CARD.visualize.pie(
  proportion = p1,
  spatial_location = layer_manual_MOB, 
  colors = my_colors, 
  radius = 0.5)

p2_plot <- CARD.visualize.pie(
  proportion = p2,
  spatial_location = layer_manual_MOB, 
  colors = my_colors, 
  radius = 0.5)

p3_plot <- CARD.visualize.pie(
  proportion = p3,
  spatial_location = layer_manual_MOB, 
  colors = my_colors, 
  radius = 0.5)

p4_plot <- CARD.visualize.pie(
  proportion = p4,
  spatial_location = layer_manual_MOB, 
  colors = my_colors, 
  radius = 0.5)

p5_plot <- CARD.visualize.pie(
  proportion = p5,
  spatial_location = layer_manual_MOB, 
  colors = my_colors, 
  radius = 0.5)

# Combine the plots into a grid
combined_plot <- plot_grid(
  p1_plot + theme(legend.position = "none"), 
  p2_plot + theme(legend.position = "none"), 
  p3_plot + theme(legend.position = "none"), 
  p4_plot + theme(legend.position = "none"), 
  p5_plot +theme(legend.position = "none"),
  labels = c("A", "B", "C", "D", "E"),
  ncol = 2
)

# Extract the legend from one of the plots


# Add the legend to the bottom of the combined plot
final_plot <- plot_grid(combined_plot , ncol = 1, rel_heights = c(1, 0.1))

# Display the final plot
print(final_plot)


###Visualize dominant cell types 




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
                                  "#c978ff"))  # Adjust colors as needed
}


####Visualize dominant cell types and noise
plot_dominant_ct(p5, layer_manual_MOB, point_size = 5)



plot_comparison <- function(proportions_sim, proportions_ref, layer_manual_MOB, point_size = 3) {
  # Helper function to calculate dominant and secondary indices
  get_types <- function(proportions) {
    dominant_idx <- apply(proportions, 1, function(x) order(x, decreasing = TRUE)[1])
    secondary_idx <- apply(proportions, 1, function(x) order(x, decreasing = TRUE)[2])
    cell_type_names <- colnames(proportions)
    
    list(
      Dominant = cell_type_names[dominant_idx],
      Secondary = cell_type_names[secondary_idx],
      Count = rowSums(proportions > 0)
    )
  }
  
  # Calculate types for both datasets
  types_sim <- get_types(proportions_sim)
  types_ref <- get_types(proportions_ref)
  
  # Prepare data frames
  df_sim <- data.frame(Spot = rownames(proportions_sim), DominantType = types_sim$Dominant,
                       SecondaryType = types_sim$Secondary, Count = types_sim$Count, 
                       Set = "Simulated", stringsAsFactors = FALSE)
  df_ref <- data.frame(Spot = rownames(proportions_ref), DominantType = types_ref$Dominant,
                       SecondaryType = types_ref$Secondary, Count = types_ref$Count, 
                       Set = "Reference", stringsAsFactors = FALSE)
  
  df_all <- rbind(df_sim, df_ref)
  
  # Ensure coordinates are present
  if (!all(c("x", "y") %in% colnames(layer_manual_MOB))) {
    stop("The input data frame 'layer_manual_MOB' must contain 'x' and 'y' columns.")
  }
  coordinates <- layer_manual_MOB[, c("x", "y")]
  coordinates$Spot <- rownames(layer_manual_MOB)
  
  df_all <- merge(df_all, coordinates, by = "Spot")
  
  # Calculate ARI and Node Purity for Dominant and Secondary
  ARI_dominant <- adjustedRandIndex(df_ref$DominantType, df_sim$DominantType)
  ARI_secondary <- adjustedRandIndex(df_ref$SecondaryType, df_sim$SecondaryType)
  
  purity_dominant <- calculateNodePurity(df_ref$DominantType, df_sim$DominantType)
  purity_secondary <- calculateNodePurity(df_ref$SecondaryType, df_sim$SecondaryType)
  
  # Generate plots
  plots <- list(
    Dominant = ggplot(df_all, aes(x = x, y = y, color = DominantType)) +
      geom_point(size = point_size) +
      facet_wrap(~Set, ncol = 1) +
      labs(title = paste("Comparison of Dominant Cell Types - ARI:", round(ARI_dominant, 2), "- Purity:", round(purity_dominant, 2))) +
      theme_minimal(),
    
    Secondary = ggplot(df_all, aes(x = x, y = y, color = SecondaryType)) +
      geom_point(size = point_size) +
      facet_wrap(~Set, ncol = 1) +
      labs(title = paste("Comparison of Secondary Cell Types - ARI:", round(ARI_secondary, 2), "- Purity:", round(purity_secondary, 2))) +
      theme_minimal()
  )
  
  return(plots)
}

calculateNodePurity <- function(true_labels, predicted_labels) {
  combined <- table(true_labels, predicted_labels)
  sum_max_per_cluster <- sum(apply(combined, 2, max))
  total_predictions <- sum(combined)
  purity <- sum_max_per_cluster / total_predictions
  return(purity)
}

# Example Usage:
result <- plot_comparison(card1, p1 , org1)
result$Dominant  # To view the plot for dominant cell types
result$Secondary  # To view the plot for secondary cell types
result$Count  # To view the plot for the count of cell types per spot














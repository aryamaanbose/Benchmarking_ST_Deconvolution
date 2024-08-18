



org1 <- read.csv("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Code/Simulation/g1.csv")
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



compute_proportions6 <- function(layer_manual_MOB, pn = 0.0, sc = 0.0, specified_dominant_proportions, neighborhood_radius, add_labels = FALSE, adjust_zeroes = FALSE) {
  
  
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
    dominant_proportion <- avg_proportions[[dominant_cell_type]] #+ runif(1, -sc, sc)  # Add noise if `sc` isn't zero
    #dominant_proportion <- max(min(dominant_proportion, 0.9), 0)  # Ensure it's within bounds
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
  #proportions <- proportions / rowSums(proportions)
  
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
compute_proportions7 <- function(layer_manual_MOB, pn = 0.0, sc = 0.0, specified_dominant_proportions, neighborhood_radius, add_labels = FALSE, adjust_zeroes = FALSE) {
  
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
    dominant_proportion <- avg_proportions[[dominant_cell_type]] # + runif(1, -sc, sc)  # Add noise if `sc` isn't zero
    #dominant_proportion <- max(min(dominant_proportion, 0.9), 0)  # Ensure it's within bounds
    proportions[i, dominant_cell_type] <- dominant_proportion
    
    # Identify neighboring spots within a given radius
    neighbors <- which(dist_matrix[i, ] <= neighborhood_radius & dist_matrix[i, ] > 0)
    remaining_percentage <- 1 - dominant_proportion
    
    if (length(neighbors) > 0) {
      # Dynamically allocate remaining proportion to neighbors based on their count
      equal_split <- remaining_percentage / length(neighbors)
      neighbor_cell_types <- layer_manual_MOB$celltypes[neighbors]
      
      for (type in unique(neighbor_cell_types)) {
        if (type != dominant_cell_type) {
          proportions[i, type] <- proportions[i, type] + equal_split * sum(neighbor_cell_types == type)
        }
      }
    }
    
    # Recalculate the proportions of non-dominant cell types
    other_cell_types <- setdiff(unique_cell_types, dominant_cell_type)
    proportions[i, other_cell_types] <- proportions[i, other_cell_types] * (1 - dominant_proportion) / sum(proportions[i, other_cell_types])
    
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
  
  # Optionally adjust rows with zeroes
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
specified_props <- c(0.90, 0.90)

specified_props2 <- c(1, 1)


names(specified_props) <- c("A", "B")



##Specify Dominant celltype prportion prportion 
##Initial Proportions if using noise parameter sc 
##If not specified will simulate between 30 to 40% for each cell type 
##Can be changed in the compute_proportions() if there are more celltypes
specified_props <- c(0.70, 0.70,  0.70, 0.70)
names(specified_props) <- c("GC", "M/TC", "OSNs", "PGC")


specified_props2 <- c(0.50, 0.50,  0.50, 0.50)
names(specified_props2) <- c("GC", "M/TC", "OSNs", "PGC")


specified_props3 <- c(0.9999, 0.99999,  0.99999, 0.99999)
names(specified_props3) <- c("GC", "M/TC", "OSNs", "PGC")


p1 <- compute_proportions6(layer_manual_MOB, pn = 0.00, sc = 0.00, specified_dominant_proportions = specified_props, neighborhood_radius = 1.69, add_labels = FALSE)


p2 <- compute_proportions6(layer_manual_MOB, pn = 0.0, sc = 0.00, specified_dominant_proportions = specified_props, neighborhood_radius = 2.6, add_labels = FALSE)


p3 <- compute_proportions6(layer_manual_MOB, pn = 0.0, sc = 0.00, specified_dominant_proportions = specified_props, neighborhood_radius = 5, add_labels = FALSE)



p4 <- compute_proportions7(layer_manual_MOB, pn = 0.0, sc = 0.00, specified_dominant_proportions = specified_props2, neighborhood_radius = 7, add_labels = FALSE)



p5 <- compute_proportions7(layer_manual_MOB, pn = 0.0, sc = 0.00, specified_dominant_proportions = specified_props3, neighborhood_radius = 7, add_labels = FALSE)



my_colors <- c( "lightblue", "orange", "red", "#9999FF")


rownames_to_keep <- c("15.076x19.996", "15.063x19.016", "15.019x18.024", "15.932x18.0", "15.893x20.047", "15.936x19.021", "12.973x18.014", "14.07x18.032", "13.059x20.044", "14.013x19.054", "14.004x20.024", "13.061x19.072", "13.969x21.05", "14.069x21.943", "13.058x22.014", "13.014x20.919", "15.89x21.016", "15.084x21.071", "15.051x21.998", "15.996x22.102", "12.08x19.06", "12.121x18.058", "12.092x20.032", "10.896x18.992", "10.877x17.97", "10.868x20.032", "12.131x21.968", "12.073x21.045", "10.868x20.971", "10.903x21.945", "10.899x23.91", "12.153x24.01", "12.115x22.986", "10.857x22.897", "12.123x25.002", "10.86x24.947", "15.949x22.955", "15.98x23.961", "15.009x22.966", "15.073x23.973", "13.048x23.951", "13.067x22.915", "14.05x23.965", "14.048x22.95", "13.111x24.919", "14.085x24.905", "15.033x24.886", "15.953x24.959", "16.987x17.981", "17.028x19.984", "16.987x18.995", "16.981x20.987", "16.965x22.101", "16.93x22.932", "16.869x23.909", "16.988x24.905", "10.906x25.928", "12.211x26.031", "13.941x25.99", "13.052x25.997", "15.065x25.942", "15.954x25.936", "16.996x25.92")


psub <- p1[rownames(p1) %in% rownames_to_keep, ]
lsub <- layer_manual_MOB[rownames(layer_manual_MOB) %in% rownames_to_keep, ]

a1 <- CARD.visualize.pie(
  proportion = psub,
  spatial_location = lsub, 
  colors = my_colors,
  radius = 0.5)



psub <- p2[rownames(p2) %in% rownames_to_keep, ]
lsub <- layer_manual_MOB[rownames(layer_manual_MOB) %in% rownames_to_keep, ]

a2 <- CARD.visualize.pie(
  proportion = psub,
  spatial_location = lsub, 
  colors = my_colors,
  radius = 0.5)



psub <- p3[rownames(p3) %in% rownames_to_keep, ]
lsub <- layer_manual_MOB[rownames(layer_manual_MOB) %in% rownames_to_keep, ]

a3 <- CARD.visualize.pie(
  proportion = psub,
  spatial_location = lsub, 
  colors = my_colors,
  radius = 0.5)




psub <- p4[rownames(p4) %in% rownames_to_keep, ]
lsub <- layer_manual_MOB[rownames(layer_manual_MOB) %in% rownames_to_keep, ]

a4 <- CARD.visualize.pie(
  proportion = psub,
  spatial_location = lsub, 
  colors = my_colors,
  radius = 0.5)


psub <- p5[rownames(p5) %in% rownames_to_keep, ]
lsub <- layer_manual_MOB[rownames(layer_manual_MOB) %in% rownames_to_keep, ]

a5 <- CARD.visualize.pie(
  proportion = psub,
  spatial_location = lsub, 
  colors = my_colors,
  radius = 0.5)

a1 + a2 + a3 + a4 + a5




CARD.visualize.pie(
  proportion = p1,
  spatial_location = layer_manual_MOB, 
  colors = my_colors,
  radius = 0.5)

















p1 <- compute_proportions6(org1, pn = 0.00, sc = 0.00, specified_dominant_proportions = specified_props, neighborhood_radius = 4, add_labels = FALSE)





p2 <- compute_proportions6(org1, pn = 0.00, sc = 0.00, specified_dominant_proportions = specified_props2, neighborhood_radius = 0, add_labels = FALSE)


p3 <- compute_proportions6(org1, pn = 0.00, sc = 0.40, specified_dominant_proportions = specified_props, neighborhood_radius = 4, add_labels = FALSE)



CARD.visualize.pie(
  proportion = p1,
  spatial_location = org1, 
  colors = c("red","blue"),
  radius = 0.5)


p1 






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
    scale_color_manual(values = c("lightblue", "orange", "red", "#9999FF")) +
    theme(plot.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),  # Increase font size and make bold
      legend.title = element_text(size = 14, face = "bold")  # Increase title font size and make bold
    )  # Adjust colors as needed
}

####Visualize dominant cell types and noise
plot_dominant_ct(p1,layer_manual_MOB, point_size = 10)

---
  title: "R Notebook"
output: html_notebook
---
  
  
  ```{r}
# install devtools if necessary
install.packages('devtools')

# install the CARD package
devtools::install_github('YingMa0107/CARD')
devtools::install_github('RGLab/MAST')
library(MAST)

# install the MuSiC package

BiocManager::install("TOAST")

devtools::install_github('xuranw/MuSiC')

if (!requireNamespace("spdep", quietly = TRUE)) install.packages("spdep")
library(spdep)
# load
library(MuSiC)
# load package
library(CARD)


library(mclust)
```


```{r}
here::i_am("Code/CARD_analysis.Rmd")

# Assigning colors to a vector
my_colors <- c( "lightblue",  
                "orange",
                "red",  # Assuming you want red for OSNs, as it's not specified
                "#9999FF")

```

```{r}

locations <- strsplit(colnames(spe), "x")
locations <- matrix(unlist(locations), ncol = 2, byrow = TRUE)
locations <- apply(locations, 2, as.numeric)
colnames(locations) <- c("x","y")

rownames(locations) <- colnames(spe)
locations <- as.data.frame(locations)

locations$x <- round(locations$x)
locations$y <- round(locations$y)
locations$Location <- rownames(locations)
dim(locations)

```

```{r}




CARD_obj = createCARDObject(
  sc_count = seurat_object1_filtered@assays$RNA$counts,
  sc_meta = seurat_object1_filtered@meta.data,
  spatial_count = spe_filtered@assays$RNA$counts,
  spatial_location = locations,
  ct.varname = "celltypes",
  ct.select = unique(seurat_object1@meta.data$celltypes),
  sample.varname = NULL,
  minCountGene = 100,
  minCountSpot = 5)

spe@assays$spatial$counts
t(Y)

CARD_obj2 = createCARDObject(
  sc_count = seurat_object1@assays$RNA$counts,
  sc_meta = seurat_object1@meta.data,
  spatial_count = spe@assays$RNA$counts,
  spatial_location = locations,
  ct.varname = "celltypes",
  ct.select = unique(seurat_object1@meta.data$celltypes),
  sample.varname = NULL,
  minCountGene = 100,
  minCountSpot = 5)

dim(seurat_object1)
dim(spe) 
```


```{r}

gene_names <- rownames(CARD_obj2@algorithm_matrix$B)
length(gene_names)

gene_names2 <- rownames(seurat_object1_filtered)

# Filter out genes starting with "Rp"
filtered_genes <- gene_names2[!grepl("^Rp", gene_names)]

# Calculate intersections
intersection_OSNS <- intersect(gene_names, genes_OSNS)
intersection_MTC <- intersect(gene_names, genes_MTC)
intersection_GC <- intersect(gene_names, genes_GC)
intersection_PGC <- intersect(gene_names, genes_PGC)



```





```{r}
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
CARD_obj2 = CARD_deconvolution(CARD_object = CARD_obj2)
CARD_proportion <- CARD_obj@Proportion_CARD
CARD_proportion2 <- CARD_obj2@Proportion_CARD
BLADE_proportion <- read.csv(here("Spatial_Blade/Results", "BLADE_autogene400_1_.csv"))
BLADE_proportion2 <- read.csv(here("Spatial_Blade/Results", "BLADE_norp_alphagroup_1_.csv"))
dim(BLADE_proportion)
rownames(BLADE_proportion) <- colnames(spe) 
rownames(BLADE_proportion2) <- colnames(spe) 

```


##MUSIC
```{r}
spe_filtered_sce <- as.SingleCellExperiment(spe_filtered)
seurat_object1_filtered_sce <- as.SingleCellExperiment(seurat_object1_filtered)

spe_sce <- data.matrix(spe_filtered@assays$RNA$counts)
colnames(spe_sce) <- colnames(spe_filtered)

Est.prop_mob = music_prop(bulk.mtx = spe_sce ,  sc.sce = seurat_object1_filtered_sce,  clusters = 'celltypes',
                          samples = 'orig.ident')


music_prop_MOB <- Est.prop_mob$Est.prop.weighted
```


```{r}
load("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/Figure3A_layer_annote (1).RData")
layer_manual_MOB

ggplot(data = layer_manual_MOB, aes(x = x, y = y)) +
  geom_point(aes(colour = Layer), size =4) + # Color points by 'Layer' # Label points by 'ID'
  theme_minimal() +
  labs(colour = "Layer")

ggplot(data = layer_manual_MOB, aes(x = x, y = y)) +
  geom_point(aes(colour = Layer), size = 3) + # Increase point size
  facet_wrap(~ Layer) + # Create separate panels for each 'Layer'
  theme_minimal() +
  labs(colour = "Layer")

intersect(colnames(spe), rownames(layer_manual_MOB))

```

```{r}
he_image <-  readJPEG(here("Data", "HE_Rep12_MOB.jpeg")) # Ensure this path is correct
he_grob <- rasterGrob(he_image, width = unit(1, "npc"), height = unit(1, "npc"))

# Assuming layer_manual_MOB is already loaded
# Step 3 is preparation, make sure layer_manual_MOB is correctly structured
names(transformed_coords_df)[names(transformed_coords_df) == "x"] <- "transformed_x"

names(transformed_coords_df)[names(transformed_coords_df) == "y"] <- "transformed_y"

# Now proceed to add the 'id' column for both dataframes, if you haven't already
transformed_coords_df$id <- rownames(transformed_coords_df)
layer_manual_MOB$id <- rownames(layer_manual_MOB)

# Merge the dataframes on 'id'
combined_df <- merge(layer_manual_MOB, transformed_coords_df, by = "id")



# Step 4: Create the overlay plot
#aspect_ratio <- dim(he_image)[2] / dim(he_image)[1] # Ensure aspect ratio is correct
spatial_plot <- ggplot() +
  annotation_custom(grob = he_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(data = combined_df, aes(x = transformed_x, y = transformed_y, color = Layer), size = 3) + # Use transformed coordinates
  coord_fixed(ratio = aspect_ratio) +
  theme_void() +
  labs(colour = "Layer")


# Display the plot
plot(spatial_plot)


```


```{r}

CARD.visualize.pie(
  proportion = CARD_proportion2,
  spatial_location = CARD_obj2@spatial_location, 
  colors = my_colors, 
  radius = 0.3)
```


```{r}
CARD.visualize.pie(
  proportion = BLADE_proportion,
  spatial_location = locations, 
  colors = my_colors, 
  radius = 0.3)
```
```{r}
CARD.visualize.pie(
  proportion = music_prop_MOB,
  spatial_location = locations, 
  colors = my_colors, 
  radius = 0.3)
```

```{r}
p2 <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = c("GC","PGC","OSNs","M/TC"),                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 4,                                 ### number of columns in the figure panel
  pointSize = 1.0)                             ### point size in ggplot2 scatterplot  
print(p2)
```


```{r}
#plot(CARD_obj@Proportion_CARD[,3], blade_result[,3])
```


```{r}
p3 <- CARD.visualize.prop(
  proportion = BLADE_proportion,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = c("GC","PGC","OSNs","M.TC"),                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 4,                                 ### number of columns in the figure panel
  pointSize = 1.0)   

### point size in ggplot2 scatterplot  
print(p3)
```

```{r}
CARD.visualize.prop(
  proportion = music_prop_MOB,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = c("GC","PGC","OSNs","M/TC"),                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 4,                                 ### number of columns in the figure panel
  pointSize = 1.0)   
```


```{r}
# Calculate the index of the dominant cell type for each spot
dominant_CARD_idx <- apply(CARD_proportion2, 1, which.max)
dominant_BLADE_idx <- apply(BLADE_proportion1, 1, which.max)

# Convert indices to cell type names
cell_type_names <- unique(seurat_object1@meta.data$celltypes)
dominant_CARD_types <- cell_type_names[dominant_CARD_idx]
dominant_BLADE_types <- cell_type_names[dominant_BLADE_idx]

# Ensure the names of the dominant types match the spot locations
names(dominant_CARD_types) <- rownames(CARD_proportion2)
names(dominant_BLADE_types) <- rownames(BLADE_proportion)

# Create a dataframe for CARD results
df_CARD <- data.frame(Spot = names(dominant_CARD_types), DominantType = dominant_CARD_types, stringsAsFactors = FALSE)
df_CARD <- merge(df_CARD, locations, by.x = "Spot", by.y = "row.names", all.x = TRUE)

# Create a dataframe for BLADE results
df_BLADE <- data.frame(Spot = names(dominant_BLADE_types), DominantType = dominant_BLADE_types, stringsAsFactors = FALSE)
df_BLADE <- merge(df_BLADE, locations, by.x = "Spot", by.y = "row.names", all.x = TRUE)

dominant_MUSIC_idx <- apply(music_prop_MOB, 1, which.max)

# Convert indices to cell type names
dominant_MUSIC_types <- cell_type_names[dominant_MUSIC_idx]

# Ensure the names of the dominant types match the spot locations
names(dominant_MUSIC_types) <- rownames(music_prop_MOB)

# Create a dataframe for MUSIC results
df_MUSIC <- data.frame(Spot = names(dominant_MUSIC_types), DominantType = dominant_MUSIC_types, stringsAsFactors = FALSE)
df_MUSIC <- merge(df_MUSIC, locations, by.x = "Spot", by.y = "row.names", all.x = TRUE)


# Function to calculate PC1 and merge it with the existing data frame
calculate_pc1_and_merge <- function(df, composition_matrix) {
  pca_result <- prcomp(composition_matrix, scale. = TRUE)
  pc1_values <- pca_result$x[, 1]
  df$PC1 <- pc1_values[match(df$Spot, rownames(composition_matrix))]
  return(df)
}

# Apply the function to each dataset
df_CARD <- calculate_pc1_and_merge(df_CARD, CARD_proportion2)
df_BLADE <- calculate_pc1_and_merge(df_BLADE, BLADE_proportion)
df_MUSIC <- calculate_pc1_and_merge(df_MUSIC, music_prop_MOB)

# Common plot theme
common_theme <- theme_minimal() + 
  theme(legend.position = "bottom")

# Common labs
common_labs <- labs(color = "Cell Type / Layer")

# Common geom_point size
point_size <- 5
```
```{r}

# Plot for CARD
p_CARD <- ggplot(df_CARD, aes(x = x, y = y, colour = PC1)) +
  geom_point(size = point_size) +
  common_theme +
  labs(title = "First Principal Component - CARD") +
  scale_colour_gradient(low = "blue", high = "red", )

p_BLADE <- ggplot(df_BLADE, aes(x = x, y = y, colour = PC1)) +
  geom_point(size = point_size) +
  common_theme +
  labs(title = "First Principal Component - BLADE") +
  scale_colour_gradient(low = "blue", high = "red")


ggplot(df_MUSIC, aes(x = x, y = y, colour = PC1)) +
  geom_point(size = point_size) +
  common_theme +
  labs(title = "First Principal Component - MUSIC") +
  scale_colour_gradient(low = "blue", high = "red")


```


```{r}
ggplot(df_CARD, aes(x = x, y = y)) +
  geom_point(aes(colour = DominantType), size = 5) + 
  common_theme +
  labs(title = "Dominant Cell Types - CARD") +
  common_labs +
  scale_color_manual(values = c("GC" = "#FF9999",    # lighter red
                                "M/TC" = "#99CC99",  # lighter green
                                "OSNs" = "#9999FF",  # lighter blue
                                "PGC" = "#D8BFD8"))  # lighter purple (Thistle)

```


```{r}
ggplot(df_BLADE, aes(x = x, y = y, color = DominantType)) +
  geom_point(aes(colour = DominantType), size = 5) + 
  common_theme +
  labs(title = "Dominant Cell Types - BLADE") +
  common_labs +
  scale_color_manual(values = c("GC" = "#FF9999",    # lighter red
                                "M/TC" = "#99CC99",  # lighter green
                                "OSNs" = "#9999FF",  # lighter blue
                                "PGC" = "#D8BFD8"))  # lighter purple (Thistle)

```
```{r}
ggplot(df_MUSIC, aes(x = x, y = y, color = DominantType)) +
  geom_point(aes(colour = DominantType), size = 5) + 
  common_theme +
  labs(title = "Dominant Cell Types - MUSIC") +
  common_labs +
  scale_color_manual(values = c("GC" = "#FF9999",    # lighter red
                                "M/TC" = "#99CC99",  # lighter green
                                "OSNs" = "#9999FF",  # lighter blue
                                "PGC" = "#D8BFD8"))
```


```{r}
layer_to_cell_type <- list(
  GCL = "GC",
  GL = "PGC",
  ONL = "OSNs",
  MCL = "M.TC"
)

layer_manual_MOB$Spot <- rownames(layer_manual_MOB)
layer_manual_MOB$celltypes <- sapply(layer_manual_MOB$Layer, function(x) layer_to_cell_type[[x]])

ggplot(data = layer_manual_MOB, aes(x = x, y = y, color = celltypes)) +
  geom_point(aes(colour = celltypes), size = 5) +
  common_theme +
  labs(title = "Mouse olfactory bulb annotated layers", color = "Layer") +
  common_labs +
  scale_color_manual(values = c("GC" = "#FF9999",    # lighter red
                                "M.TC" = "#99CC99",  # lighter green
                                "OSNs" = "#9999FF",  # lighter blue
                                "PGC" = "#D8BFD8"))  # lighter purple (Thistle)


```
```{r}
df_CARD <- df_CARD[order(df_CARD$Spot),]
df_BLADE <- df_BLADE[order(df_BLADE$Spot),]
df_MUSIC <- df_MUSIC[order(df_MUSIC$Spot),]
```


```{r}
# Reorder df_CARD to match the order of spots in layer_manual_MOB
df_CARD_ordered <- df_CARD[match(layer_manual_MOB$Spot, df_CARD$Spot),]

# Reorder df_BLADE to match the order of spots in layer_manual_MOB
df_BLADE_ordered <- df_BLADE[match(layer_manual_MOB$Spot, df_BLADE$Spot),]

# Reorder df_BLADE to match the order of spots in layer_manual_MOB
df_MUSIC_ordered <- df_MUSIC[match(layer_manual_MOB$Spot, df_MUSIC$Spot),]


# Verify alignment for df_CARD
all(layer_manual_MOB$Spot == df_CARD_ordered$Spot)  # Should return TRUE

# Verify alignment for df_BLADE
all(layer_manual_MOB$Spot == df_BLADE_ordered$Spot)  # Should return TRUE


# Verify alignment for df_CARD
all(layer_manual_MOB$Spot == df_MUSIC_ordered$Spot)

```
```{r}
predicted_CARD <- df_CARD_ordered$DominantType
predicted_BLADE <- df_BLADE_ordered$DominantType
predicted_MUSIC <- df_MUSIC_ordered$DominantType

# Assuming you're comparing against a specific mapped cell type or direct layer classification
true_labels <- layer_manual_MOB$celltypes # Or layer_manual_MOB$Layer if directly comparable

```



```{r}
# Calculate ARI
ARI_CARD <- mclust::adjustedRandIndex(true_labels, predicted_CARD)
ARI_BLADE <- mclust::adjustedRandIndex(true_labels, predicted_BLADE)

ARI_MUSIC <- mclust::adjustedRandIndex(true_labels, predicted_MUSIC)

# Print ARI results
print(paste("ARI for CARD:", ARI_CARD))
print(paste("ARI for BLADE:", ARI_BLADE))

# Assuming ARI_CARD and ARI_BLADE are single numeric values
ari_data <- data.frame(
  Method = c("CARD", "BLADE", "MUSIC"),
  ARI_Value = c(ARI_CARD, ARI_BLADE, ARI_MUSIC)
)

# Use ggplot2 to create the bar chart
library(ggplot2)

ggplot(ari_data, aes(x = Method, y = ARI_Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylim(0, 1) +  # ARI values range between 0 and 1
  labs(title = "Adjusted Rand Index (ARI) for CARD and BLADE",
       x = "Method",
       y = "ARI Value") +
  theme_minimal() +
  theme(legend.position = "none") +  # Hide the legend since colors are apparent
  geom_text(aes(label = round(ARI_Value, 2)), vjust = -0.3)  # Add value labels above bars


```

```{r}
calculateNodePurity <- function(true_labels, predicted_labels) {
  combined <- table(true_labels, predicted_labels)
  sum_max_per_cluster <- sum(apply(combined, 2, max))  # Sum of max counts per predicted cluster
  total_predictions <- sum(combined)  # Total number of predictions
  purity <- sum_max_per_cluster / total_predictions
  return(purity)
}

# Calculate Node Purity
Purity_CARD <- calculateNodePurity(true_labels, predicted_CARD)
Purity_BLADE <- calculateNodePurity(true_labels, predicted_BLADE)
Purity_MUSIC <- calculateNodePurity(true_labels, predicted_MUSIC)

# Print Node Purity results

print(paste("Node Purity for CARD:", Purity_CARD))
print(paste("Node Purity for BLADE:", Purity_BLADE))
print(paste("Node Purity for MUSIC:", Purity_MUSIC))

purity_data <- data.frame(
  Method = c("CARD", "BLADE","MUSIC"),
  Purity_Value = c(Purity_CARD, Purity_BLADE, Purity_MUSIC)
)


ggplot(purity_data, aes(x = Method, y = Purity_Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylim(0, 1) +  # Purity values theoretically range between 0 and 1
  labs(title = "Node Purity",
       x = "Method",
       y = "Purity Value") +
  theme_minimal() +
  theme(legend.position = "none") +  # Hide the legend since colors are apparent
  geom_text(aes(label = round(Purity_Value, 2)), vjust = -0.3)  # Add value labels above bars




table(layer_manual_MOB$Layer)

```


```{r}
``
```

```{r}
coords <- locations[, c("x", "y")]

nb <- knn2nb(knearneigh(coords, k = 4), sym = TRUE)

# Convert neighbors to weights
weights <- nb2listw(nb, style = "W", zero.policy = TRUE)

```

```{r}
# Convert categorical variables into dummy variables
layer_manual_MOB_dummies <- model.matrix(~ layer_manual_MOB$celltype - 1)
colnames(layer_manual_MOB_dummies) <- c("GC", "M.TC", "OSNs" ,"PGC")

coords_gt<- layer_manual_MOB[,c("x","y")]
nbgt <- knn2nb(knearneigh(coords_gt, k = 4), sym = TRUE)

# Convert neighbors to weights
weights_gt <- nb2listw(nbgt, style = "W", zero.policy = TRUE)


# Create a dataframe to store Moran's I results
moran_results_gt <- data.frame(CellType = character(), MoranI = numeric(), p.value = numeric())


# Loop through each dummy variable column in layer_manual_MOB_dummies
for (i in 1:ncol(layer_manual_MOB_dummies)) {
  # Perform Moran's I test
  moran_result_gt <- moran.test(layer_manual_MOB_dummies[, i], weights_gt, zero.policy = TRUE)
  # Add results to the dataframe
  moran_results_gt <- rbind(moran_results_gt, data.frame(
    CellType = colnames(layer_manual_MOB_dummies)[i],
    MoranI = moran_result_gt$estimate["Moran I statistic"],
    p.value = moran_result_gt$p.value
  ))
}

ggplot(moran_results_gt, aes(x = CellType, y = MoranI, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.2f", MoranI)), vjust = -0.5, color = "black") +
  labs(title = "Moran's I Statistics for Each Cell Type for Ground Truth",
       x = "Cell Type",
       y = "Moran's I") +
  scale_fill_manual(values = c("GC" = "red", "M.TC" = "green", "OSNs" = "blue", "PGC" = "grey"))+
  theme_minimal() +
  theme(legend.position = "none") 


```


```{r}
moran_results <- data.frame(CellType = character(), MoranI = numeric(), p.value = numeric())

# Loop through each cell type column in CARD_proportion
for (i in 1:ncol(CARD_proportion2)) {
  # Perform Moran's I test
  moran_result <- moran.test(CARD_proportion2[, i], weights, zero.policy = TRUE)
  
  # Add results to the dataframe
  moran_results <- rbind(moran_results, data.frame(
    CellType = colnames(CARD_proportion2)[i],
    MoranI = moran_result$estimate["Moran I statistic"],
    p.value = moran_result$p.value
  ))
}


ggplot(moran_results, aes(x = CellType, y = MoranI, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.2f", MoranI)), vjust = -0.5, color = "black") +
  labs(title = "Moran's I Statistics for Each Cell Type for CARD",
       x = "Cell Type",
       y = "Moran's I") +
  scale_fill_manual(values = c("GC" = "red", "M/TC" = "green", "OSNs" = "blue", "PGC" = "grey"))+
  theme_minimal() +
  theme(legend.position = "none") 

```

```{r}
local_morans_list <- list()

# Loop through each cell type column in CARD_proportion2
for (i in 1:ncol(CARD_proportion2)) {
  # Calculate Local Moran's I
  local_moran <- localmoran(CARD_proportion2[, i], weights)
  
  # Store results in the list
  local_morans_list[[colnames(CARD_proportion2)[i]]] <- local_moran
}

local_moran_gc <- local_morans_list[["GC"]]
coords$local_moran_I_GC <- local_moran_gc[,1]
coords$significance_GC <- local_moran_gc[,5] < 0.05

# Plotting Local Moran's I for 'GC'
plot_gc <- ggplot(coords, aes(x = x, y = y, color = significance_GC)) + 
  geom_point(aes(size = local_moran_I_GC)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for GC",
       x = "X Coordinate",
       y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")

# Extract Local Moran's I results for 'M/TC'
local_moran_m_tc <- local_morans_list[["M/TC"]]
coords$local_moran_I_M_TC <- local_moran_m_tc[,1]
coords$p_value_M_TC <- local_moran_m_tc[,5]

# Create a binary significance column for 'M/TC'
coords$significance_M_TC <- coords$p_value_M_TC < 0.05

# Plot for 'M/TC'
plot_m_tc <- ggplot(coords, aes(x = x, y = y, color = significance_M_TC)) + 
  geom_point(aes(size = local_moran_I_M_TC)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for M/TC (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")

# Extract Local Moran's I results for 'OSNs'
local_moran_osns <- local_morans_list[["OSNs"]]
coords$local_moran_I_OSNs <- local_moran_osns[,1]
coords$p_value_OSNs <- local_moran_osns[,5]

# Create a binary significance column for 'OSNs'
coords$significance_OSNs <- coords$p_value_OSNs < 0.05

# Plot for 'OSNs'
plot_osns <- ggplot(coords, aes(x = x, y = y, color = significance_OSNs)) + 
  geom_point(aes(size = local_moran_I_OSNs)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for OSNs (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")


local_moran_pgc <- local_morans_list[["PGC"]]
coords$local_moran_I_PGC <- local_moran_pgc[,1]
coords$p_value_PGC <- local_moran_pgc[,5]

# Create a binary significance column for 'OSNs'
coords$significance_PGC <- coords$p_value_PGC < 0.05

# Plot for 'OSNs'
plot_PGC <- ggplot(coords, aes(x = x, y = y, color = significance_PGC)) + 
  geom_point(aes(size = local_moran_I_OSNs)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for PGCs (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")




```



```{r}



```





```{r}
moran_results <- data.frame(CellType = character(), MoranI = numeric(), p.value = numeric())

# Loop through each cell type column in CARD_proportion
for (i in 1:ncol(BLADE_proportion)) {
  # Perform Moran's I test
  moran_result <- moran.test(BLADE_proportion[, i], weights, zero.policy = TRUE)
  
  # Add results to the dataframe
  moran_results <- rbind(moran_results, data.frame(
    CellType = colnames(BLADE_proportion)[i],
    MoranI = moran_result$estimate["Moran I statistic"],
    p.value = moran_result$p.value
  ))
}

ggplot(moran_results, aes(x = CellType, y = MoranI, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.2f", MoranI)), vjust = -0.5, color = "black") +
  labs(title = "Moran's I Statistics for Each Cell Type for BLADE",
       x = "Cell Type",
       y = "Moran's I") +
  scale_fill_manual(values = c("GC" = "red", "M.TC" = "green", "OSNs" = "blue", "PGC" = "grey"))+
  theme_minimal() +
  theme(legend.position = "none") 




```

```{r}
local_morans_list <- list()

# Loop through each cell type column in CARD_proportion2
for (i in 1:ncol(BLADE_proportion)) {
  # Calculate Local Moran's I
  local_moran <- localmoran(BLADE_proportion[, i], weights)
  
  # Store results in the list
  local_morans_list[[colnames(BLADE_proportion)[i]]] <- local_moran
}

local_moran_gc <- local_morans_list[["GC"]]
coords$local_moran_I_GC <- local_moran_gc[,1]
coords$significance_GC <- local_moran_gc[,5] < 0.05

# Plotting Local Moran's I for 'GC'
plot_gc <- ggplot(coords, aes(x = x, y = y, color = significance_GC)) + 
  geom_point(aes(size = local_moran_I_GC)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for GC",
       x = "X Coordinate",
       y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")

# Extract Local Moran's I results for 'M/TC'
local_moran_m_tc <- local_morans_list[["M.TC"]]
coords$local_moran_I_M_TC <- local_moran_m_tc[,1]
coords$p_value_M_TC <- local_moran_m_tc[,5]

# Create a binary significance column for 'M/TC'
coords$significance_M_TC <- coords$p_value_M_TC < 0.05

# Plot for 'M/TC'
plot_m_tc <- ggplot(coords, aes(x = x, y = y, color = significance_M_TC)) + 
  geom_point(aes(size = local_moran_I_M_TC)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for M/TC (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")

# Extract Local Moran's I results for 'OSNs'
local_moran_osns <- local_morans_list[["OSNs"]]
coords$local_moran_I_OSNs <- local_moran_osns[,1]
coords$p_value_OSNs <- local_moran_osns[,5]

# Create a binary significance column for 'OSNs'
coords$significance_OSNs <- coords$p_value_OSNs < 0.05

# Plot for 'OSNs'
plot_osns <- ggplot(coords, aes(x = x, y = y, color = significance_OSNs)) + 
  geom_point(aes(size = local_moran_I_OSNs)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for OSNs (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")


local_moran_pgc <- local_morans_list[["PGC"]]
coords$local_moran_I_PGC <- local_moran_pgc[,1]
coords$p_value_PGC <- local_moran_pgc[,5]

# Create a binary significance column for 'OSNs'
coords$significance_PGC <- coords$p_value_PGC < 0.05

# Plot for 'OSNs'
plot_PGC <- ggplot(coords, aes(x = x, y = y, color = significance_PGC)) + 
  geom_point(aes(size = local_moran_I_OSNs)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for PGCs (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")


```


```{r}
moran_results <- data.frame(CellType = character(), MoranI = numeric(), p.value = numeric())

# Loop through each cell type column in CARD_proportion
for (i in 1:ncol(music_prop_MOB)) {
  # Perform Moran's I test
  moran_result <- moran.test(music_prop_MOB[, i], weights, zero.policy = TRUE)
  
  # Add results to the dataframe
  moran_results <- rbind(moran_results, data.frame(
    CellType = colnames(music_prop_MOB)[i],
    MoranI = moran_result$estimate["Moran I statistic"],
    p.value = moran_result$p.value
  ))
}

ggplot(moran_results, aes(x = CellType, y = MoranI, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.2f", MoranI)), vjust = -0.5, color = "black") +
  labs(title = "Moran's I Statistics for Each Cell Type for MUSIC",
       x = "Cell Type",
       y = "Moran's I") +
  scale_fill_manual(values = c("GC" = "red", "M/TC" = "green", "OSNs" = "blue", "PGC" = "grey"))+
  theme_minimal() +
  theme(legend.position = "none") 



```


```{r}
local_morans_list <- list()

# Loop through each cell type column in CARD_proportion2
for (i in 1:ncol(music_prop_MOB)) {
  # Calculate Local Moran's I
  local_moran <- localmoran(music_prop_MOB[, i], weights)
  
  # Store results in the list
  local_morans_list[[colnames(music_prop_MOB)[i]]] <- local_moran
}

local_moran_gc <- local_morans_list[["GC"]]
coords$local_moran_I_GC <- local_moran_gc[,1]
coords$significance_GC <- local_moran_gc[,5] < 0.05

# Plotting Local Moran's I for 'GC'
plot_gc <- ggplot(coords, aes(x = x, y = y, color = significance_GC)) + 
  geom_point(aes(size = local_moran_I_GC)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for GC",
       x = "X Coordinate",
       y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")

# Extract Local Moran's I results for 'M/TC'
local_moran_m_tc <- local_morans_list[["M/TC"]]
coords$local_moran_I_M_TC <- local_moran_m_tc[,1]
coords$p_value_M_TC <- local_moran_m_tc[,5]

# Create a binary significance column for 'M/TC'
coords$significance_M_TC <- coords$p_value_M_TC < 0.05

# Plot for 'M/TC'
plot_m_tc <- ggplot(coords, aes(x = x, y = y, color = significance_M_TC)) + 
  geom_point(aes(size = local_moran_I_M_TC)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for M/TC (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")

# Extract Local Moran's I results for 'OSNs'
local_moran_osns <- local_morans_list[["OSNs"]]
coords$local_moran_I_OSNs <- local_moran_osns[,1]
coords$p_value_OSNs <- local_moran_osns[,5]

# Create a binary significance column for 'OSNs'
coords$significance_OSNs <- coords$p_value_OSNs < 0.05

# Plot for 'OSNs'
plot_osns <- ggplot(coords, aes(x = x, y = y, color = significance_OSNs)) + 
  geom_point(aes(size = local_moran_I_OSNs)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for OSNs (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")


local_moran_pgc <- local_morans_list[["PGC"]]
coords$local_moran_I_PGC <- local_moran_pgc[,1]
coords$p_value_PGC <- local_moran_pgc[,5]

# Create a binary significance column for 'OSNs'
coords$significance_PGC <- coords$p_value_PGC < 0.05

# Plot for 'OSNs'
plot_PGC <- ggplot(coords, aes(x = x, y = y, color = significance_PGC)) + 
  geom_point(aes(size = local_moran_I_OSNs)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Moran's I for PGCs (Significant Spots Highlighted)",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "none")

```
















```{r}
card_long <- CARD_proportion2 %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Location") %>%
  pivot_longer(cols = -Location, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Method = "CARD")

blade_long <- BLADE_proportion %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Location") %>%
  pivot_longer(cols = -Location, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Method = "BLADE")

# Merge with location data
card_merged <- left_join(card_long, locations, by = "Location")
blade_merged <- left_join(blade_long, locations, by = "Location")

# Combine CARD and BLADE data
combined_data <- bind_rows(card_merged, blade_merged)

music_long <- music_prop_MOB %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Location") %>%
  pivot_longer(cols = -Location, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Method = "MUSIC")

# Merge with location data
music_merged <- left_join(music_long, locations, by = "Location")

# Combine CARD, BLADE, and MUSIC data
combined_data_all <- bind_rows(combined_data, music_merged)

# Update the plot to include MUSIC data
p_all <- ggplot(combined_data_all, aes(x = x, y = y, fill = Proportion)) +
  geom_point(shape = 21, size = 3.7) +  # shape 21 is a filled circle, adjust size as needed
  facet_wrap(~Method + CellType, ncol = 4) +  # Adjust the ncol if necessary
  scale_fill_gradient2(low = "black", high = "red", mid = "lightgreen", midpoint = 0.5, 
                       limit = c(0, 1), name = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "right") +  # Include legend on the right
  labs(title = "Cell Type Proportions for CARD, BLADE, and MUSIC")

# Display the updated plot
print(p_all)
```
```{r}


CARD.visualize.Cor(BLADE_proportion,colors = NULL) # if not provide, we will use the default colors



```



Correllation with marker gene expression and extimated proportions of the cell types 

```{r}
# Extract the gene expression values for each cell type marker gene
marker_genes <- list(
  GC = genes_GC_25,
  PGC = genes_PGC_25,
  OSNs = genes_OSNS_25,
  M.TC = genes_MTC_25
)


# Calculate the correlation between the gene expression values and the estimated cell type proportions
correlation_results <- data.frame(CellType = character(), Gene = character(), Correlation = numeric())

for (cell_type in names(marker_genes)) {
  if (cell_type %in% colnames(BLADE_proportion)) {
    for (gene in marker_genes[[cell_type]]) {
      # Extract gene expression values
      if(gene %in% rownames(GetAssayData(object = spe, assay = "RNA", slot = "data"))) {
        gene_expression <- GetAssayData(object = spe, assay = "RNA", slot = "data")[gene, ]
      } else {
        next  # Skip this gene if it's not found
      }
      
      # Extract cell type proportions
      cell_type_proportions <- BLADE_proportion[, cell_type]
      
      # Check dimensions
      if (length(gene_expression) == length(cell_type_proportions)) {
        # Calculate correlation
        correlation <- cor(gene_expression, cell_type_proportions)
        
        # Store results
        correlation_results <- rbind(correlation_results, data.frame(
          CellType = cell_type,
          Gene = gene,
          Correlation = correlation
        ))
      } else {
        print(paste("Mismatched dimensions for gene", gene, "in cell type", cell_type))
      }
    }
  } else {
    print(paste("Cell type", cell_type, "not found in BLADE_proportion"))
  }
}


split_data <- split(correlation_results, correlation_results$CellType)

total_correlation <- aggregate(Correlation ~ CellType, correlation_results, mean)

# Print total correlation values
print(total_correlation)

plots <- lapply(names(split_data), function(cell_type) {
  total_corr <- total_correlation$Correlation[total_correlation$CellType == cell_type]
  ggplot(split_data[[cell_type]], aes(x = Gene, y = Correlation)) +
    geom_bar(stat = "identity", aes(fill = Correlation), position = "dodge") +  # Fill bars based on Correlation
    scale_fill_gradient(low = "blue", high = "red") +  # Gradient from low (blue) to high (red) correlation
    geom_text(aes(label = round(Correlation, 2), y = Correlation + 0.02), vjust = 0, position = position_dodge(width = 0.9), size = 3) +  # Add correlation values on top
    annotate("text", x = Inf, y = Inf, label = paste("Total Correlation:", round(total_corr, 2)), hjust = 1, vjust = 1, size = 5) + # Add total correlation
    labs(title = paste("Correlation in", cell_type),
         x = "Marker Gene", y = "Correlation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none")  # Remove legend for simplicity
})



do.call(grid.arrange, c(plots, ncol = 2))

```


```{r}
rm(marker_genes, correlation_results, scatterplots)

```

```{r}
marker_genes <- list(
  GC = genes_GC_25,
  PGC = genes_PGC_25,
  OSNs = genes_OSNS_25,
  "M/TC" = genes_MTC_25
)


# Calculate the correlation between the gene expression values and the estimated cell type proportions
correlation_results <- data.frame(CellType = character(), Gene = character(), Correlation = numeric())

for (cell_type in names(marker_genes)) {
  if (cell_type %in% colnames(music_prop_MOB)) {
    for (gene in marker_genes[[cell_type]]) {
      # Extract gene expression values
      if(gene %in% rownames(GetAssayData(object = spe, assay = "RNA", slot = "data"))) {
        gene_expression <- GetAssayData(object = spe, assay = "RNA", slot = "data")[gene, ]
      } else {
        next  # Skip this gene if it's not found
      }
      
      # Extract cell type proportions
      cell_type_proportions <- music_prop_MOB[, cell_type]
      
      # Check dimensions
      if (length(gene_expression) == length(cell_type_proportions)) {
        # Calculate correlation
        correlation <- cor(gene_expression, cell_type_proportions)
        
        # Store results
        correlation_results <- rbind(correlation_results, data.frame(
          CellType = cell_type,
          Gene = gene,
          Correlation = correlation
        ))
      } else {
        print(paste("Mismatched dimensions for gene", gene, "in cell type", cell_type))
      }
    }
  } else {
    print(paste("Cell type", cell_type, "not found in BLADE_proportion"))
  }
}


split_data <- split(correlation_results, correlation_results$CellType)


total_correlation <- aggregate(Correlation ~ CellType, correlation_results, mean)

# Print total correlation values
print(total_correlation)

plots <- lapply(names(split_data), function(cell_type) {
  total_corr <- total_correlation$Correlation[total_correlation$CellType == cell_type]
  ggplot(split_data[[cell_type]], aes(x = Gene, y = Correlation)) +
    geom_bar(stat = "identity", aes(fill = Correlation), position = "dodge") +  # Fill bars based on Correlation
    scale_fill_gradient(low = "blue", high = "red") +  # Gradient from low (blue) to high (red) correlation
    geom_text(aes(label = round(Correlation, 2), y = Correlation + 0.02), vjust = 0, position = position_dodge(width = 0.9), size = 3) +  # Add correlation values on top
    annotate("text", x = Inf, y = Inf, label = paste("Total Correlation:", round(total_corr, 2)), hjust = 1, vjust = 1, size = 5) + # Add total correlation
    labs(title = paste("Correlation in", cell_type),
         x = "Marker Gene", y = "Correlation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none")  # Remove legend for simplicity
})

do.call(grid.arrange, c(plots, ncol = 2))
```

```{r}



```





```{r warning=TRUE}
# Extract the gene expression values for each cell type marker gene
marker_genes <- list(
  GC = genes_GC_25,
  PGC = genes_PGC_25,
  OSNs = genes_OSNS_25,
  "M/TC" = genes_MTC_25
)

# Calculate the correlation between the gene expression values and the estimated cell type proportions
correlation_results <- data.frame(CellType = character(), Gene = character(), Correlation = numeric())

for (cell_type in names(marker_genes)) {
  if (cell_type %in% colnames(CARD_proportion2)) {
    for (gene in marker_genes[[cell_type]]) {
      
      # Extract gene expression values
      if(gene %in% rownames(GetAssayData(object = spe, assay = "RNA", slot = "data"))) {
        gene_expression <- GetAssayData(object = spe, assay = "RNA", slot = "data")[gene, ]
      } else {
        next  # Skip this gene if it's not found
      }
      
      # Extract cell type proportions
      cell_type_proportions <- CARD_proportion2[, cell_type]
      
      # Check dimensions
      if (length(gene_expression) == length(cell_type_proportions)) {
        # Calculate correlation
        correlation <- cor(gene_expression, cell_type_proportions)
        
        # Store results
        correlation_results <- rbind(correlation_results, data.frame(
          CellType = cell_type,
          Gene = gene,
          Correlation = correlation
        ))
      } else {
        print(paste("Mismatched dimensions for gene", gene, "in cell type", cell_type))
      }
    }
  } else {
    print(paste("Cell type", cell_type, "not found in CARD_proportions2"))
  }
}


split_data <- split(correlation_results, correlation_results$CellType)


total_correlation <- aggregate(Correlation ~ CellType, correlation_results, mean)

# Print total correlation values
print(total_correlation)

plots <- lapply(names(split_data), function(cell_type) {
  total_corr <- total_correlation$Correlation[total_correlation$CellType == cell_type]
  ggplot(split_data[[cell_type]], aes(x = Gene, y = Correlation)) +
    geom_bar(stat = "identity", aes(fill = Correlation), position = "dodge") +  # Fill bars based on Correlation
    scale_fill_gradient(low = "blue", high = "red") +  # Gradient from low (blue) to high (red) correlation
    geom_text(aes(label = round(Correlation, 2), y = Correlation + 0.02), vjust = 0, position = position_dodge(width = 0.9), size = 3) +  # Add correlation values on top
    annotate("text", x = Inf, y = Inf, label = paste("Total Correlation:", round(total_corr, 2)), hjust = 1, vjust = 1, size = 5) + # Add total correlation
    labs(title = paste("Correlation in", cell_type),
         x = "Marker Gene", y = "Correlation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none")  # Remove legend for simplicity
})


do.call(grid.arrange, c(plots, ncol = 2))



```
```{r}
```



```{r}
CARD.visualize.Cor(CARD_proportion2,colors = NULL)

CARD.visualize.Cor(BLADE_proportion, colors = NULL)

CARD.visualize.Cor(music_prop_MOB, colors = NULL)

```


```{r warning=FALSE}
gene_names <- rownames(GetAssayData(seurat_object1_filtered, assay = "RNA", slot = "data"))
# Split the Seurat object by cell type
seurat_list <- SplitObject(seurat_object1_filtered, split.by = "celltypes")

# Apply the function to each cell type subset
mean_expr_list <- lapply(seurat_list, calc_mean_expr)

# Combine into a single data matrix
mean_expr_matrix <- do.call(cbind, mean_expr_list)

# Name the columns as cell types
colnames(mean_expr_matrix) <- names(mean_expr_list)

# Assign gene names to the rows (using the same gene_names variable from before)
rownames(mean_expr_matrix) <- gene_names

mean_expr_matrix <- as.data.frame(mean_expr_matrix)
# Assuming mean_expr_list is your dataframe and it has at least 4 columns
# Scatterplot between column 2 and 3
create_scatter_plot <- function(data, x_col_name, y_col_name) {
  ggplot(data, aes(x = .data[[x_col_name]], y = .data[[y_col_name]])) +
    geom_point() +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    theme_minimal() +
    labs(title = paste("Scatterplot between", x_col_name, "and", y_col_name),
         x = x_col_name,
         y = y_col_name)
}

plot_list <- list()

# Names of the columns you want to plot
col_names <- colnames(mean_expr_matrix)

# Loop over combinations of columns and create plots
for (i in 1:(length(col_names) - 1)) {
  for (j in (i+1):length(col_names)) {
    x_col_name <- col_names[i]
    y_col_name <- col_names[j]
    plot_list[[paste(x_col_name, y_col_name)]] <- create_scatter_plot(mean_expr_matrix, x_col_name, y_col_name)
  }
}

# Arrange all plots in a grid
do.call(grid.arrange, c(plot_list, ncol = 2))




```

```{r warning=FALSE}

gene_names <- rownames(GetAssayData(seurat_object1, assay = "RNA", slot = "data"))
# Split the Seurat object by cell type
seurat_list <- SplitObject(seurat_object1, split.by = "celltypes")
# Apply the function to each cell type subset
mean_expr_list <- lapply(seurat_list, calc_mean_expr)

# Combine into a single data matrix
mean_expr_matrix <- do.call(cbind, mean_expr_list)

# Name the columns as cell types
colnames(mean_expr_matrix) <- names(mean_expr_list)

# Assign gene names to the rows (using the same gene_names variable from before)
rownames(mean_expr_matrix) <- gene_names

# Apply the function to each cell type subset
mean_expr_list <- lapply(seurat_list, calc_mean_expr)

# Combine into a single data matrix
mean_expr_matrix <- do.call(cbind, mean_expr_list)

# Name the columns as cell types
colnames(mean_expr_matrix) <- names(mean_expr_list)

# Assign gene names to the rows (using the same gene_names variable from before)
rownames(mean_expr_matrix) <- gene_names


mean_expr_matrix <- as.data.frame(mean_expr_matrix)
# Assuming mean_expr_list is your dataframe and it has at least 4 columns
# Scatterplot between column 2 and 3
# Function to create a scatterplot with a trend line
create_scatter_plot <- function(data, x_col_name, y_col_name) {
  ggplot(data, aes(x = .data[[x_col_name]], y = .data[[y_col_name]])) +
    geom_point() +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    theme_minimal() +
    labs(title = paste("Scatterplot between", x_col_name, "and", y_col_name),
         x = x_col_name,
         y = y_col_name)
}

plot_list <- list()

# Names of the columns you want to plot
col_names <- colnames(mean_expr_matrix)

# Loop over combinations of columns and create plots
for (i in 1:(length(col_names) - 1)) {
  for (j in (i+1):length(col_names)) {
    x_col_name <- col_names[i]
    y_col_name <- col_names[j]
    plot_list[[paste(x_col_name, y_col_name)]] <- create_scatter_plot(mean_expr_matrix, x_col_name, y_col_name)
  }
}

# Arrange all plots in a grid
do.call(grid.arrange, c(plot_list, ncol = 2))



```


```{r}
correlation_matrix <- cor(CARD_obj@Proportion_CARD[:,], blade_result[:,])

# Load the necessary library
library(pheatmap)

# Create the heatmap
pheatmap(correlation_matrix, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         display_numbers = TRUE, # Set to FALSE to hide numbers in cells
         fontsize_number = 8, # Adjust font size for the numbers
         legend = TRUE
)
```


```{r}
gene1 <- c('Gm26901', 'Lman2l', 'Ormdl1', 'Fam117b', 'Fzd5', 'Igfbp5', 'Scg2', 'Dock10', 'Krtap28-13', 'Itm2c', 'Tmem183a', 'B3galt2', 'Vamp4', 'Ndufs2', 'Pex19', 'Opn3', 'Olfm1', 'Zeb2', 'Ttc21b', 'Dlx1', 'Olfr1056', 'C1qtnf4', 'Ano3', 'Rmdn3', 'Exd1', 'Nusap1', 'Casc4', 'Slc20a1', 'Chgb', 'Napb', 'Rbm39', 'Snhg11', 'Rims4', 'Zfos1', 'Syp', 'Pgrmc1', 'Xpnpep2', 'Phf6', 'Dkc1', 'Xist', 'Drp2', 'Bex2', 'Ngfrap1', 'Phf8', 'Ppef1', 'Piga', 'Sec62', 'Ttc14', '4932438A13Rik', 'Ccdc169', 'Fstl5', 'Gria2', 'Tsacc', 'S100a5', 'Gm43062', 'Dennd2c', 'Wdr77', 'Slc16a4', 'Ntng1', 'Plppr4', 'Sep15', 'Gng5', 'Pnisr', 'Gm20878', 'Col15a1', 'Whrn', 'Zmym4', 'Zbtb40', 'Cda', 'Ttc34', 'Mterf1a', 'Reln', 'Lhfpl3', 'Fosl2', 'Epha5', 'Adamts3', 'Pf4', 'Tmem150cos', 'Rabgef1', 'Gnb2', 'Prkar1b', 'Gm43597', 'Met', 'Nap1l5', 'Stambp', 'Pcyox1', 'March8', 'Etnk1', 'Mboat7', 'Chmp2a', 'Zfp658', 'Slc17a7', 'Ube3a', 'Nmb', 'Omp', 'Plekhb1', 'Spon1', 'RP23-44H21.1', 'Grm1', 'Themis', 'Fam229b', 'Serinc1', 'Atp5d', 'Gm15608', 'Syt1', 'Rdh16', 'Wrn', 'Cpe', 'Hmgxb4', 'Nrn1l', 'Calb2', 'Rfwd3', 'Zcchc14', 'Gm26759', 'Gm2974', 'Fhit', 'Synpr', 'Il3ra', 'Cdhr1', 'Ipo4', 'Ppp2r2a', 'Kbtbd3', 'Trpc6', 'Zfp426', 'Slc37a2', '2010007H06Rik', 'Nptn', 'B930082K07Rik', 'Map2k1', 'Rwdd2a', 'Mthfsl', 'Tcta', 'Cck', 'Nktr', 'Rtn4', 'Gabra1', 'Gabrb2', 'Adra1b', 'Rpl26', 'Aurkb', 'Gm11205', 'Crlf3', 'Bzrap1', 'Pdk2', 'Ikzf3', 'Krt222', 'Map3k3', 'Hn1', 'H3f3b', 'Srsf2', 'Tnrc6c', 'Tk1', 'Tha1', 'Elmo1', 'Hist1h3e', 'Nrsn1', 'Nedd9', 'Hk3', 'Zfp808', 'Ankrd55', 'Mrps30', 'Paip1', '4921508M14Rik', 'Tspan13', 'Etv1', 'Rtn1', 'Six1', 'Meg3', 'Amn', 'Rictor', 'Rpl30', 'Oxr1', 'Adcy8', '1700109K24Rik', 'Syngr1', 'Fbln1', 'Selo', '5330439K02Rik', 'Tns2', 'Alg1', 'Gm15738', 'Chrd', 'Gsk3b', 'Tagln3', 'St3gal6', 'Son', 'Sft2d1', 'Slc22a3', 'Ndufb10', 'Gng13', 'Gm28043', 'H2-D1', 'Safb2', 'Wdr43', 'Taf4b', 'Pcdhac2', 'Ablim1', '6430562O15Rik', 'Gm9821', 'Lrrc71'
)


gene2 = c(
  "1110017D15Rik", "1700012B09Rik", "3632451O06Rik", "4430402I18Rik", "Acbd7",
  "Adcyap1", "AI593442", "Ano3", "Apold1", "Atf5", "Atp2b1", "Atp6ap1l",
  "B3galt2", "Cacna2d2", "Cacnb2", "Cbln1", "Cbln4", "Cck",
  "Cd63", "Cdhr1", "Cdkn1b", "Chga", "Chgb", "Chst8",
  "Cited2", "Cnr1", "Coro6", "Cpa6", "Cpne6", "Crim1",
  "Cyth1", "Dkkl1", "Dner", "Dnm3", "Efcab10", "Eml5",
  "Enho", "Eomes", "Etv1", "Fam213b", "Fstl5", "Gabra1",
  "Gm11549", "Gng13", "Gpr22", "Gpsm1", "Grem1", "Gria3",
  "Grin2b", "Grm1", "Gucy1a3", "H3f3b", "Icam5", "Igf1",
  "Il1rapl1", "Kank3", "Kcna2", "Kcnb2", "Kcnip2", "Kitl",
  "L1cam", "Lhfp", "Lingo1", "Lrrc55", "Lrrtm1", "Magi1",
  "Marcksl1", "Mgll", "Ms4a15", "Nap1l5", "Ndnf", "Nfib",
  "Nmb", "Nppa", "Npr1", "Nptx1", "Nptx2", "Nrn1",
  "Nrn1l", "Nrsn1", "Nxph1", "Olfm1", "Omp", "Opn3",
  "Otop2", "Pcp4", "Pde1c", "Pde5a", "Penk", "Pik3r1",
  "Prkca", "Pth2", "Ptprd", "Rab3c", "Reln", "Resp18",
  "Rgs10", "Rgs7", "Riiad1", "Scn2b", "Scube1", "Sdc3",
  "Sema3c", "Serpine2", "Shisa3", "Slc17a7", "Slc1a3", "Slc20a1",
  "Slc6a7", "Sox4", "Sp8", "Spock1", "Spon1", "Spp1",
  "Stoml3", "Sv2b", "Syndig1l", "Syt6", "Syt7", "Tagln3",
  "Thsd7a", "Tmem163", "Tmem176b", "Tpt1", "Trh", "Tspan13",
  "Tspan9", "Tstd1", "Ttc28", "Whrn")


intersect(gene1, gene2)

```



```{r}
alpha_data <- read.csv(here("Spatial_Blade/Results/Alpha_values", "alpha_norp.csv"))

head(alpha_data)



alpha_data <- BLADE_proportion2

if (nrow(alpha_data) == nrow(locations)) {
  # Append columns 'x' and 'y' from 'locations' dataframe to 'alpha_data'
  alpha_data <- cbind(alpha_data, locations[, c('x', 'y')])
} else {
  print("Number of rows in 'alpha_data' and 'locations' do not match!")
}

alpha_correlation <- apply(alpha_data[, c('GC', 'PGC', 'OSNs', 'M.TC')], 1, function(row) {
  cor(row, colMeans(alpha_data[, c('GC', 'PGC', 'OSNs', 'M.TC')], na.rm = TRUE), use = "complete.obs")
})

# Create a dataframe with correlation values and coordinates
correlation_df <- data.frame(Correlation = alpha_correlation, x = locations$x, y = locations$y)

# Plot the correlations against coordinates
ggplot(correlation_df, aes(x = x, y = y, color = Correlation)) +
  geom_point(size =3) +
  scale_color_gradient(low = "blue", high = "red", name = "Correlation") +
  labs(title = "Correlation of Alpha Data with Coordinates",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()




```



```{r}
alpha_data <- read.csv(here("Spatial_Blade/Results", "BLADE_norp_alpha.csv"))
alpha_data <- read.csv(here("Spatial_Blade/Results", "BLADE_norp_alphagroup_1_.csv"))
beta_data <-  read.csv(here("Spatial_Blade/Results/Alpha_values", "beta_norp_group.csv"))

alpha_data <- CARD_proportion2


if (nrow(alpha_data) == nrow(locations)) {
  # Append columns 'x' and 'y' from 'locations' dataframe to 'alpha_data'
  alpha_data <- cbind(alpha_data, locations[, c('x', 'y')])
} else {
  print("Number of rows in 'alpha_data' and 'locations' do not match!")
}

alpha_correlation <- apply(alpha_data[, c('GC', 'PGC', 'OSNs', 'M/TC')], 1, function(row) {
  cor(row, colMeans(alpha_data[, c('GC', 'PGC', 'OSNs', 'M/TC')], na.rm = TRUE), use = "complete.obs")
})

# Create a dataframe with correlation values and coordinates
correlation_df <- data.frame(Correlation = alpha_correlation, x = locations$x, y = locations$y)

# Plot the correlations against coordinates
ggplot(correlation_df, aes(x = x, y = y, color = Correlation)) +
  geom_point(size =3.8) +
  scale_color_gradient(low = "blue", high = "red", name = "Correlation") +
  labs(title = "Correlation of spots with average expression across all spots",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()
```


```{r}
# Calculate correlation coefficient for each cell type
correlation <- apply(beta_data[, -1], 2, function(x) cor(x, alpha_data[, colnames(beta_data[, -1])], use = "complete.obs"))

# Plotting correlation coefficient
cell_types <- colnames(beta_data[, -1])
correlation_df <- data.frame(Cell_Type = cell_types, Correlation = correlation)

# Reshape data for plotting
correlation_df_long <- melt(correlation_df, id.vars = "Cell_Type", variable.name = "Metric")

# Plot correlation coefficient
ggplot(correlation_df_long, aes(x = Cell_Type, y = value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cell Type", y = "Correlation Coefficient", fill = "Metric") +
  ggtitle("Comparison of Alpha and Beta for Each Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



```

```{r}
# Transpose the alpha_data dataframe
alpha_data_transposed <- t(alpha_data[, -1])  # Exclude the first column if it's not needed

# Calculate correlation coefficient between rows
row_correlation <- cor(alpha_data_transposed, use = "complete.obs")

# Plotting correlation coefficient
heatmap(row_correlation, 
        Rowv = NA, Colv = NA,      # Turn off row and column clustering
        col = colorRampPalette(c("blue", "white", "red"))(100),  # Define color scheme
        scale = "none",           # Do not scale rows and columns separately
        xlab = "Row", ylab = "Row",
        main = "Correlation Between Rows in Alpha Data")


```

```{r}

```


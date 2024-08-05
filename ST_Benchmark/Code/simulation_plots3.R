
blade_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Medium', 'Hard'), times = 4),
  Metric = rep(c('GC', 'MTC', 'OSNs', 'PGC'), each = 5),
  Value = c(
    blade_baseline1$Moran_Differences$MoransI_GC, blade_baseline2$Moran_Differences$MoransI_GC, blade_baseline3$Moran_Differences$MoransI_GC, blade_medium$Moran_Differences$MoransI_GC, blade_hard$Moran_Differences$MoransI_GC,
    blade_baseline1$Moran_Differences$MoransI_M.TC, blade_baseline2$Moran_Differences$MoransI_M.TC, blade_baseline3$Moran_Differences$MoransI_M.TC, blade_medium$Moran_Differences$MoransI_M.TC, blade_hard$Moran_Differences$MoransI_M.TC,
    blade_baseline1$Moran_Differences$MoransI_OSNs, blade_baseline2$Moran_Differences$MoransI_OSNs, blade_baseline3$Moran_Differences$MoransI_OSNs, blade_medium$Moran_Differences$MoransI_OSNs, blade_hard$Moran_Differences$MoransI_OSNs,
    blade_baseline1$Moran_Differences$MoransI_PGC, blade_baseline2$Moran_Differences$MoransI_PGC, blade_baseline3$Moran_Differences$MoransI_PGC, blade_medium$Moran_Differences$MoransI_PGC, blade_hard$Moran_Differences$MoransI_PGC
  )
)


bladesp_metrics <- data.frame(
  Condition =rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Medium', 'Hard'), times = 4),
  Metric = rep(c('GC', 'MTC', 'OSNs', 'PGC'), each = 5),
  Value = c(
    bladesp_baseline1$Moran_Differences$MoransI_GC, bladesp_baseline2$Moran_Differences$MoransI_GC, bladesp_baseline3$Moran_Differences$MoransI_GC, bladesp_medium$Moran_Differences$MoransI_GC, bladesp_hard$Moran_Differences$MoransI_GC,
    bladesp_baseline1$Moran_Differences$MoransI_M.TC, bladesp_baseline2$Moran_Differences$MoransI_M.TC, bladesp_baseline3$Moran_Differences$MoransI_M.TC, bladesp_medium$Moran_Differences$MoransI_M.TC, bladesp_hard$Moran_Differences$MoransI_M.TC,
    bladesp_baseline1$Moran_Differences$MoransI_OSNs, bladesp_baseline2$Moran_Differences$MoransI_OSNs, bladesp_baseline3$Moran_Differences$MoransI_OSNs, bladesp_medium$Moran_Differences$MoransI_OSNs, bladesp_hard$Moran_Differences$MoransI_OSNs,
    bladesp_baseline1$Moran_Differences$MoransI_PGC, bladesp_baseline2$Moran_Differences$MoransI_PGC, bladesp_baseline3$Moran_Differences$MoransI_PGC, bladesp_medium$Moran_Differences$MoransI_PGC, bladesp_hard$Moran_Differences$MoransI_PGC
  )
)




card_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Medium', 'Hard'), times = 4),
  Metric = rep(c('GC', 'MTC', 'OSNs', 'PGC'), each = 5),
  Value = c(
    card_baseline1$Moran_Differences$MoransI_GC, card_baseline2$Moran_Differences$MoransI_GC, card_baseline3$Moran_Differences$MoransI_GC, card_medium$Moran_Differences$MoransI_GC, card_hard$Moran_Differences$MoransI_GC,
    card_baseline1$Moran_Differences$MoransI_M.TC, card_baseline2$Moran_Differences$MoransI_M.TC, card_baseline3$Moran_Differences$MoransI_M.TC, card_medium$Moran_Differences$MoransI_M.TC, card_hard$Moran_Differences$MoransI_M.TC,
    card_baseline1$Moran_Differences$MoransI_OSNs, card_baseline2$Moran_Differences$MoransI_OSNs, card_baseline3$Moran_Differences$MoransI_OSNs, card_medium$Moran_Differences$MoransI_OSNs, card_hard$Moran_Differences$MoransI_OSNs,
    card_baseline1$Moran_Differences$MoransI_PGC, card_baseline2$Moran_Differences$MoransI_PGC, card_baseline3$Moran_Differences$MoransI_PGC, card_medium$Moran_Differences$MoransI_PGC, card_hard$Moran_Differences$MoransI_PGC
  )
)




music_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Medium', 'Hard'), times = 4),
  Metric = rep(c('GC', 'MTC', 'OSNs', 'PGC'), each = 5),
  Value = c(
    music_baseline1$Moran_Differences$MoransI_GC, music_baseline2$Moran_Differences$MoransI_GC, music_baseline3$Moran_Differences$MoransI_GC, music_medium$Moran_Differences$MoransI_GC, music_hard$Moran_Differences$MoransI_GC,
    music_baseline1$Moran_Differences$MoransI_M.TC, music_baseline2$Moran_Differences$MoransI_M.TC, music_baseline3$Moran_Differences$MoransI_M.TC, music_medium$Moran_Differences$MoransI_M.TC, music_hard$Moran_Differences$MoransI_M.TC,
    music_baseline1$Moran_Differences$MoransI_OSNs, music_baseline2$Moran_Differences$MoransI_OSNs, music_baseline3$Moran_Differences$MoransI_OSNs, music_medium$Moran_Differences$MoransI_OSNs, music_hard$Moran_Differences$MoransI_OSNs,
    music_baseline1$Moran_Differences$MoransI_PGC, music_baseline2$Moran_Differences$MoransI_PGC, music_baseline3$Moran_Differences$MoransI_PGC, music_medium$Moran_Differences$MoransI_PGC, music_hard$Moran_Differences$MoransI_PGC
  )
)




rctd_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Medium', 'Hard'), times = 4),
  Metric = rep(c('GC', 'MTC', 'OSNs', 'PGC'), each = 5),
  Value = c(
    rctd_baseline1$Moran_Differences$MoransI_GC, rctd_baseline2$Moran_Differences$MoransI_GC, rctd_baseline3$Moran_Differences$MoransI_GC, rctd_medium$Moran_Differences$MoransI_GC, rctd_hard$Moran_Differences$MoransI_GC,
    rctd_baseline1$Moran_Differences$MoransI_M.TC, rctd_baseline2$Moran_Differences$MoransI_M.TC, rctd_baseline3$Moran_Differences$MoransI_M.TC, rctd_medium$Moran_Differences$MoransI_M.TC, rctd_hard$Moran_Differences$MoransI_M.TC,
    rctd_baseline1$Moran_Differences$MoransI_OSNs, rctd_baseline2$Moran_Differences$MoransI_OSNs, rctd_baseline3$Moran_Differences$MoransI_OSNs, rctd_medium$Moran_Differences$MoransI_OSNs, rctd_hard$Moran_Differences$MoransI_OSNs,
    rctd_baseline1$Moran_Differences$MoransI_PGC, rctd_baseline2$Moran_Differences$MoransI_PGC, rctd_baseline3$Moran_Differences$MoransI_PGC, rctd_medium$Moran_Differences$MoransI_PGC, rctd_hard$Moran_Differences$MoransI_PGC
  )
)




cibersort_metrics <- data.frame(
  Condition =rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Medium', 'Hard'), times = 4),
  Metric = rep(c('GC', 'MTC', 'OSNs', 'PGC'), each = 5),
  Value = c(
    cibersort_baseline1$Moran_Differences$MoransI_GC, cibersort_baseline2$Moran_Differences$MoransI_GC, cibersort_baseline3$Moran_Differences$MoransI_GC, cibersort_medium$Moran_Differences$MoransI_GC, cibersort_hard$Moran_Differences$MoransI_GC,
    cibersort_baseline1$Moran_Differences$MoransI_M.TC, cibersort_baseline2$Moran_Differences$MoransI_M.TC, cibersort_baseline3$Moran_Differences$MoransI_M.TC, cibersort_medium$Moran_Differences$MoransI_M.TC, cibersort_hard$Moran_Differences$MoransI_M.TC,
    cibersort_baseline1$Moran_Differences$MoransI_OSNs, cibersort_baseline2$Moran_Differences$MoransI_OSNs, cibersort_baseline3$Moran_Differences$MoransI_OSNs, cibersort_medium$Moran_Differences$MoransI_OSNs, cibersort_hard$Moran_Differences$MoransI_OSNs,
    cibersort_baseline1$Moran_Differences$MoransI_PGC, cibersort_baseline2$Moran_Differences$MoransI_PGC, cibersort_baseline3$Moran_Differences$MoransI_PGC, cibersort_medium$Moran_Differences$MoransI_PGC, cibersort_hard$Moran_Differences$MoransI_PGC
  )
)

my_palette <- c("#FF0000", "#FFFFFF","#0000FF")  # Red to White to Blue

# Function to create a heatmap for each method

set_condition_levels <- function(df) {
  df$Condition <- factor(df$Condition, levels = c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Medium', 'Hard'))
  return(df)
}

# Apply the function to each data frame
blade_metrics <- set_condition_levels(blade_metrics)
bladesp_metrics <- set_condition_levels(bladesp_metrics)
card_metrics <- set_condition_levels(card_metrics)
music_metrics <- set_condition_levels(music_metrics)
rctd_metrics <- set_condition_levels(rctd_metrics)
cibersort_metrics <- set_condition_levels(cibersort_metrics)
plot_heatmap <- function(method_data, method_name, show_x_labels = TRUE) {
  plot <- ggplot(method_data, aes(x = Metric, y = Condition, fill = Value)) +
    geom_tile() +
    geom_text(aes(label = round(Value, 2)), size = 4) +  # Add metric values to the tiles
    scale_fill_gradientn(colors = my_palette, limits = c(-1, 1)) +  # Adjust color scale limits
    labs(title = method_name, x = NULL, y = NULL) +  # Add method name as title
    theme_minimal() +
    theme(axis.text = element_text(size = 12, face = "bold"),   # Set axis text size for both x and y axes
          axis.title = element_blank(),   # Remove axis titles
          panel.grid = element_blank(),   # Remove gridlines
          legend.position = "none",       # Remove the legend
          plot.title = element_text(size = 15, face = "bold"))
  
  if (!show_x_labels) {
    plot <- plot + theme(axis.text.y = element_blank(), axis.ticks.x = element_blank())
  }
  
  return(plot)
}

# List of data frames for each method
methods <- list(
  BLADE = blade_metrics,
  "SPATIAL BLADE" = bladesp_metrics,
  CARD = card_metrics,
  MuSiC = music_metrics,
  RCTD = rctd_metrics,
  CIBERSORT = cibersort_metrics
)

# Create a list to store heatmaps
heatmaps <- list()

# Create heatmaps for each method and arrange them
for (i in seq_along(methods)) {
  method_name <- names(methods)[i]
  method_data <- methods[[i]]
  show_x_labels <- !(method_name %in% c("SPATIAL BLADE", "CARD", "RCTD", "CIBERSORT"))
  heatmap <- plot_heatmap(method_data, method_name, show_x_labels)
  heatmaps[[i]] <- heatmap
}

# Arrange the heatmaps in a grid
grid.arrange(
  grobs = heatmaps, 
  ncol = 3)
)







blade_metrics_rmse <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 1),
  Metric = rep(c('Spot-Wise RMSE'), each = 5),
  Value = c(
    blade_baseline1$error_results$result$rowwise_rmse, blade_baseline2$error_results$result$rowwise_rmse, blade_baseline3$error_results$result$rowwise_rmse, blade_medium$error_results$result$rowwise_rmse
    , blade_hard$error_results$result$rowwise_rmse
  )
)



bladesp_metrics_rmse <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 1),
  Metric = rep(c('Spot-Wise RMSE'), each = 5),
  Value = c(
    bladesp_baseline1$error_results$result$rowwise_rmse, bladesp_baseline2$error_results$result$rowwise_rmse, bladesp_baseline3$error_results$result$rowwise_rmse, bladesp_medium$error_results$result$rowwise_rmse
    , bladesp_hard$error_results$result$rowwise_rmse
  )
)



card_metrics_rmse <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 1),
  Metric = rep(c('Spot-Wise RMSE'), each = 5),
  Value = c(
    card_baseline1$error_results$result$rowwise_rmse, card_baseline2$error_results$result$rowwise_rmse, card_baseline3$error_results$result$rowwise_rmse, card_medium$error_results$result$rowwise_rmse
    , card_hard$error_results$result$rowwise_rmse
  )
)


music_metrics_rmse <- data.frame(
  Condition =  rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 1),
  Metric = rep(c('Spot-Wise RMSE'), each = 5),
  Value = c(
    music_baseline1$error_results$result$rowwise_rmse, music_baseline2$error_results$result$rowwise_rmse, music_baseline3$error_results$result$rowwise_rmse, music_medium$error_results$result$rowwise_rmse
    , music_hard$error_results$result$rowwise_rmse
  )
)



cibersort_metrics_rmse <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 1),
  Metric = rep(c('Spot-Wise RMSE'), each = 5),
  Value = c(
    cibersort_baseline1$error_results$result$rowwise_rmse, cibersort_baseline2$error_results$result$rowwise_rmse, cibersort_baseline3$error_results$result$rowwise_rmse, cibersort_medium$error_results$result$rowwise_rmse
    , cibersort_hard$error_results$result$rowwise_rmse
  )
)



rctd_metrics_rmse <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 1),
  Metric = rep(c('Spot-Wise RMSE'), each = 5),
  Value = c(
    rctd_baseline1$error_results$result$rowwise_rmse, rctd_baseline2$error_results$result$rowwise_rmse, rctd_baseline3$error_results$result$rowwise_rmse, rctd_medium$error_results$result$rowwise_rmse
    , rctd_hard$error_results$result$rowwise_rmse
  )
)


plot_violin <- function(method_data, method_name) {
  ggplot(method_data, aes(x = Condition, y = Value, fill = Condition)) +
    geom_violin(trim = FALSE) +  # Violin plot to show data distribution
    labs(title = method_name, x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),  # Remove x-axis labels
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 12, face = "bold"),  # Bold and bigger text for y-axis
          legend.position = "none",
          plot.title = element_text(size = 14, face = "bold"),  # Bold and bigger text for title
          plot.subtitle = element_text(size = 12, face = "bold"))  # Bold and bigger text for subtitle
}

# Combine all metrics data into one data frame with an additional column for the method
combine_metrics <- function(df, method_name) {
  df$Method <- method_name
  df$Condition <- factor(df$Condition, levels = c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'))
  return(df)
}

all_metrics_rmse <- rbind(
  combine_metrics(blade_metrics_rmse, "BLADE"),
  combine_metrics(card_metrics_rmse, "CARD"),
  combine_metrics(music_metrics_rmse, "MuSiC"),
  combine_metrics(bladesp_metrics_rmse, "Spatial BLADE"),
  combine_metrics(cibersort_metrics_rmse, "CIBERSORT"),
  combine_metrics(rctd_metrics_rmse, "RCTD")
)





# Define custom colors for the methods
custom_colors <- c(
  "BLADE" = "#1f77b4",       # Blue
  "SPATIAL BLADE" = "#aec7e8",# Light Blue
  "CARD" = "#ff7f0f",        # Orange
  "MuSiC" = "#2ca02c",       # Green
  "CIBERSORT" = "#d62728",   # Red
  "RCTD" = "#9467bd"         # Purple
)

combined_plot <- ggplot(all_metrics_rmse, aes(x = Method, y = Value, fill = Method)) +
  geom_violin(trim = FALSE) +  # Violin plot to show data distribution
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", size = 0.5) +  # Mean line
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), vjust = -0.5, color = "black") +  # Mean text
  facet_grid(Condition ~ Method, scales = "free") +  # Create grid of plots
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text for better readability
    axis.title.y = element_blank(),
    legend.position = "left",  # Position the legend on the left
    plot.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 15, face = "bold")
  ) + 
  labs(title = "Spot-Wise RMSE for Each Method Across Scenarios", y = "Spot-Wise RMSE") +  # Adjust title and labels
  coord_cartesian(ylim = c(0, 1))  # Set y-axis limits to 0-1 for all plots



all_metrics_rmse$Method <- factor(all_metrics_rmse$Method, levels = c("BLADE","CARD","CIBERSORT","MuSiC", "RCTD", "Spatial BLADE"))

# Filter data for Scenario 1
scenario_1_data <- all_metrics_rmse %>% filter(Condition == "Scenario 1")

# Create the violin plot for Scenario 1
violin_plot_scenario_1 <- ggplot(scenario_1_data, aes(x = Method, y = Value, fill = Method)) +
  geom_violin(trim = FALSE) +  # Violin plot to show data distribution
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", size = 0.5) +  # Mean line
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), vjust = -1.5, color = "black") +  # Mean text
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # Rotate x-axis text for better readability
    axis.title.y = element_text(size = 10),  # Y-axis title size
    legend.position = "none",  # Remove legend
    plot.title = element_text(size = 10, face = "bold")
  ) +
  labs(title = "Spot-Wise RMSE for Different Methods in Scenario 1", y = "Spot-Wise RMSE")  # Adjust title and labels


# Print the violin plot for Scenario 1
print(violin_plot_scenario_1)

final_plot <- (ari_plot | violin_plot_scenario_1) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 15, face = "bold"))

# Print the final plot
print(final_plot)


# Print the combined plot
print(combined_plot)

boxplot_rmse <- ggplot(all_metrics_rmse, aes(x = Condition, y = Value, fill = Method)) +
  geom_boxplot() +  # Boxplot to show RMSE distribution
  facet_wrap(~ Method, scales = "free_x") +  # Facet by method
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis text
    axis.title.y = element_text(size = 16),  # Y-axis title size
    axis.text.y = element_text(size = 14),   # Y-axis text size
    legend.position = "none",  # Hide legend
    plot.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 15, face = "bold")
  ) +
  labs(title = "Spot-Wise RMSE for Each Method Across Scenarios", y = "Spot-Wise RMSE", x = "Scenario")

# Print the boxplot
print(boxplot_rmse)


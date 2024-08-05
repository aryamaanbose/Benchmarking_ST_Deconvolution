

blade_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 3),
  Metric = rep(c('ARI', 'ARI2', 'Overall RMSE'), each = 5),
  Value = c(
    blade_baseline1$ARI_and_Purity$ARI, blade_baseline2$ARI_and_Purity$ARI, blade_baseline3$ARI_and_Purity$ARI, blade_medium$ARI_and_Purity$ARI, blade_hard$ARI_and_Purity$ARI,
    blade_baseline1$ARI_and_Purity$ARI2, blade_baseline2$ARI_and_Purity$ARI2, blade_baseline3$ARI_and_Purity$ARI2, blade_medium$ARI_and_Purity$ARI2, blade_hard$ARI_and_Purity$ARI2,
    blade_baseline1$error_results$result$overall_rmse, blade_baseline2$error_results$result$overall_rmse, blade_baseline3$error_results$result$overall_rmse, blade_medium$error_results$result$overall_rmse,blade_hard$error_results$result$overall_rmse
  )
)




bladesp_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 3),
  Metric = rep(c('ARI', 'ARI2', 'Overall RMSE'), each = 5),
  Value = c(
    bladesp_baseline1$ARI_and_Purity$ARI, bladesp_baseline2$ARI_and_Purity$ARI, bladesp_baseline3$ARI_and_Purity$ARI, bladesp_medium$ARI_and_Purity$ARI, bladesp_hard$ARI_and_Purity$ARI,
    bladesp_baseline1$ARI_and_Purity$ARI2, bladesp_baseline2$ARI_and_Purity$ARI2, bladesp_baseline3$ARI_and_Purity$ARI2, bladesp_medium$ARI_and_Purity$ARI2, bladesp_hard$ARI_and_Purity$ARI2,
    bladesp_baseline1$error_results$result$overall_rmse, bladesp_baseline2$error_results$result$overall_rmse, bladesp_baseline3$error_results$result$overall_rmse, bladesp_medium$error_results$result$overall_rmse,bladesp_hard$error_results$result$overall_rmse
  )
)




card_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 3),
  Metric = rep(c('ARI', 'ARI2', 'Overall RMSE'), each = 5),
  Value = c(
    card_baseline1$ARI_and_Purity$ARI, card_baseline2$ARI_and_Purity$ARI, card_baseline3$ARI_and_Purity$ARI, card_medium$ARI_and_Purity$ARI, card_hard$ARI_and_Purity$ARI,
    card_baseline1$ARI_and_Purity$ARI2, card_baseline2$ARI_and_Purity$ARI2, card_baseline3$ARI_and_Purity$ARI2, card_medium$ARI_and_Purity$ARI2, card_hard$ARI_and_Purity$ARI2,
    card_baseline1$error_results$result$overall_rmse, card_baseline2$error_results$result$overall_rmse, card_baseline3$error_results$result$overall_rmse, card_medium$error_results$result$overall_rmse,card_hard$error_results$result$overall_rmse
  )
)



music_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 3),
  Metric = rep(c('ARI', 'ARI2', 'Overall RMSE'), each = 5),
  Value = c(
    music_baseline1$ARI_and_Purity$ARI, music_baseline2$ARI_and_Purity$ARI, music_baseline3$ARI_and_Purity$ARI, music_medium$ARI_and_Purity$ARI, music_hard$ARI_and_Purity$ARI,
    music_baseline1$ARI_and_Purity$ARI2, music_baseline2$ARI_and_Purity$ARI2, music_baseline3$ARI_and_Purity$ARI2, music_medium$ARI_and_Purity$ARI2, music_hard$ARI_and_Purity$ARI2,
    music_baseline1$error_results$result$overall_rmse, music_baseline2$error_results$result$overall_rmse, music_baseline3$error_results$result$overall_rmse, music_medium$error_results$result$overall_rmse,music_hard$error_results$result$overall_rmse
  )
)





cibersort_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 3),
  Metric = rep(c('ARI', 'ARI2', 'Overall RMSE'), each = 5),
  Value = c(
    cibersort_baseline1$ARI_and_Purity$ARI, cibersort_baseline2$ARI_and_Purity$ARI, cibersort_baseline3$ARI_and_Purity$ARI, cibersort_medium$ARI_and_Purity$ARI, cibersort_hard$ARI_and_Purity$ARI,
    cibersort_baseline1$ARI_and_Purity$ARI2, cibersort_baseline2$ARI_and_Purity$ARI2, cibersort_baseline3$ARI_and_Purity$ARI2, cibersort_medium$ARI_and_Purity$ARI2, cibersort_hard$ARI_and_Purity$ARI2,
    cibersort_baseline1$error_results$result$overall_rmse, cibersort_baseline2$error_results$result$overall_rmse, cibersort_baseline3$error_results$result$overall_rmse, cibersort_medium$error_results$result$overall_rmse,cibersort_hard$error_results$result$overall_rmse
  )
)





rctd_metrics <- data.frame(
  Condition = rep(c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5'), times = 3),
  Metric = rep(c('ARI', 'ARI2', 'Overall RMSE'), each = 5),
  Value = c(
    rctd_baseline1$ARI_and_Purity$ARI, rctd_baseline2$ARI_and_Purity$ARI, rctd_baseline3$ARI_and_Purity$ARI, rctd_medium$ARI_and_Purity$ARI, rctd_hard$ARI_and_Purity$ARI,
    rctd_baseline1$ARI_and_Purity$ARI2, rctd_baseline2$ARI_and_Purity$ARI2, rctd_baseline3$ARI_and_Purity$ARI2, rctd_medium$ARI_and_Purity$ARI2, rctd_hard$ARI_and_Purity$ARI2,
    rctd_baseline1$error_results$result$overall_rmse, rctd_baseline2$error_results$result$overall_rmse, rctd_baseline3$error_results$result$overall_rmse, rctd_medium$error_results$result$overall_rmse,rctd_hard$error_results$result$overall_rmse
  )
)


desired_order <- c('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5')

# Set factor levels for each metrics data frame
set_order <- function(metrics_df) {
  metrics_df$Condition <- factor(metrics_df$Condition, levels = desired_order)
  return(metrics_df)
}

blade_metrics <- set_order(blade_metrics)
card_metrics <- set_order(card_metrics)
music_metrics <- set_order(music_metrics)
cibersort_metrics <- set_order(cibersort_metrics)
rctd_metrics <- set_order(rctd_metrics)
bladesp_metrics <- set_order(bladesp_metrics)


# Function to create individual metric plots
create_metric_plot <- function(metrics_df, metric_name, remove_x_label = FALSE, remove_y_label = TRUE, add_title = FALSE, title = "") {
  p <- ggplot(metrics_df, aes(x = Condition, y = Value, fill = Value)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_text(aes(label = sprintf("%.2f", Value)), vjust = ifelse(metrics_df$Value < 0.4, -0.5, 1.5), size = 3.5, position = position_dodge(width = 0.8)) +
    scale_fill_gradient(low = "orange",high = "red", guide = "none") +  # Remove the color scale legend
    scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
    labs(x = NULL, y = if (remove_y_label) NULL else metric_name) +
    theme_minimal() +
    theme(
      axis.text.x = if (remove_x_label) element_blank() else element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_blank(),  # Hide y-axis labels
      axis.ticks.y = element_blank(),  # Hide y-axis ticks
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = if (remove_y_label) element_blank() else element_text(size = 10),
      plot.title = element_blank()
    )
  
  if (add_title) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}




create_combined_plot <- function(blade_metrics, bladesp_metrics, card_metrics, music_metrics, cibersort_metrics, rctd_metrics) {
  # Define plot groups with titles using plot_annotation correctly
  blade_group <- plot_grid(
    create_metric_plot(blade_metrics %>% filter(Metric == "ARI"), "ARI", TRUE, TRUE, TRUE, "Blade ARI"),
    create_metric_plot(blade_metrics %>% filter(Metric == "ARI2"), "ARI2", FALSE, TRUE, TRUE, "Blade ARI2"),
    ncol = 1, labels = c("A", "B")
  ) / plot_annotation(title = "BLADE", theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  bladesp_group <- plot_grid(
    create_metric_plot(bladesp_metrics %>% filter(Metric == "ARI"), "ARI", TRUE, TRUE, TRUE, "Bladesp ARI"),
    create_metric_plot(bladesp_metrics %>% filter(Metric == "ARI2"), "ARI2", FALSE, TRUE, TRUE, "Bladesp ARI2"),
    ncol = 1
  ) / plot_annotation(title = "SPATIAL BLADE", theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  card_group <- plot_grid(
    create_metric_plot(card_metrics %>% filter(Metric == "ARI"), "ARI", TRUE, TRUE, TRUE, "Card ARI"),
    create_metric_plot(card_metrics %>% filter(Metric == "ARI2"), "ARI2", FALSE, TRUE, TRUE, "Card ARI2"),
    ncol = 1
  ) / plot_annotation(title = "CARD", theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  music_group <- plot_grid(
    create_metric_plot(music_metrics %>% filter(Metric == "ARI"), "ARI", TRUE, TRUE, TRUE, "Music ARI"),
    create_metric_plot(music_metrics %>% filter(Metric == "ARI2"), "ARI2", FALSE, TRUE, TRUE, "Music ARI2"),
    ncol = 1, labels =c("A", "B")
  ) / plot_annotation(title = "MUSIC", theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  cibersort_group <- plot_grid(
    create_metric_plot(cibersort_metrics %>% filter(Metric == "ARI"), "ARI", TRUE, TRUE, TRUE, "Cibersort ARI"),
    create_metric_plot(cibersort_metrics %>% filter(Metric == "ARI2"), "ARI2", FALSE, TRUE, TRUE, "Cibersort ARI2"),
    ncol = 1
  ) / plot_annotation(title = "CIBERSORT", theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  rctd_group <- plot_grid(
    create_metric_plot(rctd_metrics %>% filter(Metric == "ARI"), "ARI", TRUE, TRUE, TRUE, "RCTD ARI"),
    create_metric_plot(rctd_metrics %>% filter(Metric == "ARI2"), "ARI2", FALSE, TRUE, TRUE, "RCTD ARI2"),
    ncol = 1
  ) / plot_annotation(title = "RCTD", theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  # Combine the groups into one plot
  combined_plot <- plot_grid(
    blade_group, bladesp_group, card_group, music_group, cibersort_group, rctd_group,
    ncol = 3
  )
  
  return(combined_plot)
}

# Create and print the combined plot
combined_plot <- create_combined_plot(blade_metrics,bladesp_metrics, card_metrics, music_metrics, cibersort_metrics, rctd_metrics)
print(combined_plot)


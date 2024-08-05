genes <- c(200, 400, 600)
blade_ari <- c(0.25, 0.19, 0.30)
music_ari <- c(0.28, 0.26, 0.30)
card_ari <- c(0.27, 0.28, 0.27)

# Create a data frame
data <- data.frame(
  Genes = rep(genes, 3),
  ARI = c(blade_ari, music_ari, card_ari),
  Method = rep(c("BLADE", "MuSiC", "CARD"), each = 3)
)

# Plot
ggplot(data, aes(x = Genes, y = ARI, color = Method)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  labs(title = "Adjusted Rand Index (ARI) vs Number of Genes",
       x = "Number of Genes",
       y = "Adjusted Rand Index (ARI)") +
  theme_minimal() +
  theme(legend.title = element_blank())




# Create data frame
data <- data.frame(
  Method = rep(c("BLADE", "CARD", "MuSiC"), each = 4),
  DEG_Type = rep(c("DEG 200", "DEG No Ribosomal", "DEG No Significantly Different", "DEG Upregulated Genes"), times = 3),
  Value = c(
    0.25, 0.36, 0.37, 0.35,  # BLADE
    0.27, 0.27, 0.26, 0.02,  # CARD
    0.28, 0.21, 0.04, 0.29   # MuSiC
  )
)

bar_plot <- ggplot(data, aes(x = DEG_Type, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +  # Bar plot with dodge position
  geom_text(aes(label = Value), vjust = -0.3, color = "black", 
            position = position_dodge(width = 0.9), size = 3.5) +  # Add values on top of bars
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis text
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 14),  # Y-axis title size
    legend.title = element_text(size = 14),  # Legend title size
    legend.text = element_text(size = 12),   # Legend text size
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(title = "DEG Values for BLADE, CARD, and MuSiC", y = "Value", fill = "Method")  # Adjust title and labels

# Print the bar plot
print(bar_plot)






all_metrics <- bind_rows(
  blade_metrics %>% mutate(Method = 'BLADE'),
  card_metrics %>% mutate(Method = 'CARD'),
  music_metrics %>% mutate(Method = 'MuSiC'),
  bladesp_metrics %>% mutate(Method = 'Spatial BLADE'),
  cibersort_metrics %>% mutate(Method = 'CIBERSORT'),
  rctd_metrics %>% mutate(Method = 'RCTD')
)



all_metrics_rmse$Method <- factor(all_metrics_rmse$Method, levels = c("BLADE", "SPATIAL BLADE", "CARD", "MuSiC", "CIBERSORT", "RCTD"))

# Filter data for ARI
ari_data <- all_metrics %>% filter(Metric == 'ARI')

# Filter data for ARI2
ari2_data <- all_metrics %>% filter(Metric == 'ARI2')

# Plot for ARI
ari_plot <- ggplot(ari_data, aes(x = Condition, y = Value, color = Method)) +
  geom_point(size = 3) +
  geom_line(aes(group = Method)) +
  labs(title = "ARI for Estimated Dominant Cell Types in Different Simulated Scenarios", x = "Scenarios", y = "ARI") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 10),  # Adjust the size here
    axis.title.x = element_text(size = 10), # Adjust the title size
    axis.title.y = element_text(size = 10), # Adjust the y-axis title size
    axis.text.y = element_text(size = 10)# Adjust the y-axis text size
  )

# Plot for ARI2
ari2_plot <- ggplot(ari2_data, aes(x = Condition, y = Value, color = Method)) +
  geom_point(size = 3) +
  geom_line(aes(group = Method)) +
  labs(title = "ARI for Secondary Dominant Cell Types vs Scenario", x = "Scenarios", y = "ARI") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),  # Adjust the size here
    axis.title.x = element_text(size = 16), # Adjust the title size
    axis.title.y = element_text(size = 16), # Adjust the y-axis title size
    axis.text.y = element_text(size = 14)   # Adjust the y-axis text size
  )

# Print the plots
print(ari_plot)
print(ari2_plot)

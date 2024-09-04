
blade_ARI <- data.frame(
  Method = c("Baseline", "30 %", "70 %", "97 %"),
  ARI = c(results_BLADE1$ARI, results_30_SP$ARI, results_50_SP$ARI, results_90_SP$ARI)
)

# Convert Method to factor with specified levels to maintain order
blade_ARI$Method <- factor(blade_ARI$Method, levels = unique(blade_ARI$Method))


# Plot ARI with connected points and formatted labels
ARI <- ggplot(blade_ARI, aes(x = Method, y = ARI)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1), size = 1) +
  geom_text(aes(label = sprintf("%.3f", ARI)), vjust = -1, size = 5) +  # Format ARI to 3 decimal places, increase text size
  theme_minimal() +
  labs(title = "ARI for spatialBLADE Tests",
       x = "Results",
       y = "Adjusted Rand Index (ARI)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    plot.title = element_text( size = 18),  # Increase plot title size
    legend.position = "none"
  )


ARI + scatter_plot

additional_points <- data.frame(
  Method = c("30 %", "70 %", "97 %"),
  ARI = c(0.31, 0.65, 0.83)
)

# Convert Method to factor with specified levels to maintain order
additional_points$Method <- factor(additional_points$Method, levels = c("30 %", "70 %", "97 %"))

ggplot(additional_points, aes(x = Method, y = ARI, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ARI), vjust = -0.5) + # Adjust vjust for better label positioning
  scale_fill_manual(values = c("grey", "grey", "grey")) + # Manual colors to keep consistency
  theme_minimal() +
  labs(title = "PCC of Estimated Proportions with Expectation Matrix",
       x = "Expectation Matrix",
       y = "Pearson Correlation Coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none")




blade_metrics <- data.frame(
  Method = c("Baseline", "30%", "70%", "97%"),
  GC = c(results_BLADE1$GC, results_30_SP$GC, results_50_SP$GC, results_90_SP$GC),
  PGC = c(results_BLADE1$PGC, results_30_SP$PGC, results_50_SP$PGC, results_90_SP$PGC),
  MTC = c(results_BLADE1$MTC, results_30_SP$MTC, results_50_SP$MTC, results_90_SP$MTC),
  OSN = c(results_BLADE1$OSN, results_30_SP$OSN, results_50_SP$OSN, results_90_SP$OSN)
)



blade_metrics_melted <- melt(blade_metrics, id.vars = "Method", variable.name = "Metric", value.name = "Value")

# Ensure the Method column follows the specified order
blade_metrics_melted$Method <- factor(blade_metrics_melted$Method, levels = c("Baseline", "30%", "70%", "97%"))

# Create a boxplot for each method with all metrics combined
Corr <- ggplot(blade_metrics_melted, aes(x = Method, y = Value, fill = Metric)) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = "white") +
  geom_jitter(position = position_jitter(0.2), size = 2, aes(color = Metric)) +
  scale_fill_manual(values = c("GC" = "lightblue", "PGC" = "#ed78ff", "MTC" = "orange", "OSN" = "red")) +
  scale_color_manual(values = c("GC" = "lightblue", "PGC" = "#ed78ff", "MTC" = "orange", "OSN" = "red")) +
  theme_minimal() +
  labs(title = "PCC of Estimated Proportions with Cell Type Marker Gene Expression",
       x = "Results",
       y = "Pearson Correlation Coefficient") +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        legend.position = "right")



scatter_plot <- ggplot(data = layer_manual_MOB, aes(x = x, y = y, color = celltypes)) +
  geom_point(aes(colour = celltypes), size = 3) +
  labs(title = "Ground Truth Annotation", color = "Cell Types") +
  scale_color_manual(values = c("GC" = "lightblue", "M.TC" = "orange", "OSNs" = "red", "PGC" = "#9999FF")) +
  theme(
    plot.title = element_text(size = 20),  # Increase title font siz
    legend.text = element_text(size = 17),  # Increase legend text font size
    legend.title = element_text(size = 20)  # Increase legend title font size
  )
# Modify individual plots with increased title font size
p1 <- results_BLADE1$plot + theme_minimal() + labs(title = "Baseline, ARI: 0.249") + theme(
  plot.title = element_text(hjust = 0.5, size = 20),  # Center and increase title font size
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none"
)
p2 <- results_30_SP$plot + theme_minimal() + labs(title = "30%, ARI: 0.161",) + theme(
  plot.title = element_text(hjust = 0.5, size = 20),  # Center and increase title font size
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none"
)
p3 <- results_50_SP$plot + theme_minimal() + labs(title = "70%, ARI: 0.266") + theme(
  plot.title = element_text(hjust = 0.5, size = 20),  # Center and increase title font size
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none"
)

p4 <- results_90_SP$plot + theme_minimal() + labs(title = "97%, ARI: 0.428") + theme(
  plot.title = element_text(hjust = 0.5, size = 20),  # Center and increase title font size
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.text = element_text(size = 17),  # Increase legend text size
  legend.title = element_text(size = 20)  # Increase legend title size
)

# Combine the individual plots
combined_plot <- (p1 | p2 | p3 | p4) +
  plot_annotation(
    title = "Inferred Dominant Cell Types",
    theme = theme(plot.title = element_text(size = 20))  # Increase title font size
  )

scatter_plot + combined_plot 

# Create a label for the combined plot
combined_plot_labeled <- combined_plot 

# Combine the ARI and Corr plots with the labeled combined plot
final_plot <- ((ARI | Corr) / combined_plot_labeled) + 
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag = element_text(size = 15, face = "bold"),  # Increase tag size
    plot.title = element_text(size = 18, face = "bold")  # Increase font size of titles in the combined plot
  )

# Print the final plot
print(final_plot)

# Print the combined plot
print(combined_plot)

scatter_plot + combined

results_90_SP$prop_matrix

results_50_SP$prop_matrix

results_30_SP$prop_matrix



read_and_set_rownames <- function(filepath) {
  data <- read_csv(filepath)
  rownames(data) <- data[[1]]
  data <- data[,-1]
  return(data)
}


em_50 <- read_and_set_rownames(here("Data/Processed_BLADE", "Expected70.csv"))
em_90 <- read_and_set_rownames(here("Data/Processed_BLADE", "Expected97.csv"))
em_30 <- read_and_set_rownames(here("Data/Processed_BLADE", "Expected30.csv"))


results_30_vec <- as.vector(as.matrix(results_30_SP$prop_matrix))
results_50_vec <- as.vector(as.matrix(results_50_SP$prop_matrix))
results_90_vec <- as.vector(as.matrix(results_90_SP$prop_matrix))

em_30_vec <- as.vector(as.matrix(em_30))
em_50_vec <- as.vector(as.matrix(em_50))
em_90_vec <- as.vector(as.matrix(em_90))

# Calculate overall Pearson correlations
cor_30 <- cor(results_30_vec, em_30_vec, method = "pearson")
cor_50 <- cor(results_50_vec, em_50_vec, method = "pearson")
cor_90 <- cor(results_90_vec, em_90_vec, method = "pearson")

# Print correlations
correlations <- data.frame(
  Method = c("30% Confidence", "50% Confidence", "90% Confidence"),
  Correlation = c(cor_30, cor_50, cor_90)
)
print(correlations)

# Function to create scatter plot for overall PCC
create_scatter_plot <- function(data1, data2, method_name) {
  ggplot(data = data.frame(x = data1, y = data2), aes(x = x, y = y)) +
    geom_point(color = "blue") +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    ggtitle(paste("Scatterplot:", method_name, "(Pearson =", round(cor(data1, data2, use = "complete.obs"), 2), ")")) +
    xlab(paste("Results", method_name)) +
    ylab(paste("Expected", method_name))
}

# Create scatterplots to visualize the overall correlations
plot_30 <- create_scatter_plot(results_30_vec, em_30_vec, "30% Confidence")
plot_50 <- create_scatter_plot(results_50_vec, em_50_vec, "70% Confidence")
plot_90 <- create_scatter_plot(results_90_vec, em_90_vec, "90% Confidence")

# Print the scatterplots
print(plot_30)
print(plot_50)
print(plot_90)














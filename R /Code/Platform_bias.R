# Assuming spe and seurat_object1 are already loaded in your environment
spe <- subset(spe, features = rownames(seurat_object1))
# Calculate the average expression for both datasets
avg_exp_spe <- rowMeans(spe@assays$spatial$data)
avg_exp_seurat_object1 <- rowMeans(seurat_object1@assays$RNA$data)

# Convert to data frames for merging
avg_exp_spe_df <- data.frame(gene = names(avg_exp_spe), avg_exp_spe = avg_exp_spe)
avg_exp_seurat_object1_df <- data.frame(gene = names(avg_exp_seurat_object1), avg_exp_seurat_object1 = avg_exp_seurat_object1)

# Merge the average expression data frames
merged_exp <- merge(avg_exp_spe_df, avg_exp_seurat_object1_df, by = "gene", all = TRUE)

# Identify top genes in both datasets
top_genes_spe <- head(merged_exp[order(-merged_exp$avg_exp_spe), "gene"], 10)
top_genes_seurat_object1 <- head(merged_exp[order(-merged_exp$avg_exp_seurat_object1), "gene"], 10)

# Combine the genes to label
genes_to_label <- union(top_genes_spe, top_genes_seurat_object1)

# Assuming 'avg_exp_spe' and 'avg_exp_seurat_object1' are already log-scaled
merged_exp$mean <- (merged_exp$avg_exp_spe + merged_exp$avg_exp_seurat_object1) / 2
merged_exp$difference <- merged_exp$avg_exp_spe - merged_exp$avg_exp_seurat_object1

# Compute the mean and standard deviation of the differences
mean_diff <- mean(merged_exp$difference, na.rm = TRUE)
sd_diff <- sd(merged_exp$difference, na.rm = TRUE)

# Identify significantly different genes based on the cutoff of 1.96 times the standard deviation from the mean difference
upper_cutoff <- mean_diff + 1.96 * sd_diff
lower_cutoff <- mean_diff - 1.96 * sd_diff
significantly_diff_genes <- merged_exp$gene[!is.na(merged_exp$gene) & (merged_exp$difference > upper_cutoff | merged_exp$difference < lower_cutoff)]

#significantly_diff_genes <- merged_exp$gene[!is.na(merged_exp$gene) & (merged_exp$difference < lower_cutoff)]

# Create labels for these significantly different genes
merged_exp$label <- ifelse(merged_exp$gene %in% significantly_diff_genes, as.character(merged_exp$gene), NA)

merged_exp <- merged_exp %>%
  mutate(is_Rp = grepl("^Rp", label))

# Print the list of significantly different genes
print("Significantly different genes:")
print(significantly_diff_genes)

# Create a Bland-Altman plot
bland_altman_plot <- ggplot(merged_exp, aes(x = mean, y = difference)) +
  geom_point(aes(color = is_Rp), alpha = 0.5) +  # Plot points with color based on the new column
  geom_hline(yintercept = mean_diff, color = "blue", linetype = "dashed") +  # Mean difference line
  geom_hline(yintercept = upper_cutoff, color = "red", linetype = "dashed") + 
  geom_hline(yintercept = lower_cutoff, color = "red", linetype = "dashed") +  # Upper and lower limit lines
  geom_text(aes(label = label), vjust = 1.5, size = 3, check_overlap = TRUE, na.rm = TRUE) +  # Add labels for significantly different genes only
  xlab("Mean Log Expression (Spatial and scRNA)") +
  ylab("Difference in Log Expression (Spatial - scRNA)") +
  ggtitle("Bland-Altman Plot for Discrepancies in Gene Expression") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove the legend
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "black"))  # Highlight "Rp" genes in red

# Print the plot
print(bland_altman_plot)


significantly_diff_genes <- significantly_diff_genes[!is.na(significantly_diff_genes)]


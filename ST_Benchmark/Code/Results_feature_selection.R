




blade_ARI <- data.frame(
  Method = c("200 Genes", "400 Genes", "600 Genes", "200 Remove Ribosomal Genes" , "200 Remove Platform Bias Genes", "200 More Upregulated Genes"),
  ARI = c(results_BLADE1$ARI, results_BLADE2$ARI, results_BLADE3$ARI, results_BLADE_norp$ARI, results_BLADE_sigdiff$ARI, results_BLADE_cutoff$ARI))

blade_sp_ARI <- data.frame(
  Method = c("200 Genes", "400 Genes", "600 Genes", "200 Remove Ribosomal Genes", "200 Remove Platform Bias Genes", "200 More Upregulated Genes"),
  ARI = c(bladesp200$ARI, bladesp400$ARI, bladesp600$ARI, bladespnorp$ARI, bladespsigdiff$ARI, bladespfc$ARI)
)


card_ARI <- data.frame(
  Method = c("200 Genes", "400 Genes", "600 Genes", "200 Remove Ribosomal Genes" , "200 Remove Platform Bias Genes", "200 More Upregulated Genes", "Full Gene Set"),
  ARI = c(results_card200$ARI, results_card400$ARI, results_card600$ARI, results_CARDnorp$ARI, results_CARD_sigdiff$ARI, results_CARD_cutoff$ARI, results_card_full$ARI))

music_ARI <- data.frame(
  Method = c("200 Genes", "400 Genes", "600 Genes", "200 Remove Ribosomal Genes" , "200 Remove Platform Bias Genes", "200 More Upregulated Genes", "Full Gene Set"),
  ARI = c(results_music200$ARI, results_music400$ARI, results_music600$ARI, results_MUSICnorp$ARI, results_MUSIC_sigdiff$ARI, results_MUSIC_cutoff$ARI, results_music_full$ARI))










combined_ARI <- rbind(
  data.frame(Method = blade_ARI$Method, ARI = blade_ARI$ARI, MethodType = "BLADE"),
  data.frame(Method = blade_sp_ARI$Method, ARI = blade_sp_ARI$ARI, MethodType = "spatialBLADE"),
  data.frame(Method = card_ARI$Method, ARI = card_ARI$ARI, MethodType = "CARD"),
  data.frame(Method = music_ARI$Method, ARI = music_ARI$ARI, MethodType = "MuSiC")
)

# Convert Method to factor with specified levels to maintain order
combined_ARI$Method <- factor(combined_ARI$Method, levels = unique(combined_ARI$Method))


ari <- ggplot(combined_ARI, aes(x = Method, y = ARI, color = MethodType, group = MethodType)) +
  geom_point(size = 3) +
  geom_line(aes(group = MethodType), size = 1) +
  theme_minimal() +
  labs(title = "ARI for Feature Selection Tests",
       x = "Gene Set",
       y = "Adjusted Rand Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.title = element_blank())








blade_correlation <- data.frame(
  Method = c("200 Genes", "400 Genes", "600 Genes", "200 Remove Ribosomal Genes", "200 Remove Platform Bias Genes", "200 More Upregulated Genes"),
  GC = c(results_BLADE1$GC, results_BLADE2$GC, results_BLADE3$GC, results_BLADE_norp$GC, results_BLADE_sigdiff$GC, results_BLADE_cutoff$GC),
  PGC = c(results_BLADE1$PGC, results_BLADE2$PGC, results_BLADE3$PGC, results_BLADE_norp$PGC, results_BLADE_sigdiff$PGC, results_BLADE_cutoff$PGC),
  MTC = c(results_BLADE1$MTC, results_BLADE2$MTC, results_BLADE3$MTC, results_BLADE_norp$MTC, results_BLADE_sigdiff$MTC, results_BLADE_cutoff$MTC),
  OSN = c(results_BLADE1$OSN, results_BLADE2$OSN, results_BLADE3$OSN, results_BLADE_norp$OSN, results_BLADE_sigdiff$OSN, results_BLADE_cutoff$OSN)
)

blade_sp_correlation <- data.frame(
  Method = c("200 Genes", "400 Genes", "600 Genes", "200 Remove Ribosomal Genes", "200 Remove Platform Bias Genes", "200 More Upregulated Genes"),
  GC = c(bladesp200$GC, bladesp400$GC, bladesp600$GC, bladespnorp$GC, bladespsigdiff$GC, bladespfc$GC),
  PGC = c(bladesp200$PGC, bladesp400$PGC, bladesp600$PGC, bladespnorp$PGC, bladespsigdiff$PGC, bladespfc$PGC),
  MTC = c(bladesp200$MTC, bladesp400$MTC, bladesp600$MTC, bladespnorp$MTC, bladespsigdiff$MTC, bladespfc$MTC),
  OSN = c(bladesp200$OSN, bladesp400$OSN, bladesp600$OSN, bladespnorp$OSN, bladespsigdiff$OSN, bladespfc$OSN)
)


card_correlation <- data.frame(
  Method = c("200 Genes", "400 Genes", "600 Genes", "200 Remove Ribosomal Genes", "200 Remove Platform Bias Genes", "200 More Upregulated Genes"),
  GC = c(results_card200$GC, results_card400$GC, results_card600$GC, results_CARDnorp$GC, results_CARD_sigdiff$GC, results_CARD_cutoff$GC),
  PGC = c(results_card200$PGC, results_card400$PGC, results_card600$PGC, results_CARDnorp$PGC, results_CARD_sigdiff$PGC, results_CARD_cutoff$PGC),
  MTC = c(results_card200$MTC, results_card400$MTC, results_card600$MTC, results_CARDnorp$MTC, results_CARD_sigdiff$MTC, results_CARD_cutoff$MTC),
  OSN = c(results_card200$OSN, results_card400$OSN, results_card600$OSN, results_CARDnorp$OSN, results_CARD_sigdiff$OSN, results_CARD_cutoff$OSN)
)

music_correlation <- data.frame(
  Method = c("200 Genes", "400 Genes", "600 Genes", "200 Remove Ribosomal Genes", "200 Remove Platform Bias Genes", "200 More Upregulated Genes"),
  GC = c(results_music200$GC, results_music400$GC, results_music600$GC, results_MUSICnorp$GC, results_MUSIC_sigdiff$GC, results_MUSIC_cutoff$GC),
  PGC = c(results_music200$PGC, results_music400$PGC, results_music600$PGC, results_MUSICnorp$PGC, results_MUSIC_sigdiff$PGC, results_MUSIC_cutoff$PGC),
  MTC = c(results_music200$MTC, results_music400$MTC, results_music600$MTC, results_MUSICnorp$MTC, results_MUSIC_sigdiff$MTC, results_MUSIC_cutoff$MTC),
  OSN = c(results_music200$OSN, results_music400$OSN, results_music600$OSN, results_MUSICnorp$OSN, results_MUSIC_sigdiff$OSN, results_MUSIC_cutoff$OSN)
)



# Combine the data frames into a single data frame
method_correlation <- rbind(
  data.frame(Method = blade_correlation$Method, GC = blade_correlation$GC, PGC = blade_correlation$PGC, MTC = blade_correlation$MTC, OSN = blade_correlation$OSN, MethodType = "BLADE"),
  data.frame(Method = blade_sp_correlation$Method, GC = blade_sp_correlation$GC, PGC = blade_sp_correlation$PGC, MTC = blade_sp_correlation$MTC, OSN = blade_sp_correlation$OSN, MethodType = "spatialBLADE"),
  
  data.frame(Method = card_correlation$Method, GC = card_correlation$GC, PGC = card_correlation$PGC, MTC = card_correlation$MTC, OSN = card_correlation$OSN, MethodType = "CARD"),
  data.frame(Method = music_correlation$Method, GC = music_correlation$GC, PGC = music_correlation$PGC, MTC = music_correlation$MTC, OSN = music_correlation$OSN, MethodType = "MuSiC")
)

# Convert Method to factor with specified levels to maintain order
method_correlation$Method <- factor(method_correlation$Method, levels = unique(method_correlation$Method))

# Display the data frame
print(method_correlation)



# Melt the combined data frame for easy plotting with ggplot2
combined_melted <- melt(method_correlation, id.vars = c("Method", "MethodType"), variable.name = "Metric", value.name = "Value")

# Filter to include only the 200 Genes method
filtered_data <- combined_melted[combined_melted$Method == "200 Remove Ribosomal Genes", ]
filtered_data2 <- combined_melted[combined_melted$Method == "200 Genes", ]

# Create a single boxplot for each method type with all metrics combined
cor <- ggplot(filtered_data, aes(x = MethodType, y = Value)) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = "white") +
  geom_jitter(position = position_jitter(0.2), size = 2, aes(color = Metric)) +
  scale_color_manual(values = c("GC" = "lightblue", "PGC" = "#ed78ff", "MTC" = "orange", "OSN" = "red")) +
  theme_minimal() +
  labs(title = "PCC for 200 DEG and Removal of Ribosomal Genes",
       x = "Method Type",
       y = "Value") +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        legend.position = "right")

cor2 <- ggplot(filtered_data2, aes(x = MethodType, y = Value)) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = "white") +
  geom_jitter(position = position_jitter(0.2), size = 2, aes(color = Metric)) +
  scale_color_manual(values = c("GC" = "lightblue", "PGC" = "#ed78ff", "MTC" = "orange", "OSN" = "red")) +
  theme_minimal() +
  labs(title = "PCC for 200 Genes",
       x = "Method Type",
       y = "Value") +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        legend.position = "right")

cor



final_plot <- (ari | cor) + plot_annotation(tag_levels = ("A")) & 
  theme(plot.tag = element_text(size = 15, face = "bold"))


# Create a function to generate the plot for each method and result
create_plot <- function(plot_object, title) {
  plot_object + theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),  # Center the title and make it smaller
      axis.title = element_blank(),  # Remove axis titles
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks = element_blank(),  # Remove axis ticks
      legend.position = "none"
    ) +
    labs(title = title) +
    geom_point(size = 2)  # Make points smaller
}

blade_plots <- list(
  create_plot(results_BLADE1$plot, "200 Genes"),
  create_plot(results_BLADE2$plot, "400 Genes"),
  create_plot(results_BLADE3$plot, "600 Genes"),
  create_plot(results_BLADE_norp$plot, "200 Remove Ribosomal Genes"),
  create_plot(results_BLADE_sigdiff$plot, "200 Remove Platform Bias Genes"),
  create_plot(results_BLADE_cutoff$plot, "200 More Upregulated Genes")
)

blade_sp_plots <- list(
  create_plot(bladesp200$plot, ""),
  create_plot(bladesp400$plot, ""),
  create_plot(bladesp600$plot, ""),
  create_plot(bladespnorp$plot, ""),
  create_plot(bladespsigdiff$plot, ""),
  create_plot(bladespfc$plot, "")
)

card_plots <- list(
  create_plot(results_card200$plot, ""),
  create_plot(results_card400$plot, ""),
  create_plot(results_card600$plot, ""),
  create_plot(results_CARDnorp$plot, ""),
  create_plot(results_CARD_sigdiff$plot, ""),
  create_plot(results_CARD_cutoff$plot, "")
)

music_plots <- list(
  create_plot(results_music200$plot, ""),
  create_plot(results_music400$plot, ""),
  create_plot(results_music600$plot, ""),
  create_plot(results_MUSICnorp$plot, ""),
  create_plot(results_MUSIC_sigdiff$plot, ""),
  create_plot(results_MUSIC_cutoff$plot, "")
)

# Combine plots for each method into rows
blade_combined <- wrap_plots(blade_plots, ncol = 6)
blade_sp_combined <- wrap_plots(blade_sp_plots, ncol = 6)
card_combined <- wrap_plots(card_plots, ncol = 6)
music_combined <- wrap_plots(music_plots, ncol = 6)

# Combine all method rows into a single plot
combined_plot <- (blade_combined / blade_sp_combined / card_combined / music_combined) +
  plot_layout(guides = "collect")

# Add a custom annotation to label the combined plot as "C"
final_plot <- combined_plot
# Print the final combined plot
print(final_plot)




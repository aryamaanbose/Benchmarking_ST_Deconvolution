

genes_to_plot <- c("Rps23", "Rpl23a", "Rpl21")  # Replace with your actual genes

# Define colors for idents in seurat_object1
idents_colors <- c("lightblue", "orange", "red", "#9999FF")

# Set the colors for Idents
Idents(seurat_object1) <- factor(Idents(seurat_object1), levels = levels(Idents(seurat_object1)))
seurat_object1@misc$colors <- idents_colors

# Set the identity in spe as "Spatial Data"

Idents(spe) <- "Spatial Data"
# Create a function to generate violin plots for a gene from both Seurat object



plot_violin_for_gene <- function(gene, seurat_object1, spe) {
  plot1 <- VlnPlot(seurat_object1, features = gene, pt.size = 0) +
    scale_fill_manual(values = idents_colors) +
    NoLegend() +  # Remove the legend
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10)) +
    ggtitle(paste(gene, "in scRNA Data"))
  
  plot2 <- VlnPlot(spe, features = gene, pt.size = 0) +
    scale_fill_manual(values = "palegreen3") +  # Set color to pastel green
    NoLegend() +  # Remove the legend
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10)) +
    ggtitle(paste(gene, "in Spatial Data"))
  
  combined_plot <- plot1 + plot2 + plot_layout(ncol = 2)  # Arrange plots side by side
  return(combined_plot)
}

# Generate and combine violin plots for all genes
combined_plots <- lapply(genes_to_plot, function(gene) {
  plot_violin_for_gene(gene, seurat_object1, spe)
})

# Combine all plots into one layout with 3 rows
final_plot <- wrap_plots(combined_plots, ncol = 1) + 
  plot_annotation(title = "Ribosomal Gene Expression in scRNA and Spatial Datasets") &
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 13)))

# Print the final combined plot
print(final_plot)


label_grob <- textGrob("C", x = unit(0.02, "npc"), y = unit(0.98, "npc"), 
                       gp = gpar(fontsize = 20, fontface = "bold"))

# Print the combined plot and add the label
grid.newpage()
grid.draw(final_plot)
grid.draw(label_grob)


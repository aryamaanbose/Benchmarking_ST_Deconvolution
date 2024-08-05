


####Cell type specific markers 



spe <- FindVariableFeatures(spe, selection.method = "vst", nfeatures = 2000)

spe <- ScaleData(spe, features = rownames(spe))

spe <- RunPCA(spe, features = VariableFeatures(object = spe))


spe <- FindNeighbors(spe, dims = 1:10)
spe<- FindClusters(spe,resolution = 0.5)



de_markers0 <- FindMarkers(spe, ident.1 = 0)
de_markers1 <- FindMarkers(spe, ident.1 = 1)
de_markers2 <- FindMarkers(spe, ident.1 = 2)
de_markers3 <- FindMarkers(spe, ident.1 = 3)

sorted_markers_GC <- de_markers1%>%
  filter(avg_log2FC >= 0) %>%
  arrange(p_val_adj)


sorted_markers_PGC <- de_markers0%>%
  filter(avg_log2FC >= 0) %>%
  arrange(p_val_adj)


sorted_markers_MTC <- de_markers3%>%
  filter(avg_log2FC >= 0) %>%
  arrange(p_val_adj)


sorted_markers_OSN <- de_markers2%>%
  filter(avg_log2FC >= 0) %>%
  arrange(p_val_adj)



GC_genes <- rownames(head(sorted_markers_GC, 10))

PGC_genes <- rownames(head(sorted_markers_PGC, 10))


MTC_genes <- rownames(head(sorted_markers_MTC, 10))


OSN_genes <- rownames(head(sorted_markers_OSN, 10))

SpatialFeaturePlot(spe, features = OSN_genes, pt.size = 10)


# Extract the gene expression values for each cell type marker gene
marker_genes <- list(
  GC = GC_genes,
  PGC = PGC_genes,
  OSNs =OSN_genes,
  M.TC = MTC_genes
)


calculate_correlations <- function(spe, marker_genes, proportions) {
  # Initialize results data frame
  correlation_results <- data.frame(CellType = character(), Gene = character(), Correlation = numeric())
  
  for (cell_type in names(marker_genes)) {
    if (cell_type %in% colnames(proportions)) {
      for (gene in marker_genes[[cell_type]]) {
        # Extract gene expression values
        if (gene %in% rownames(GetAssayData(object = spe, assay = "spatial", slot = "data"))) {
          gene_expression <- GetAssayData(object = spe, assay = "spatial", slot = "data")[gene, ]
        } else {
          next  # Skip this gene if it's not found
        }
        
        # Extract cell type proportions
        cell_type_proportions <- proportions[, cell_type]
        
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
      print(paste("Cell type", cell_type, "not found in proportions"))
    }
  }
  
  split_data <- split(correlation_results, correlation_results$CellType)
  
  total_correlation <- aggregate(Correlation ~ CellType, correlation_results, mean)
  mean_total_correlation <- mean(total_correlation$Correlation)
  
  # Create a list to store results
  results <- list()
  
  # Store total correlations per cell type
  for (i in 1:nrow(total_correlation)) {
    cell_type <- total_correlation$CellType[i]
    results[[cell_type]] <- total_correlation$Correlation[i]
  }
  
  # Store the mean of all total correlations
  results$total <- mean_total_correlation
  
  
  # Print total correlation values
  print("Total correlation values per cell type:")
  print(total_correlation)
  print(paste("Mean of all total correlations:", mean_total_correlation))
  
  return(results)
}

correlation <- calculate_correlations(spe, marker_genes, results_BLADE1$prop_matrix)


























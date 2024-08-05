load("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/PDAC/spatial_count.RData")

load("/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/PDAC/Figure4A_layer_annote.RData")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)

library(enrichplot)


# Create a Seurat object with spatial count data
sp <- CreateSeuratObject(counts = spatial_count)


common_cells <- intersect(rownames(layer_manual_PDAC), rownames(sp@meta.data))

# Subset Seurat object to include only cells with corresponding metadata
sp <- subset(sp, cells = common_cells)



sp <- AddMetaData(object = sp, metadata = layer_manual_PDAC, col.name = "Region")


genes_per_spot <- Matrix::colSums(GetAssayData(sp, assay = "RNA", slot = "counts") > 0)



VlnPlot(sp, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3 ,pt.size =0.2)


sp <- NormalizeData(sp, normalization.method = "LogNormalize", scale.factor = 10000)



sp <- ScaleData(sp)
sp<- FindVariableFeatures(sp)
sp<- RunPCA(sp)


sp <- FindNeighbors(sp, dims = 1:10)
sp<-RunUMAP(sp, dims = 1:10)


DimPlot(sp,reduction = "umap", group.by = "Region" )



Idents(sp) <- "Region"

# Find markers for each cluster
region_markers <- FindAllMarkers(sp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


perform_gsea <- function(gene_list, orgdb, keytype = "SYMBOL") {
  gsea_result <- enrichGO(gene          = gene_list,
                          OrgDb         = orgdb,
                          keyType       = keytype,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.001,
                          qvalueCutoff  = 0.05)
  return(gsea_result)
}


gsea_results <- list()
regions <- unique(region_markers$cluster)

for (region in regions) {
  genes <- region_markers %>% filter(cluster == region) %>% pull(gene)
  gsea_result <- perform_gsea(genes, org.Hs.eg.db)
  gsea_results[[as.character(region)]] <- gsea_result
}

# View the GSEA results
str(gsea_results)



gsea_results_long <- do.call(rbind, lapply(names(gsea_results), function(region) {
  result <- as.data.frame(gsea_results[[region]]@result)
  result$Region <- region
  result
}))

plot_gsea_results <- function(gsea_result, title) {
  plot <- dotplot(gsea_result, showCategory = 20) + ggtitle(title)
  print(plot)
}

for (region in names(gsea_results)) {
  plot_gsea_results(gsea_results[[region]], paste("GSEA for Region", region))
}













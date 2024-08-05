


here::i_am("Code/scrna_expressionprofile.R")


MOBSC_sce_raw <- MOBSC_sce_raw <- fread(here("Data/scrna_mob/Raw", "GSE121891_OB_6_runs.raw.dge.csv"))

meta_data <- read.csv(gzfile(here("Data/scrna_mob/Meta","GSE121891_OB_metaData_seurat.csv.gz")))

meta_data2 <- read.csv(here("Data/scrna_mob/Meta", "GSE121891_Figure_2_metadata.txt"), sep = "\t")

MOBSC_sce <- readRDS((url("https://figshare.com/ndownloader/files/40581983")))

raw_counts <- read_tsv(here("Data/spatial_mob/Raw","Spatial_MOB_Raw_counts.tsv"))



########################################################

raw_counts <- as.data.frame(raw_counts)

rownames(raw_counts) <- raw_counts[, 1]

raw_counts <- raw_counts[, -1]
spe <- CreateSeuratObject(counts = t(raw_counts), assay = "spatial")

#Extracting spatial coordinates from row names
coords <- strsplit(colnames(spe), "x")
coords <- matrix(unlist(coords), ncol = 2, byrow = TRUE)
coords <- apply(coords, 2, as.numeric) # Convert to numeric


transformation_matrix <- matrix(c(290.4, 0, 0, 0, 291.1, 0, -290.4, -291.1, 1), nrow = 3, byrow = TRUE)


# Adding a third column of 1's for affine transformation
homogeneous_coords <- cbind(coords, rep(1, nrow(coords)))

# Apply the transformation
transformed_coords <- homogeneous_coords %*% transformation_matrix

# Keeping only the transformed x and y coordinates
transformed_coords <- transformed_coords[, 1:2]


transformed_coords_df <- as.data.frame(transformed_coords)
names(transformed_coords_df) <- c("x", "y")
rownames(transformed_coords_df) <- rownames(raw_counts)


spe@images$image = new(Class = "SlideSeq", assay = "spatial", key = "he_image",coordinates = transformed_coords_df)


col_sums <- colSums(GetAssayData(spe))

# Identify spots with zero expression
spots_with_zero_expression <- which(col_sums == 0)


##Filterinf for number of unique genes in each spot check violin plot 
##Spots with low expression count 
spe <- subset(spe, nCount_spatial > 500)

length(colnames(spe))

spe <- NormalizeData(spe, normalization.method = "LogNormalize", scale.factor = 10000)




#################################################################################################################################

##scRNA Expressio profile generation for iterative removal of regulatory genes scenario 



rownames(MOBSC_sce_raw) <- MOBSC_sce_raw[[1]]  # Set the first column as row names


seurat_object1 <- CreateSeuratObject(counts = MOBSC_sce_raw)



seurat_object2 <- as.Seurat(MOBSC_sce) ## 182 genes 


seurat_object1 <- subset(seurat_object1, cells = colnames(seurat_object2))

# Extract cell names from both objects
cellnames_object1 <- rownames(seurat_object1@meta.data)
cellnames_object2 <- rownames(seurat_object2@meta.data)

# Check if the cell names match (you can use set operations to find mismatches)
mismatched_cells <- setdiff(cellnames_object1, cellnames_object2)

celltypes_data <- seurat_object2@meta.data$cell_type

# Add the 'celltypes' column to seurat_object1's metadata
seurat_object1@meta.data$celltypes <- celltypes_data[match(cellnames_object1, cellnames_object2)]

# Add FinalIds to meta_data
seurat_object1@meta.data$celltypes2 <- meta_data2$FinalIds[match(rownames(seurat_object1@meta.data), rownames(meta_data2))]

# Add ClusterName column to seurat_object1's metadata
##After subset only neuronal cell types are selected 
seurat_object1@meta.data$celltypes3 <- meta_data$ClusterName[match(rownames(seurat_object1@meta.data), meta_data$X)]

###GENES####
Mt_genes <-  grep("^Mt", rownames(seurat_object1), value = TRUE)

mt_genes <-  grep("^mt", rownames(seurat_object1), value = TRUE)

Rp_genes <- grep("^Rp", rownames(seurat_object1), value = TRUE)

rp_genes <- grep("^rp", rownames(seurat_object1), value = TRUE)

Malat1 <- grep("Malat1",  rownames(seurat_object1), value = TRUE)


# Add specific genes that are very highly expressed across all cells


genes_to_exclude <- c("Malat1",mt_genes)

# Invert the exclusion to get a list of genes to keep
genes_to_keep <- setdiff(rownames(seurat_object1), genes_to_exclude)

# Subset the Seurat object to keep only the desired genes
seurat_object1 <- subset(seurat_object1, features = genes_to_keep)

seurat_object1[["percent.mt"]] <- PercentageFeatureSet(seurat_object1, pattern = "^mt-")

seurat_object1<- subset(seurat_object1, subset = nFeature_RNA > 200 & percent.mt < 60)

##Filtering for genes with expression above zero, rowsums of counts data should be above zero 
genes_with_expression <- rownames(seurat_object1)[rowSums(GetAssayData(seurat_object1, slot = "counts")) > 0]

# Subset the Seurat object to only include these genes
seurat_object1 <- subset(seurat_object1, features = genes_with_expression)


###Normalize and scale####
seurat_object1 <- NormalizeData(seurat_object1, normalization.method = "LogNormalize", scale.factor = 10000)

###Find most variable genes 


seurat_object1 <- ScaleData(seurat_object1, features = rownames(seurat_object1))
# Run PCA
seurat_object1 <- FindVariableFeatures(seurat_object1, selection.method = "vst", nfeatures = 2000)
#
top10 <- head(VariableFeatures(seurat_object1, ),10)

plot1 <- VariableFeaturePlot(seurat_object1)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)




variable_features <- VariableFeatures(seurat_object1)
file_path <- "/Users/aryamaanbose/My_Drive/BLADE_for_SpatialDeconvolution/Processing/BLADE/Data/Variable_names.csv"
csv_string <- paste0('"', paste(variable_features, collapse = '","'), '"')


writeLines(csv_string, file_path)


seurat_object1_variable <- subset(seurat_object1, features = variable_features)
dim(seurat_object1_variable)


###DE gene selection###

Idents(seurat_object1) <- "celltypes"

##Take only positive fold change values 
markers_all <- FindAllMarkers(seurat_object1, only.pos = TRUE, min.pct = 0.1)

markers <- markers_all[markers_all$cluster == 'GC', ]
markersPGC <- markers_all[markers_all$cluster == 'PGC', ]
markersOSNS <- markers_all[markers_all$cluster == 'OSNs', ]
markersMTC <- markers_all[markers_all$cluster == 'M/TC', ]

create_volcano_plot <- function(data, log2fc_cutoff, pval_cutoff, dataset_name, n_labels = 10) {
  
  data$logP <- -log10(data$p_val_adj)
  data$Significant <- ifelse(data$p_val_adj < pval_cutoff & abs(data$avg_log2FC) > log2fc_cutoff, "Yes", "No")
  
  top_genes <- data %>%
    filter(Significant == "Yes") %>%
    arrange(desc(logP), desc(abs(avg_log2FC))) %>%
    slice_head(n = n_labels) %>%
    mutate(rank = row_number())  # Add ranking numbers
  
  p <- ggplot(data, aes(x = avg_log2FC, y = logP, color = Significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
    geom_text(data = top_genes, aes(label = paste(rank, gene, sep = ". ")), vjust = 2, size = 3, check_overlap = TRUE, color = "blue") + # Add ranking numbers to gene labels
    theme_minimal() +
    labs(title = paste("Volcano Plot for", dataset_name), x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), color = "blue", linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_cutoff), color = "blue", linetype = "dashed")
  return(p)
}

gc_volcano_plot <- create_volcano_plot(markers, log2fc_cutoff = 0, pval_cutoff = 0.01, dataset_name = "GC", n_labels = 55)
PGC_volcano_plot <- create_volcano_plot(markersPGC, log2fc_cutoff = 1, pval_cutoff = 0.01, dataset_name = "PGC", n_labels = 55)
OSN_volcano_plot <- create_volcano_plot(markersOSNS, log2fc_cutoff = 1, pval_cutoff = 0.01, dataset_name = "OSNs", n_labels = 55)
MTC_volcano_plot <- create_volcano_plot(markersMTC, log2fc_cutoff = 0.5, pval_cutoff = 0.01, dataset_name = "M/TC", n_labels = 55)


gc_volcano_plot + PGC_volcano_plot + OSN_volcano_plot + MTC_volcano_plot


VlnPlot(seurat_object1, features = c("Fth1","2010107E04Rik", "Pfdn2"))

VlnPlot(seurat_object1, features = c("Kif5b", "Zfos1", "Gm11808"))


#########Filtering##################




##
log2FC_cutoff <- 2

##By P value
sorted_markers_GC <- markers %>%
  filter(avg_log2FC >= 0) %>%
  arrange(p_val_adj) %>%
  head(100)


sorted_markers_PGC <- markersPGC %>%
  filter(avg_log2FC >= 0) %>%
  arrange(p_val) %>%
  head(100)

sorted_markers_OSNS <- markersOSNS %>%
  filter(avg_log2FC >= 0) %>%
  arrange(p_val) %>%
  head(100)

sorted_markers_MTC <- markersMTC %>%
  filter(avg_log2FC >= 0) %>%
  arrange(p_val) %>%
  head(100)


genes_PGC <- rownames(sorted_markers_PGC)
genes_GC <- rownames(sorted_markers_GC)
genes_OSNS <- rownames(sorted_markers_OSNS)
genes_MTC <- rownames(sorted_markers_MTC)



# Find common genes across all cell types
common_genes = Reduce(intersect, list(genes_GC, genes_PGC, genes_OSNS, genes_MTC))
# Combine all top markers into one list
gene_pool <- unique(c(genes_GC, genes_PGC, genes_OSNS, genes_MTC))

intersect(rownames(spe), gene_pool_autogene_variable)


ribo_genepool <- grep("^Rp", gene_pool, value = TRUE)

gene_pool_filtered <- gene_pool[!gene_pool %in% ribo_genepool]


gene_pool_autogene200 <- c('Slco5a1', 'Kcnb2', 'Gtf3c3', 'Pgap1', 'Resp18', 'Trip12', 'Tmem185b', 'Ipo9', 'Stx6', 'Ankrd45', 'Gorab', 'Pou2f1', 'Ildr2', 'Tstd1', 'Alyref2', 'Slc30a10', 
                           'Kcnh1', 'Acbd7', 'Camk1d', 'Stam', 'Man1b1', 'Abca2', 'Camsap1', 'Rexo4', 'Cacfd1', 'Setx', 'Gle1', 'Ccdc148', 'Slc25a12', 'Cdca7', 'Itgav', '1110051M20Rik',
                           'Lin7c', 'Lpcat4', 'Grem1', 'Trp53bp1', 'Anapc1', 'Gpcpd1', 'Cep250', 'Cnbd2', 'Zfp334', 'Cbln4', 'Gm6710', 'Zfp931', 'Oprl1', 'Wdr13', 'B630019K06Rik', 'Gria3',
                           'Mcf2', 'Mamld1', 'L1cam', 'Atp7a', 'Bex4', 'Rragb', 'Rab9', 'Ythdf3', 'Crh', 'Mfn1', 'Trpc3', 'Pgrmc2', 'Stoml3', 'Gucy1a3', 'Riiad1', 
                           'Man1a2', 'Kcna2', 'Hiat1', 'Pkn2', 'St6galnac3', 'Rmdn1', 'Gba2', 'Tmod1', 'Tmem246', 'AI314180', 'Astn2', 'Aldoart1', 'Mysm1', 'Magoh', '0610037L13Rik', 
                           'Txndc12', 'Map7d1', 'Trnau1ap', 'Nppa', 'Clcn6', 'Mxra8', 'Klhl17', 'Sema3e', 'Sema3c', 'Tmem60', 'Reln', 'Eif2b4', 'Gpn1', 'Slc4a1ap', 'Pisd', 'Sepsecs', 
                           'Paics', 'Wdfy3', 'Spp1', 'Tmed5', 'Fbrsl1', 'Ccdc60', 'Nos1', 
                           'Caln1', 'Slc12a9', 'Zfp12', 'Cpsf4', 'Hsph1', 'Thsd7a', 'D630045J12Rik', 'Chn2', 'Pde1c', 'Ndnf', 'Rnf103',
                           'Nagk', '1810044D09Rik', 'Trh', 'Foxp1', 'Rybp', 'Spsb2', 'Tspan9', 'Apold1', 'Lmo3', 'Cmas', 'St8sia1', 'Gm15706', 'Fam60a', 'Zscan18', '2310022A10Rik', 
                           'Sirt2', 'Zfp383', 'Arhgap33', 'Uri1', 'Atf5', 'Slc17a7', 'Kcnc1', 'Sv2b', 'Sh3gl3', 'Tmem126a', 'Omp', 'Prkrir', 'Dgat2', 'Tub', 'Copb1', 'Plk1', 'Msx3',
                           'Sirt3', 'Epm2a', 'Rspo3', 'Gopc', 
                           'Diras1', 'Pip5k1c', 'BC025920', 'Mettl25', 'Ppm1h', 'Ndufa4l2', 'Dlgap2', 'Alg11', 'Erlin2', 'Fcho1',
                           'Orc6', 'Phkb', 'Zfp612', 'Ccsap', 'Ttc13', 'Gm17296', 'Nrp1', 'Cacna1d', 'Cdhr1', 'Arhgef40', 'Mettl3', 'Haus4', 'Fgf9', 'Nefm', 'Dct', 'Cep164', 'AI593442', 'Igfbp3', 
                           'Pex13', 'Slu7', 'Slc16a11', 'Prpf8', 'B230217C12Rik', 'Jup', 'Mettl2', 'Pecam1', 'Sox9', 'Rnf157', 'Arl16', 'Rfng', 'Ryr2', 'Scgn', 'Msh3', 'Serf1', 'Asxl2', 'Nbas', 'Zfyve1', 'Rictor', 'Csmd3', 'Mtss1', 'Syt10', 'Slc4a8')


gene_pool_autogene400 <- c('Rrs1', 'Cpa6', 'Kcnb2', 'Fam135a', 'Zfp451', 'Tmem131', 'Rnf149', 'Pms1', 'Abcb6', 'Dnpep', 'Mrpl44', 'Psmd1', 'Arl4c', 'Gpc1', 'Dtymk', 'Mgat5', '2900009J06Rik', 'Rab3gap1', 'Ubxn4', 'Dars', 'Rbbp5', 'Rabif', 'Zbtb41', 'Scyl3', 'Nme7', 'Tiprl', 'Creg1', 'Pogk', 'Tomm40l', 'Tstd1', 'Pex19', 'Rgs7', 'Opn3', 'C130074G19Rik', '9630028B13Rik', 'Smyd2', 'Ppp2r5a', 'Rcor3', 'Acbd7', 'Cacnb2', 'Nebl', 'Otud1', 'Col5a1', 'Ralgds', 'Sh3glb2', 'Asb6', 'Ralgps1', 'Mapkap1', 'Acvr1', 'Csrnp3', 'Stk39', 'Dlx2', 'Dnajc10', 'Serping1', 'Lmo2', 'Eif3m', 'Lgr4', 'Ano3', 'Slc12a6', 'Grem1', 'Plcb2', 'Usp8', '1810024B03Rik', 'Crls1', 'Tbc1d20', 'Kif3b', 'Pxmp4', 'Rbl1', 'Ttpal', 'Dbndd2', 'Ncoa5', 'Cbln4', 'Gid8', 'Slc17a9', 'Arfgap1', 'Arfrp1', 'Akap17b', 'Gria3', 'Gabra3', 'G6pdx', 'Nhsl2', 'Magee1', 'Magt1', 'Hdx', 'Nap1l3', 'Armcx6', 'Kdm5c', 'Glra2', 'Nceh1', 'Cldn11', 'Kcnmb2', 'Zmat3', 'Slc7a11', 'Maml3', 'Stoml3', 'P2ry12', 'Golim4', 'Smg5', 'Rab13', 'Igsf3', 'Cept1', 'Gnat2', '5330417C22Rik', 'Rnpc3', 'Snx7', 'Pkn2', 'Sh3glb1', 'St6galnac5', 'Acadm', 'Penk', 'Decr1', 'Fbxl4', 'Mob3b', '1110017D15Rik', 'Sigmar1', 'Fam166b', 'Creb3', 'Gng10', 'Inip', 'Astn2', 'Adamtsl1', 'Caap1', 'Yipf1', '3110021N24Rik', 'Ppih', 'Hivep3', 'Rims3', 'Rlf', 'Hpcal4', 'Sf3a3', 'Ago1', 'Laptm5', 'Alpl', 'Ece1', 'Ubr4', 'Gm13157', 'Nppa', 'Srm', 'Slc45a1', 'Ski', 'Sema3c', '2700038G22Rik', 'Fam126a', 'Agap3', '1700096K18Rik', 'Hadha', 'Ppm1g', 'Bre', 'Yes1', 'Rgs12', 'Wfs1', 'Ccdc149', 'Chic2', 'Tmem165', 'Cep135', 'Spp1', 'Kctd10', 'Unc119b', 'Pxn', 'Nos1', 'Ptpn11', 'Naa25', 'Acad12', 'Ift81', 'Rnf34', 'Rilpl2', 'Fzd10', 'Phkg1', 'Wbscr17', 'Rasa4', 'Agfg2', 'Nudt1', 'Amz1', 'Radil', 'Pms2', 'Gng11', 'Thsd7a', 'Cald1', 'Zc3hav1', 'Hibadh', 'Wipf3', 'Vopp1', 'Krcc1', 'Mrpl35', 'Vamp8', 'Lrrtm4', 'Dguok', '2310040G24Rik', 'Aplf', 'H1fx', 'Zxdc', 'Trh', 'Lrig1', 'Cntn4', 'Hrh1', 'Kcna1', 'Prmt8', 'Tspan9', 'Rhno1', 'Cdkn1b', 'Apold1', 'Gm8994', 'Slco1a4', 'Cacng8', 'Brsk1', 'Isoc2b', 'Zfp580', 'Zfp772', 'Ehd2', 'Nanos2', 'Mypop', 'Ppp1r37', 'Map3k10', 'Zfp60', 'Fbxo17', 'Capn12', 'Sipa1l3', 'Sdhaf1', 'Lrfn3', 'Gramd1a', 'Atf5', 'Slc17a7', 'Pth2', 'Grwd1', 'Lysmd4', 'Mrpl46', 'Polg', 'Nmb', 'Gm2115', 'Crebzf', 'Ccdc90b', 'Omp', 'Prkrir', 'Rnf169', 'Inppl1', 'Nup98', 'Trim3', 'Stk33', 'Arntl', 'Spns1', 'Asphd1', 'Kif22', 'Sephs2', 'Bckdk', 'Mettl10', 'Zfp941', 'Katna1', 'Grm1', 'Nt5dc1', 'Cdc40', 'Cd164', 'Ostm1', 'Gopc', 'Egr2', 'Pcnt', 'Apc2', 'Dapk3', 'Zfr2', 'Sirt6', 'Igf1', 'Tmcc3', 'Kitl', 'Slc6a15', 'Osbpl8', 'Kcnc2', 'Cct2', 'Usp15', 'Mettl1', 'March9', 'Pip4k2c', 'Ctxn1', 'Whsc1l1', 'Spcs3', 'Naf1', 'Ints10', 'Gmip', '1700030K09Rik', 'Rasd2', 'Nr3c2', 'Gipc1', 'Dcaf15', 'Nacc1', 'Cbln1', 'Gm2694', 'Fam192a', 'Pllp', 'Cdh11', 'Acd', 'Nrn1l', 'Slc7a6os', 'Cyb5b', 'Atxn1l', 'Kcng4', 'Zc3h18', 'Nup133', 'Flnb', 'Psmd6', 'Fam149b', '1700112E06Rik', 'Rft1', 'Cdhr1', 'Exoc5', '3632451O06Rik', 'Apex1', 'Emc9', 'Shisa2', 'Sacs', 'Ptk2b', 'Gm6878', 'Fndc3a', 'Ccdc82', '1700012B09Rik', 'Panx1', 'Tmed1', 'Dpy19l1', 'Snx19', 'Tirap', 'Spa17', 'Sc5d', 'Tmem136', 'Hmbs', 'Scn2b', 'Sidt2', '2900052N01Rik', 'Pts', 'Dlat', 'AI593442', 'Cib2', 'Etfa', 'Mpi', 'Thsd4', 'Ibtk', 'Dopey1', 'Rbp1', 'Slc35g2', 'Aste1', 'Parp3', 'Qrich1', 'Eomes', 'Csrnp1', 'Trak1', 'Tmem42', 'Ap1b1', 'Grb10', 'Lgalsl', 'Acyp2', 'Pwwp2a', 'Med7', 'Trim7', 'Nhp2', 'Cdkn2aipnl', 'Sept8', 'Med9', 'Wscd1', 'Hic1', 'Blmh', 'Spag5', 'Pdk2', 'Stard3', 'Nbr1', 'Tmub2', 'Pitpnc1', 'Unk', 'P4hb', 'Notum', 'Fn3k', 'Ero1lb', 'Tbce', 'Cdk13', 'Hist1h1c', 'Cap2', 'Spock1', 'Ptdss1', 'Adcy2', 'Gm9776', 'Fam169a', 'Taf9', 'Mrps30', 'Nol10', 'Hectd1', 'Daam1', 'Max', 'Ift43', 'Itpk1', 'Plcxd3', 'Sdc2', 'Tspyl5', 'Ncald', 'Emc2', 'Ebag9', 'Mtss1', 'Slc52a2', 'Arhgap39', 'Nol12', 'Cacna1i', 'Zc3h7b', 'Tcf20', 'Scube1', 'Syt10', 'Arid2')


gene_pool_autogene400_1 <- c('Oprk1', 'Cspp1', 'Cpa6', 'Ly96', 'Mrpl30', 'Wdr75', 'Mfsd6', 'Hecw2', 'Pgap1', 'Trak2', 'Wdr12', 'Raph1', 'Pard3b', 'Klf7', 'Mettl21a', 'March4', 'Ndufa10', 'Atg4b', 'Gin1', '2310035C23Rik', 'Tsn', 'Klhdc8a', 'Mdm4', 'Tnni1', 'Gm19705', 'Swt1', 'Rnf2', 'Uhmk1', 'Dusp12', 'Usf1', 'Tstd1', 'Kcnj10', 'Opn3', 'Hnrnpu', 'Smyd3', 'Lpgat1', 'Acbd7', 'Hspa14', 'Cacnb2', 'Gad2', 'Nsmf', 'Rexo4', 'Olfm1', 'Dolk', 'Fam73b', 'Scai', 'Tnfaip6', 'Kcnh7', 'Nostrin', 'Tlk1', 'Sp9', 'Cir1', 'Agps', 'Ypel4', 'Arfgap2', 'Cry2', 'Ext2', 'Fjx1', 'Cstf3', 'Ccdc73', 'Mpped2', 'Ano3', 'Grem1', 'Rasgrp1', 'Ankrd63', 'Disp2', 'Ehd4', 'Secisbp2l', 'Mal', 'Ndufaf5', 'Polr3f', 'Entpd6', 'Nsfl1c', 'Pdrg1', 'BC029722', 'Stk4', 'Dpm1', 'Bcas1', 'Cbln4', 'Fam210b', 'Stx16', 'Osbpl2', 'Suv39h1', 'Ebp', 'Cybb', 'Uxt', 'Sept6', 'Upf3b', 'Gabra3', 'L1cam', 'Vbp1', 'Il1rapl1', 'Pola1', 'Arhgef9', 'Efnb1', '5330434G04Rik', 'Tceal6', 'Tceal8', 'Morc4', 'Nxt2', 'Prps2', 'Ralyl', 'Crh', 'Nlgn1', 'Il2', 'Fat4', 'Lhfp', 'Stoml3', 'Nbea', 'Rnf13', 'Tiparp', 'Fstl5', 'Ctso', 'Gucy1a3', 'Fcrls', 'S100a4', 'Riiad1', 'Scnm1', 'Syt6', 'Kcna2', 'Gpr61', 'Ntng1', 'Vcam1', 'Agl', 'Rwdd3', 'Abcd3', 'Cxxc4', 'Metap1', 'Gtf2b', 'Odf2l', 'Lpar3', 'Spata1', 'Penk', 'Triqk', 'Otud6b', 'Rngtt', '1110017D15Rik', 'Enho', 'Igfbpl1', 'Stx17', 'Rad23b', 'Rgs3', 'Kdm4c', 'Ptprd', 'Lurap1l', 'Nfib', 'Slc24a2', 'Nfia', 'Alg6', 'Dab1', 'Nasp', 'Mycbp', 'Meaf6', 'Ubxn11', 'Syf2', 'Nipal3', 'Kif17', 'Emc1', 'Nppa', 'Srm', 'Tmem201', 'Fam213b', 'Isg15', 'Sema3a', 'Pclo', 'Tmem60', 'Reln', 'Lhfpl3', 'Eif2b4', 'Ppm1g', 'Zfp512', 'Jakmip1', 'Drd5', 'Ube2k', 'Sgcb', 'Epha5', 'Rufy3', 'Prkg2', 'Spp1', 'Zfp951', 'Pxmp2', 'Gltp', 'Cit', 'Aldh2', 'Glt1d1', 'Vkorc1l1', 'Clip2', 'Rfc2', 'Pdgfa', '0610040B10Rik', 'Tmem130', 'Nxph1', 'Thsd7a', 'Ube2h', 'Parp12', 'Cntnap2', 'Lancl2', 'Gng12', '1600020E01Rik', 'Mgll', 'Zxdc', 'Nup210', 'Trh', 'Magi1', 'Crbn', 'Arl8b', 'Lhfpl4', 'Tmcc1', 'Tapbpl', 'Apold1', 'Grin2b', 'Dennd5b', 'Bicd1', 'Zfp628', 'Fiz1', 'Sae1', 'Slc1a5', 'Nova2', 'Eml2', 'Zfp111', 'Cic', 'Tmem91', 'Ltbp4', 'Pld3', 'Map4k1', 'Fxyd1', 'Atf5', 'Med25', 'Bcl2l12', 'Slc17a7', 'Pth2', 'Kdelr1', 'Igf1r', 'Nmb', 'Me3', 'Picalm', 'Tmem126a', 'Omp', 'Rnf121', 'Art5', 'Trim34a', 'Mrvi1', 'Spon1', 'Ubfd1', 'Prkcb', 'Zfp771', 'Pycard', 'Rgs10', 'Inpp5f', 'Nsmce4a', 'Plekha1', 'Bccip', 'Fuom', 'Rspo3', 'Frk', 'Fyn', 'H2afy2', 'Kiss1r', 'Apba3', 'Pip5k1c', 'Timp3', 'Arl1', 'Kitl', 'Slc6a15', 'Mak16', 'Nrg1', 'Trappc11', 'Hmgb2', 'Npy1r', 'Ssbp4', 'Zfp827', 'Inpp4b', 'Ccdc130', 'Chd9', 'Dhodh', 'Taf1c', 'Zc3h18', 'Rnf166', 'Spg7', 'Gnpat', 'Cdhr1', '3632451O06Rik', 'Prmt5', 'Ap1g2', 'Cpne6', 'Xpo4', 'Kcnrg', 'Pnma2', 'Wbp4', 'Uchl3', 'Gpr180', 'Mbnl2', '1700012B09Rik', 'Zfp558', 'Icam5', 'Scn3b', 'Pvrl1', 'Cadm1', 'AI593442', 'Lingo1', 'Cln6', 'Pias1', 'Tpm1', 'Gtf2a2', 'Tcf12', 'Nedd4', 'Ccpg1', 'Rsl24d1', 'Fbxo9', 'Cyb5r4', 'Bcl2a1b', 'Zic1', 'Chst2', 'Atr', 'Copb2', 'Mrps22', 'Pik3cb', 'Dock3', 'Hemk1', 'Hyal2', 'Slc38a3', 'Ip6k1', 'Slc25a20', 'Cck', 'Nefh', 'Aftph', 'Peli1', 'Psme4', 'Gabra1', 'Itk', 'Tbc1d9b', '0610009B22Rik', 'Igtp', 'Efnb3', 'Camkk1', 'Inpp5k', 'Vps53', 'Coro6', 'Tmem199', 'Psmd11', 'Myo1d', 'Aatf', 'Mrps23', 'Xylt2', 'Ngfr', 'Zfp652', 'Prkca', 'Rgs9', 'Btbd17', 'Cygb', 'Cyth1', 'Rbfox3', 'Baiap2', 'Mrpl12', 'Fam195b', 'Alyref', 'Metrnl', 'Hecw1', 'Zfp322a', 'Tmem170b', 'Bicd2', 'Spin1', 'Pdlim7', 'Habp4', 'Papd7', 'Lysmd3', 'Rasa1', 'Atp6ap1l', 'Zcchc9', 'Dhfr', 'Foxd1', 'Bdp1', 'Mocs2', 'Hmgcs1', 'Stxbp6', 'Cfl2', 'Srp54b', 'Klhdc1', 'Atp5s', 'Jkamp', 'Sel1l', 'Eml5', 'Ppp4r4', 'Papola', 'Ppp2r5c', 'Amn', 'Ctnnd2', 'Baalc', 'Col22a1', 'Lynx1', 'Ly6e', 'Exosc4', 'Cyhr1', 'Tst', 'Pick1', 'Tmem184b', 'Mgat3', 'Scube1', 'Kif21a', 'Kcnh3')


gene_pool_autogene_variable <- c('Gls', 'Unc80', 'Igfbp2', 'Igfbp5', 'Resp18', 'Serpine2', 'Pid1', 'Dner', 'Pam', 'Tmem163', 'R3hdm1', 'Dnm3', 'Pcp4l1', 'Pea15a', 'Atp1a2', 'Rgs7', 'Acbd7', 'Gad2', 'Gpsm1', 'Olfm1', 'Dnm1', 'Stxbp1', 'Strbp', 'Cacnb4', 'Mettl5os', 'Slc1a2', 'Grem1', 'Disp2', 'Chgb', 'Napb', 'Map1lc3a', 'Ptprt', 'Ywhab', 'Slc12a5', 'Tshz2', 'Eef1a2', 'Dnajc5', 'Syp', 'Stmn2', 'Crh', 'Pcdh10', 'Lhfp', 'Dclk1', 'Pfn2', 'Serpini1', 'Fstl5', 'Glrb', 'Gucy1a3', 'Syt11', 'Npr1', 'S100a5', 'Syt6', 'Kcna2', 'Ntng1', 'Prss12', 'Ppp3ca', 'Rap1gds1', 'Ttll7', 'Penk', 'Fut9', 'Lingo2', '1110017D15Rik', 'Enho', 'Ptprd', 'Nfib', 'Dnajc6', 'Rab3b', 'Marcksl1', 'C1qb', 'Camk2n1', 'Nppa', 'Clstn1', 'Pclo', 'Sema3c', 'Reln', 'Srpk2', 'Rheb', 'Ywhah', 'Nat8l', 'Mrfap1', 'Rasl11b', 'Spp1', 'Ywhag', 'Dync1i1', 'Nxph1', 'Thsd7a', 'Cntnap2', 'Cbx3', 'Nap1l5', 'Ndnf', 'Trh', 'Slc6a1', 'Syn2', 'Clstn3', 'Tpi1', 'Pianp', 'Gabarapl1', 'Cdkn1b', 'Apold1', 'Peg3', 'Scn1b', 'Atf5', 'Slc17a7', 'Kcnc1', 'Mef2a', 'Akap13', 'Nmb', 'Pak1', 'Fam168a', 'Hbb-bs', 'Spon1', 'Rgs10', 'Cend1', 'Ctsd', 'Th', 'Vip', 'Rgs17', 'Grm1', 'Cited2', 'Rspo3', 'Tspyl4', 'Spock2', 'Arid5b', 'Zwint', 'Nt5dc3', 'Btg1', 'Mcf2l', 'Cpe', 'Rab3a', 'Gm2694', 'Mt2', '6430548M08Rik', 'Tubb3', 'Cdhr1', 'Ghitm', 'Cpne6', 'Ctsb', 'Clu', 'Pcdh17', 'Thy1', 'Fxyd6', 'AI593442', 'Map2k1', 'Tmem30a', 'Elovl4', 'Tpbg', 'Zic1', 'Clstn2', 'Cpne4', 'Arpp21', 'Cck', 'Gabrg2', 'Gabra1', 'Pafah1b1', 'Cltc', 'Dynll2', 'Msi2', 'Atp6v0a1', 'Mapt', 'Limd2', 'Prkar1a', 'Grb2', 'Cyth1', 'Timp2', 'Gng4', 'Nrsn1', 'Sox4', 'Nrn1', 'Cltb', 'Sncb', 'Cxcl14', 'Map1b', 'Rab3c', 'Vsnl1', 'Id2', 'Tspan13', 'Etv1', 'Stxbp6', 'Syndig1l', 'Chga', 'Sp8', 'Sepp1', 'Slc1a3', 'Basp1', 'Ywhaz', 'Syngr1', 'Sept3', 'Slc38a1', 'Faim2', 'Apod', 'Gsk3b', 'Zbtb20', 'Atp6v1a', 'Tagln3', 'Abi3bp', 'Robo2', 'Pcp4', 'Gng13', 'Ly6g6e', 'Atp6v1g2', 'Tubb5', 'Pja2', 'Dlgap1', 'Kif5b', 'Reep5', 'Slc6a7', 'Atp5a1', 'Doc2g', 'Prdx5', 'Rtn3', 'Gng3', 'Syt7', 'Ms4a15', 'Prune2', 'Pgam1', 'Got1', 'Kcnip2', 'Pdcd4')
length(gene_pool_autogene_variable)


grep("^Rp", variable_features, value = TRUE)

intersect(rownames(seurat_object1), gene_pool_autogene_variable)


intersect(rownames(spe),gene_pool_autogene400)

setdiff(gene_pool_autogene400, rownames(spe))



seurat_object1_filtered <- subset(seurat_object1, features = deg)
rownames(seurat_object1_filtered)


dim(seurat_object1_filtered)

spe_filtered <- subset(spe, features = rownames(seurat_object1_filtered))
seurat_object1_filtered <- subset(seurat_object1_filtered, features = rownames(spe_filtered))

#####BLAND ALTMAN####







avg_exp_spe <- rowMeans(spe_filtered@assays$spatial$data)
avg_exp_seurat_object1 <- rowMeans(seurat_object1_filtered@assays$RNA$data)

grep("^Rp", rownames(seurat_object1_filtered), value = TRUE)

dim(seurat_object1_filtered_noribo)

merged_exp <- merge(avg_exp_spe, avg_exp_seurat_object1, by = "row.names", all = TRUE)
names(merged_exp) <- c("gene", "avg_exp_spe", "avg_exp_seurat_object1")


merged_exp$mean <- (merged_exp$avg_exp_spe + merged_exp$avg_exp_seurat_object1) / 2
merged_exp$difference <- merged_exp$avg_exp_spe - merged_exp$avg_exp_seurat_object1

# Compute the mean and standard deviation of the differences
mean_diff <- mean(merged_exp$difference, na.rm = TRUE)
sd_diff <- sd(merged_exp$difference, na.rm = TRUE)

# Identify significantly different genes based on the cutoff of 1.96 times the standard deviation FROM THE MEAN DIFFERENCE
upper_cutoff <- mean_diff + 1.96 * sd_diff
lower_cutoff <- mean_diff - 1.96 * sd_diff
significantly_diff_genes <- merged_exp$gene[merged_exp$difference > upper_cutoff | merged_exp$difference < lower_cutoff]

# Create labels for these significantly different genes
merged_exp$label <- ifelse(merged_exp$gene %in% significantly_diff_genes, as.character(merged_exp$gene), NA)

# Identify ribosomal genes
merged_exp$ribosomal <- ifelse(merged_exp$gene %in% ribo_genepool, as.character(merged_exp$gene), NA)

# Create a Bland-Altman plot
bland_altman_plot <- ggplot(merged_exp, aes(x = mean, y = difference)) +
  geom_point(alpha = 0.5) +  # Plot all points with some transparency
  geom_hline(yintercept = mean_diff, color = "blue", linetype = "dashed") +  # Mean difference line
  geom_hline(yintercept = upper_cutoff, color = "red", linetype = "dashed") + 
  geom_hline(yintercept = lower_cutoff, color = "red", linetype = "dashed") + # Upper and lower limit lines
  geom_text(aes(label = label), vjust = 1.5, size = 3, check_overlap = TRUE, na.rm = TRUE) + # Add labels for significantly different genes
  geom_text(aes(label = ribosomal), vjust = -1.5, size = 3, check_overlap = TRUE, na.rm = TRUE) +# Add labels for ribosomal genes
  xlab("Mean Log Expression (Spatial and scRNA)") +
  ylab("Difference in Log Expression (Spatial - scRNA)") +
  ggtitle("Bland-Altman Plot for Discrepancies in Gene Expression") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(bland_altman_plot)


VlnPlot(seurat_object1, features = c("Rpl30","Rps28", "Rpl23a"))

VlnPlot(spe, features = c("Rpl30","Rps28", "Rpl23a"))


VlnPlot(seurat_object1, features = significantly_diff_genes)

VlnPlot(spe, significantly_diff_genes)

#seurat_object1_filtered <- subset(seurat_object1_filtered, features = setdiff(rownames(seurat_object1_filtered), significantly_diff_genes))












##Only works on seurat_object1_filtered
seurat_object1_filtered <- subset(seurat_object1, features = deg)

dim(seurat_object1_filtered)

length(deg)


####SAVING#####

####Making expression profiles####
# Function to calculate standard deviation for each gene
calc_std_dev <- function(seurat_obj) {
  data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  return(rowSds(as.matrix(data)))
}

# Function to calculate mean expression for each gene
calc_mean_expr <- function(seurat_obj) {
  data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  return(rowMeans(as.matrix(data)))
}

# Function to calculate variance for each gene
calc_variance <- function(seurat_obj) {
  data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  return(rowVars(as.matrix(data)))
}

calc_corrected_std <- function(mean_matrix, variance_matrix) {
  gene_list <- rownames(mean_matrix)
  ct_list <- colnames(mean_matrix)
  New_std <- matrix(nrow = length(gene_list), ncol = length(ct_list))
  
  for (i in seq_along(ct_list)) {
    trend <- fitTrendVar(mean_matrix[, i], variance_matrix[, i])$trend
    New_std[, i] <- sqrt(trend(mean_matrix[, i]))
  }
  
  rownames(New_std) <- gene_list
  colnames(New_std) <- ct_list
  return(New_std)
}

# Split the Seurat object by cell type
seurat_list <- SplitObject(seurat_object1_filtered, split.by = "celltypes")

# Apply the function to each cell type subset
std_dev_list <- lapply(seurat_list, calc_std_dev)

# Combine into a single data matrix
std_dev_matrix <- do.call(cbind, std_dev_list)

# Name the columns as cell types
colnames(std_dev_matrix) <- names(std_dev_list)

# Extract gene names
gene_names <- rownames(GetAssayData(seurat_object1_filtered, assay = "RNA", slot = "data"))

# Assign gene names to the rows
rownames(std_dev_matrix) <- gene_names

# Apply the function to each cell type subset
mean_expr_list <- lapply(seurat_list, calc_mean_expr)

# Combine into a single data matrix
mean_expr_matrix <- do.call(cbind, mean_expr_list)

# Name the columns as cell types
colnames(mean_expr_matrix) <- names(mean_expr_list)

# Assign gene names to the rows (using the same gene_names variable from before)
rownames(mean_expr_matrix) <- gene_names


# Apply the function to each cell type subset
variance_list <- lapply(seurat_list, calc_variance)

# Combine into a single data matrix
variance_matrix <- do.call(cbind, variance_list)

# Name the columns as cell types
colnames(variance_matrix) <- names(variance_list)

# Assign gene names to the rows (using the same gene_names variable from before)
rownames(variance_matrix) <- gene_names

New_std <- calc_corrected_std(mean_expr_matrix, variance_matrix)







####SAVING#####

data_slot <- GetAssayData(spe, slot = "data")
data_slot <- GetAssayData(spe, slot = "data")
data_frame <- as.data.frame(as.matrix(data_slot))


#scrna_fullvariable <- GetAssayData(seurat_object1_variable, slot = "counts")

datasets <- list(
  mean_expr_matrix = "deg_mean600.csv",
  New_std = "deg_sd600.csv"
)

# Base directory within your project structure
base_dir <- here("Data", "Processed_BLADE")


# Iterate over the datasets and their file names
for (data_name in names(datasets)) {
  # Construct the file path
  file_path <- file.path(base_dir, datasets[[data_name]])
  
  # Ensure the directory exists
  if (!dir.exists(dirname(file_path))) {
    dir.create(dirname(file_path), recursive = TRUE)
  }
  
  # Write the dataset to a CSV file
  write.csv(get(data_name), file = file_path, row.names = TRUE)
}







dim(mean_expr_matrix)




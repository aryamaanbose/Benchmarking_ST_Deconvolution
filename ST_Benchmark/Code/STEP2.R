#####STEP2: Simulate scRNA DATA based on a reference data

seurat_object1_markers <- subset(seurat_object1, features = deg_norp)

dim(seurat_object1_markers)


sc <-  as.SingleCellExperiment(seurat_object1_markers)

set.seed(42)

MOBSC_data <- construct_data(
  sce = sc,
  assay_use = "counts",
  celltype = "celltypes",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "1"
)

MOBSC_marginal <- fit_marginal(
  data = MOBSC_data,
  predictor = "gene",
  mu_formula = "celltypes",
  sigma_formula = "celltype",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  parallelization = "pbmcmapply"
  
)

MOBSC_copula <- fit_copula(
  sce = sc,
  assay_use = "counts",
  marginal_list = MOBSC_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 2,
  input_data = MOBSC_data$dat
)

MOBSC_para <- extract_para(
  sce = sc,
  marginal_list = MOBSC_marginal,
  n_cores = 2,
  family_use = "nb",
  new_covariate = MOBSC_data$newCovariate,
  data = MOBSC_data$dat
)

MOBSC_newcount <- simu_new(
  sce = sc,
  mean_mat = MOBSC_para$mean_mat,
  sigma_mat = MOBSC_para$sigma_mat,
  zero_mat = MOBSC_para$zero_mat,
  quantile_mat = NULL,
  copula_list = MOBSC_copula$copula_list,
  n_cores = 2,
  family_use = "nb",
  input_data = MOBSC_data$dat,
  new_covariate = MOBSC_data$newCovariate,
  filtered_gene = MOBSC_data$filtered_gene
)



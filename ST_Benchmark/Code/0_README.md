### CODE

<u>*Explanation*</u>:

*This directory contains the code executed on R, notebooks are available as .Rmd files. 

The `here` package was used to automate file loading and saving.

This repository is linked to the Data subrepository, which contains the following directories:

- `scrna_mob`: Contains scRNA data.
- `spatial_mob`: Contains spatial data.
- `PDAC`: Contains files for PDAC analysis.
- `Processed_BLADE`: Contains processed files for the execution of BLADE.


**Documentation files**: 

- `Readme.MD` in 'ST_Benchmark' Contains instructions for executing the scripts in the ST_Benchmark project.
The code was tested using a combination of real and simulated data, with accuracy metrics calculated to evaluate performance.

**Script 1:** `Preprocess_ST.R`  
*Description:* This script handles the loading and preprocessing of MOB spatial data, including metadata of ground truth annotations for spots.

**Script 2:** `Preprocess_scRNA.R`  
*Description:* This script manages the loading and preprocessing of MOB scRNA data, including metadata in the form of cell type annotations.

**Script 3:** `Real_data_analysis`  
*Description:* Implements the evaluation of real ST data with accuracy metrics such as Adjusted Rand Index (ARI) and Pearson Correlation Coefficient (PCC), executed within the `process_data()` function.

**Script 4:** `Simulation Framework Scripts`  
*Description:* For simulated data, these scripts must be run in sequence as `STEP1`, `STEP2`, `STEP3`, and `STEP4`. The deconvolution tools are applied to this data in `STEP4`.

**Script 5:** `Simulation_Results.R`  
*Description:* Includes the `run_analysis()` function for calculating accuracy metrics for simulated data, such as Adjusted Rand Index, Node Purity, Morans I Difference, and spot-wise and cell-type-wise RMSE.

**Script 6:** `PDAC analysis scripts`  
*Description:* The `Deconvolution_PDAC` script executes the analysis and visualization for PDAC data. The deconvolution results are validated using `GSEA.R` on regions of interest.

Execute the analysis as per the steps outlined in the `Readme_ST_Benchmark`.

Packages and Versioning of the analysis is saved in SessionInfo.txt

### Benchmarking Spatial Transcriptomics Tools Deconvolution Tools that use a Reference scRNA Profile, Using Real and Simulated Data

**Main contact:**
Aryamaan Bose, a.bose@student.vu.nl

#### Directories:

- **ST_Benchmark**:  
  Contains all items related to the execution of frameworks based on R, including the preprocessing scripts and analysis files.

- **ST_Benchmark/spatial_BLADE**:  
  Holds the Python analysis scripts and framework files for BLADE. The `runBLADE.py` script is used to implement specific functionalities based on the functions defined.

#### Dependencies:
The Python scripts depend on the outputs from the R scripts. Specifically, the R scripts must be executed first to generate the necessary files for the Python-based BLADE deconvolution analysis. 

The execution for deconvolution tools based on R mainly (`MuSiC`, `CARD`, `RCTD`, and `CibersortX`) are present in the 'Code' directory.

#### How to Execute the Code:

1. **Running R Scripts**:
   - Navigate to the `ST_Benchmark` directory.
   - Execute `SETUP.R` to install and load necessary packages.
   - Execute `preprocess_ST.Rmd` followed by `preprocess_scRNA.Rmd` to generate the required files for deconvolution analysis. `preprocess_scRNA.Rmd` contains the function `process`.
   - To remove genes with platform bias, first run `Platform_bias.R`.
   - To create expectation matrices for BLADE, run `Expectation_Matrix_BLADE.R` and specify the dominant cell type proportions for each cluster in the script.

2. **Real Data Analysis**
   - For applications of `CARD` and `MuSiC` on real data, run `CARD_Music_realdata.R`.
   - For real data analysis, run the `Real_data_analysis.Rmd` notebook. Files for BLADE (from Python) should be processed with `process_dataset`, and R-based tools should be processed with `process_dataset2`. 
   - To generate graphs for expectation matrix tests, run the file `Results_SP.R`.
   - To generate graphs for the feature selection test, run the file `Results_feature_selection.R`.

3. **Simulated Data Analysis**
   - Simulation is performed using `scDesign3` to estimate cell type-specific distributions from reference scRNA data. In our case, we have defined the feature set in `preprocess_scRNA` and then scaled the sampled counts from these distributions with a predefined proportions matrix (emulating spatial autocorrelation) to simulate synthetic ST data. 
   - `STEP1.R`: For defining proportions. Calculate `calculate_optimal_radius()` using an input as the number of neighbors, followed by `compute_proportions6()` where the user can define parameters `specified_dominant_proportions`, `neighbourhood_radius` (radius where nearby spots influence each other), `sc` (random noise, subtract `+sc`, `-sc` from the defined dominant cell type proportions), and `pn` (`pn%` of spots will be shuffled for their dominant cell type annotation).
   - `STEP2.R`: Estimate cell type-specific parameters and simulate synthetic scRNA counts data using `scDesign3`.
   - `STEP3.R`: Simulate synthetic data from the results of `STEP1` and `STEP2`.
   - `STEP4.R`: Perform deconvolution on simulated data with `MuSiC`, `CARD`, `RCTD`, and `CibersortX`.
   - For BLADE, analysis has to be run in `ST_Python`, using `runBLADE.py`.

4. **Running BLADE**:
   - Ensure all R scripts have been run as described above.
   - Navigate to the `ST_Benchmark/spatial_BLADE` directory.
   - Activate the environment using `source BLADE/bin/activate`.
   - Use `requirements.txt` to install the necessary packages.
   - The framework for BLADE is present in `BLADE.py`.
   - Run `python runBLADE.py` to start the BLADE analysis.

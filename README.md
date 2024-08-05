### Processing

**Short title + project description:**
Benchmarking Spatial Transcriptomics Tools using Real and Simulated Data

**Main contact:**
[Contact Name - Email]

**Team:**
[List team members and roles]

#### Directories:

- **0_SoftwareEnvironment**: 
  Contains general and project-specific information about the computing environments used, including programming languages, integrated development environments, and package managers.

- **Data**: 
  Stores raw, meta, and pre-processed data relevant to the project.

- **ST_Benchmark_R**:  
  Contains all items related to the execution of frameworks based on R, including the preprocessing scripts and analysis files.

- **ST_Benchmark_Python**:  
  Holds the Python analysis scripts and framework files for BLADE. The `runBLADE.py` script is used to implement specific functionalities based on the functions defined.

#### Dependencies:
The Python scripts depend on the outputs from the R scripts. Specifically, the R scripts must be executed first to generate the necessary files for the Python-based BLADE deconvolution analysis.

#### How to execute the code:

1. **Running R Scripts**:
   - Navigate to the `ST_Benchmark_R` directory.
   - Execute `preprocess_ST.Rmd` followed by `preprocess_scRNA.Rmd` to generate the required files for deconvolution analysis.

2. **Running Python Scripts**:
   - Ensure all R scripts have been run as described above.
   - Navigate to the `ST_Benchmark_Python` directory.
   - Run `python runBLADE.py` to start the BLADE analysis.

**Note**: More details on parameters and settings will be added later.

### GitHub
The Processing directory should be synchronized to the GitHub repository excluding data and results, which reside in the FSS. Use the provided `.gitignore` templates to exclude specific subdirectories and files from the repository.

- The repository name and URL are listed in the `github.txt` file for reference by the FSS Navigator.


**Software Environment**

* **Operating System(s) / Version(s) Used During Development (and Testing):**  
  The software was developed and tested on Linux and macOS environments. It should also be compatible with Windows, but this has not been explicitly tested.

* **Specific Hardware Requirements:**  
 The software was tested on RAM: 132 Gb CPUs: 20 – Intel(R) Xeon(R) CPU E5-2690 v2 @ 3.00GHz. 
* **Software Environment (e.g., Programming Language + Version):**  
  - Python 3.8 or later

**Conceptual Description of Methodology**

This project uses the BLADE framework to perform deconvolution analyses on spatial transcriptomics data. The primary analysis is conducted using the `Framework_Iterative()` function in the `BLADE.py` script. The methodology includes selecting genes specifically for deconvolution experiments using the `autogene.py` script, which implements the method described in the paper with DOI: 10.1016/j.cels.2021.05.006.

**Random Seed:**  
A random seed is used for reproducibility in the `autogene.py` script, ensuring consistent gene selection across runs.

**Description of Manual Steps:**

1. Follow the instructions in the main README file to start the environment.
2. Download additional dependencies using the `requirements.txt` file.
3. Execute `autogene.py` to select genes for deconvolution. The selected genes are stored in `/spatial_BLADE/Results`.
4. Run `runBLADE.py` to execute the main BLADE analysis. Ensure that input files are placed in `/Processed_BLADE` within the `/Data` directory.

**Key Scripts and Their Functions:**

- **`BLADE.py`:** Contains the entire framework for BLADE analysis. Utilizes the `Framework_Iterative()` function for executing the analysis.

- **`autogene.py`:** Runs a gene selection tool specifically designed for deconvolution experiments. Implementation follows the documentation at [AutoGenes Documentation](https://autogenes.readthedocs.io/en/latest/getting_started.html).

- **`runBLADE.py`:** Main script to run BLADE analysis. Automatically determines the project root directory and connects input files from `/Processed_BLADE` in the `/Data` directory. The script uses spatial transcriptomics (ST) data, single-cell RNA-seq (scRNA-seq) mean and standard deviations, and the expectation matrix as input. The output is stored in the `/Results` directory in `/ST_benchmark`.

  - **Hyperparameters:** These are predefined for the analysis within the `run_BLADE()` function and can be modified as needed. The function takes input of the necessary files as specified.

**Intermediate Results and Output:**

Intermediate results and the hierarchical output are stored in the `/Results` directory. These outputs include detailed results from the BLADE analysis, allowing for inspection of different layers of the analysis.

This README provides an overview of the computation framework and essential steps for conducting the analyses within this project.
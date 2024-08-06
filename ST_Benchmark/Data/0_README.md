### Data Description

<u>Datasets</u>

**Dataset 1: MOB scRNA Data**  
- **Directory:** `scrna_MOB`
- **Description:** Mouse Olfactory Bulb (MOB) data taken from GSE121891.
- **Raw:** Contains raw count files for scRNA-seq.
- **Meta:** Contains metadata files for cell type annotations.
- **Data Owner:** GSE121891
- **License:** Refer to the GEO database for license details.

**Dataset 2: MOB Spatial Transcriptomics Data**  
- **Directory:** `scRNA_spatial`
- **Description:** Spatial transcriptomics data of the Mouse Olfactory Bulb.
- **Raw:** Raw counts file taken from [Spatial Research Resource](https://www.spatialresearch.org/resources-published-datasets/doi-10-1126science-aaf2403/).
- **Meta:** Contains ground truth annotations for the spots in the MOB ST data, taken from the authors of CARD (DOI: 10.1038/s41587-022-01273-7).
- **Image:** Contains the H&E section #12 for corresponding MOB ST data, also taken from [Spatial Research Resource](https://www.spatialresearch.org/resources-published-datasets/doi-10-1126science-aaf2403/).
- **Data Owner:** Authors of the study published in Science with DOI: 10.1126/science.aaf2403
- **License:** Refer to the original publication for license details.

**Dataset 3: Processed BLADE Data**  
- **Directory:** `Processed_BLADE`
- **Description:** Contains all the processed ST and scRNA data and expectation matrix files for BLADE execution.
- **Purpose:** For use in the `\spatialBLADE` directory for BLADE analysis.

**Dataset 4: PDAC Deconvolution Data**  
- **Directory:** `PDAC_deconvolution`
- **Description:** Contains scRNA and ST counts data, x and y coordinate information for ST spots, and ground truth annotations done by pathologists to highlight regions of interest.
- **Source:** The data is taken from [GSE111672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672) and downloaded from the [CARD tutorial](https://yma-lab.github.io/CARD/documentation/03_data.html).
- **Processing:** The data is processed using the `Deconvolution_PDAC.Rmd` script.


**Licenses:**

- For all data obtained from public databases, the respective data usage licenses should be consulted and complied with as per the source requirements.

This README file provides a structured overview of the datasets used in the project, their sources, processing steps, and how they relate to the overall workflow.

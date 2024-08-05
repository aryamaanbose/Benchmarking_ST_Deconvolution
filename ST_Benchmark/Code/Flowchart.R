install.packages("DiagrammeR")
library(DiagrammeR)

library(DiagrammeR)

grViz("
digraph flowchart {
  
  node [fontname = Helvetica, shape = rectangle, style=filled, fillcolor=white, fontsize = 35]
  edge [fontname = Helvetica, penwidth=2]

  // Define nodes and clusters
  subgraph cluster_cvi {
    label = 'COLLAPSED VARIATIONAL INFERENCE'
    fontsize = 40
    fontname = Helvetica
    style = 'filled, bold'
    fontcolor = black
    bgcolor = '#FFFFE0'
    
    SetPriors [label=' INITIALIZE HIDDEN VARIABLES AND HYPERPARAMETERS \\n SET PRIORS : Dirichlet for Cell Type Fractions per Sample, Log-Normal informed by Mean and SD of scRNA for Per Gene per Cell Type \\n TRUE POSTERIOR: Observed bulk RNA-seq data', fillcolor='#ADD8E6']
    ComputeELBO [label='COMPUTE ELBO\\n Measure of how well the variational distribution approximates the true posterior\\n Log-likelihood of the observed data (given the variational parameters) + KL divergence between the variational distribution and the prior', fillcolor='#90EE90']
    Optimization [label='OPTIMIZATION\\nL-BFGS Algorithm', shape=circle, fillcolor='#FFB6C1']
    UpdateParams [label='UPDATE VARIATIONAL PARAMETERS\\nAdjust to maximize ELBO', fillcolor='#90EE90']
    ConvergeCheck [label='Check Convergence Criteria for ELBO', shape=diamond, fillcolor='#FFA07A']
    OutputResults [label='Output Estimated Fractions per Cell Type', fillcolor='orange']
    Continue [label='Continue Inference', style='dotted, bold']
    

    SetPriors -> Optimization
    Optimization -> ComputeELBO
    ComputeELBO -> UpdateParams
    UpdateParams -> ConvergeCheck
    ConvergeCheck -> OutputResults [label='Yes', color='green', fontsize = 30]
    ConvergeCheck -> Continue [label='No', color='red', fontsize = 30]
    Continue -> Optimization
  }

  // Main Flow
  Input [label='INPUT DATA: Bulk RNA-seq and scRNA-seq (Annotated for Cell Types)', fillcolor='#D3D3D3']
  Preprocess [label='Quality Control, Log-Normalization and Filtration for Informative Genes', fillcolor='#D3D3D3']
  BulkRNA [label='Bulk RNA-seq Data', shape=diamond, fillcolor='#FFDAB9']
  scRNA [label='scRNA-seq Data', shape=diamond, fillcolor='#FFDAB9']
  EstimateParameters [label='Estimate Parameters\\nPer Gene Per Cell Type\\n(Mean and SD)', fillcolor='#ADD8E6']
  CreateExpression [label='Create Expression Profiles', fillcolor='#ADD8E6']
  
  // Connections
  Input -> Preprocess
  Preprocess -> BulkRNA
  Preprocess -> scRNA
  scRNA -> EstimateParameters
  EstimateParameters -> CreateExpression
  BulkRNA -> SetPriors
  CreateExpression -> SetPriors
}
")



grViz("
digraph flowchart { // Layout direction: left to right
  node [fontname = 'Helvetica', fontsize=10, shape = rectangle, style=filled, fillcolor=white, width=1.5, height=0.5]
  edge [fontname = 'Helvetica', penwidth=1.5]

  // Main Inputs
  
  Input2 [label='Input:  Real Tissue Coordinates' fillcolor='#D3D3D3']
  Input1 [label='Input: Reference scRNA data', fillcolor='#D3D3D3']

  // Proportion Matrix Design Flow
  subgraph cluster_proportion {
    label = 'Proportion Matrix Design'
    color = '#B3E2F2'
    style = 'filled, rounded'
    bgcolor = '#E6F4F1'


    DefineDominant [label='Define Dominant Cell Type Proportions', fillcolor='#B3E2F2']
    DefineNeighbourhood [label='Define Neighbourhood and Radius', fillcolor='#B3E2F2']
    Noise [label='Noise Parameter', fillcolor='#B3E2F2']
    SimulateProp [label='Simulate Realistic Proportion Design', fillcolor='#B3E2F2']
    
    DefineDominant -> DefineNeighbourhood -> Noise -> SimulateProp
  }

  // scDesign3 Model Simulation Framework
  subgraph cluster_scDesign3 {
    label = 'scDesign3 Model Simulation Framework'
    color = '#FFCC99'
    style = 'filled, rounded'
    bgcolor = '#FFEEDD'

    EstimateParams [label='Estimate Parameters from scRNA Data', fillcolor='#FFCC99']
    CopulaFramework [label='Use Copula Framework to Fit Joint Distribution of gene expression counts', fillcolor='#FFCC99']
    SyntheticData [label='Generate Synthetic Data per Cell Type', fillcolor='#FFCC99']
    
    EstimateParams -> CopulaFramework -> SyntheticData
  }

  // Simulating Synthetic Spatial Data
  subgraph cluster_spatial {
    label = 'Simulating Synthetic Spatial Data'
    color = '#90EE90'
    style = 'filled, rounded'
    bgcolor = '#E6FFCC'

    SpatialDataGeneration [label='Simulating Synthetic Spatial Data', fillcolor='#90EE90']
    SampleCells [label='Sample Fixed Number of Cells from Counts', fillcolor='#B7E9AD']
    MultiplyExpression [label='Multiply Expression by Proportions', fillcolor='#B7E9AD']
    AverageExpressions [label='Average Expressions for Each Spot', fillcolor='#B7E9AD']
    
    SpatialDataGeneration -> SampleCells -> MultiplyExpression -> AverageExpressions
  }

  // Connections between main input and clusters
  Input2 -> DefineDominant
  Input1 -> EstimateParams  
  SimulateProp -> SpatialDataGeneration
  SyntheticData -> SpatialDataGeneration
}
")







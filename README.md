# SingleCellAnalysis

### 1. SingleCellScanpyAutomatedDataAnalysisWorkflow.ipynb
This is an end-to-end automated Python 3.7 script using Scanpy for analysis of Single Cell transcriptomic data. It runs in multiple steps for performing discrete tasks in the pipeline. The script is very-well commented and explained with step-by-step instructions on inputs and outputs and functions performed. It can also use input as excel file for filetring based on QC metrics for multiple datasets (Example File included: QCmetricsFilter.xlsx). The section for automated integration is under progress and will be updated soon.

### 2. Seurat3_to_Monocle3ObjectConversion_forGeneModuleAnalysis.R
This script converts Seurat object to a Mnocle 3 cds object and does Gene Module Analysis.

### 3. AutomatedSeuratObjectIntegration.R
This script does automated integration of multiple Seurat objects in a memory efficient manner not overloading the environment. It also save the objects and creates DotPlots for multiple files containing marker list for cell types.

# LRP1-a-novel-IP-biomarker

# LRP1 IPF biomarker project – scRNAseq and bulk RNA analysis

This repository contains the R scripts used for the analysis of:
- Single-cell RNA-seq (Human Lung Cell Atlas)
- Fibroblast sub-clustering
- LRP1 expression analysis
- CytoTRACE
- CellChat signaling analysis
- Pathway scoring (PROGENy / GSVA)

## Requirements
- R version 4.4+
- Seurat 4.4+
- monocle3
- CellChat
- GSVA
- tidyverse

## Running the analysis
Each script in the `scripts/` directory corresponds to one step of the workflow:
- `01_QC_filtering.R` : Data loading and QC
- `02_Integration.R` : SCTransform & Integration
- `03_UMAP_clustering.R` : Clustering and annotation
- `04_LRP1_analysis.R` : LRP1 expression and visualization

## Contact
For questions, please contact: [Makoto Yamamoto]

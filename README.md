# longcovid_metaanalysis

This repository contains analysis code and results from our meta-analysis of bulk RNA-seq datasets profiling immune cells in COVID-19 and Long COVID.

## Folder Structure

- `results/`  
  Contains output files from differential expression analysis, Gene Ontology (GO) enrichment, KEGG pathway enrichment, and statistical tests across various immune cell comparisons.

- `code/`  
  Includes shell and R scripts for:
  - Preprocessing RNA-seq data (SRA download, trimming, STAR alignment)
  - Differential expression analysis using DESeq2
  - Gene annotation via biomaRt
  - Functional enrichment analysis (GO and KEGG)
  - Overlap and regression-based statistical tests

  **Some scripts are templates** and require user modification (e.g., file paths or contrast groups), while others are final versions used in the manuscript.

## Usage

1. Clone or download the repository.
2. Modify scripts under `code/` as needed for your dataset and environment.
3. Run scripts in the following general order:
   - `SRA_download_template.sh`
   - `FastQC_Trimming_template.sh`
   - `STAR_alignment.sh`
   - `DESeq2_differential_expression_template.R`
   - `annotate_DEG_results_with_biomaRt.R`
   - `GO_KEGG_functional_analyses.R`
   - Optional statistical analyses (e.g., `fisher_overlap_DEG_template.R`, `multiple_linear_regression.R`)

> For full details on each script and its input/output, see the associated README files inside the `code/` folder.

## Citation

If you use or reference this repository, please cite:

**Maison DP, Khadka VS, Mohd-Ibrahim I, Peluso MJ, Henrich TJ, Deng Y, Gerschenson M.**  
*Peripheral Immune Progression to Long COVID is Associated with Mitochondrial Gene Transcription: a Meta-Analysis.* 2025.

---

Statistical Tests
=================
Author: Maison et al., 2025

This directory contains R scripts for statistical evaluation of:
- Overlap between differential gene expression (DEG) results from different datasets.
- Covariate-adjusted expression effects via multiple linear regression.

--------------------------------------------------------------------------------
Contents:
---------

1. fisher_overlap_DEG_template.R
   - Description:
     Performs overlap analysis between two DEG sets using:
     * Fisher's exact test
     * Hypergeometric test
   - Input:
     * Two CSV files from DESeq2, each containing columns: ENSEMBL, log2FoldChange, padj
   - Output:
     * Overlap statistics and significance
     * Contingency table for overlap
     * CSV file classifying genes into "mono only", "PBMC only", or "both"
     * Results printed to console and optionally saved

2. multiple_linear_regression.R
   - Description:
     Performs multiple linear regression on each gene to assess the relationship between gene expression and sample-level covariates.
   - Input:
     * Expression data file (rows = genes, columns = samples)
     * Metadata file (must include sample_id column matching expression column names)
   - Output:
     * Data table with regression coefficients, p-values, and adjusted p-values for age and condition (or other covariates)
     * Optionally filtered tables of significant genes

--------------------------------------------------------------------------------
Usage Instructions:
-------------------
1. Open each script in R or RStudio.
2. Modify the file paths and filenames at the top of each script:
   * Input file paths
   * Output file names
   * Variable names (e.g., for conditions, covariates)
3. Run the script.

Ensure for `multiple_linear_regression.R`:
- The metadata and expression matrix use matching sample IDs.
- Categorical variables are converted to factors.

Ensure for `fisher_overlap_DEG_template.R`:
- Both input files use the same gene identifier format (e.g., ENSEMBL).
- The DEGs are defined using padj < 0.05 by default; modify if needed.

--------------------------------------------------------------------------------
Dependencies:
-------------
These scripts use base R packages only. No external libraries are required.

--------------------------------------------------------------------------------
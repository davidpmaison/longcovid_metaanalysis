===================================================================
README: Differential Expression Analysis Workflow
===================================================================
Author: Maison et al., 2025

This directory provides a complete workflow for differential gene expression analysis using R and DESeq2. The pipeline proceeds from raw STAR `ReadsPerGene.out.tab` outputs to fully annotated and visualized results.

-----------------
WORKFLOW OVERVIEW
-----------------

This analysis is performed in three main steps:

1.  **Differential Expression (DE) Analysis**
    * **Script:** `DESeq2_differential_expression_template.R`
    * **Purpose:** Performs the core DESeq2 analysis from STAR-aligned count files.
    * **Input:** Multiple STAR-generated `ReadsPerGene.out.tab` files (column 2 = unstranded) and a metadata file.
    * **Output:**
        - DE results CSV file with log2 fold changes, p-values, and adjusted p-values.
        - Normalized count matrix (optional).
        - Shrunk log2FC using `apeglm`.

2.  **Annotation of DE Results**
    * **Script:** `annotate_DEG_results_with_biomaRt.R`
    * **Purpose:** Adds gene annotations such as gene symbols to the DESeq2 result file using Ensembl BioMart.
    * **Input:** DE results file from Step 1.
    * **Output:** Annotated CSV file with gene symbol and description columns.

3.  **Visualization**
    * **Script:** `volcano_plot_template.R`
    * **Purpose:** Generates a volcano plot for visualizing differential expression.
    * **Input:** Annotated DE results file from Step 2.
    * **Output:** Volcano plot image (PDF or PNG).

-----------------
HOW TO USE
-----------------

1. **Prepare Input Files:**
    * STAR alignment outputs named `ReadsPerGene.out.tab` from each sample.
    * Metadata CSV file with:
        - `sample_id`: must match STAR filenames
        - `condition`: group labels (e.g., CR, LC)
        - Additional columns (e.g., batch, age) optional

2. **Configure DESeq2 Script:**
    * Open `DESeq2_differential_expression_template.R`
    * Edit the following:
        - `dir`: working directory path
        - `ref`: reference group name (e.g., "CR")
        - `files`: pattern-matched STAR output
        - Output file name

3. **Run Analysis:**
    a. Run `DESeq2_differential_expression_template.R`
    b. Run `annotate_DEG_results_with_biomaRt.R` on the resulting DE CSV
    c. Run `volcano_plot_template.R` on the annotated CSV to generate plots

-----------------
NOTES
-----------------
- The DESeq2 script extracts column 2 from STAR outputs, which corresponds to unstranded counts. Modify if your experiment used stranded protocols.
- The volcano plot uses default cutoffs (`log2FC > 1`, `p < 0.05`) but can be customized.
- Gene annotations are pulled from the Ensembl database via `biomaRt`, so an active internet connection is required during Step 2.

===================================================================
# Gene Annotation Script (Template)
# Author: David Maison
# Purpose: Annotate ENSEMBL IDs with gene symbols and gene types using biomaRt.

# Load necessary libraries
library(biomaRt)
library(readr)

# === USER INPUT REQUIRED ===
# Set the path to the directory containing your DEG results CSV file
setwd("path/to/your/results_directory")  # <-- Update this path

# Specify the input CSV file containing differential expression results
input_file <- "your_DEG_results.csv"  # <-- Update with your actual filename

# Specify the output file name for the annotated results
output_file <- "annotated_DEG_results.csv"  # <-- Desired output file name

# === STEP 1: Load DEG Results ===
deg_data <- read_csv(input_file)

# Ensure the first column contains ENSEMBL gene IDs; rename it for clarity
names(deg_data)[1] <- "ENSEMBL"

# === STEP 2: Query Ensembl BioMart for Annotations ===
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract gene symbols and gene types
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = deg_data$ENSEMBL,
  mart = mart
)

# === STEP 3: Merge Annotations ===
deg_data$GeneID <- gene_info$external_gene_name[match(deg_data$ENSEMBL, gene_info$ensembl_gene_id)]

# Optionally fill missing GeneID with gene biotype
for (i in seq_len(nrow(deg_data))) {
  if (is.na(deg_data$GeneID[i]) || deg_data$GeneID[i] == "") {
    deg_data$GeneID[i] <- gene_info$gene_biotype[match(deg_data$ENSEMBL[i], gene_info$ensembl_gene_id)]
  }
}

# === STEP 4: Reorder Columns (Optional) ===
deg_data <- deg_data[, c("ENSEMBL", "GeneID", names(deg_data)[-1])]

# === STEP 5: Save Annotated Results ===
write_csv(deg_data, output_file)

# === END ===
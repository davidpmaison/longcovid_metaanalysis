# Volcano Plot Template for DESeq2 Differential Expression Results
# Author: David Maison
# Usage: Customize file paths, comparison name, and sample sizes below.

# Load required package
library(EnhancedVolcano)

# === User Input Section ===

# Set the working directory where your differential expression result CSV file is located
setwd("YOUR/WORKING/DIRECTORY/HERE")

# Specify the input CSV file with differential expression results (must include GeneID, log2FoldChange, and padj columns)
deg_file <- "YOUR_DEG_RESULTS_FILE.csv"

# Specify output PDF file name for the volcano plot
output_pdf <- "volcano_plot_DESeq2.pdf"

# Specify plot title and subtitle
plot_title <- "Comparison: CONDITION_A vs CONDITION_B"
plot_subtitle <- "CELL TYPE: CONDITION_A n = XX; CONDITION_B n = YY"

# === End of User Input Section ===

# Load differential expression results
res <- read.csv(deg_file)

# Generate volcano plot
pdf(file = output_pdf, width = 7, height = 7)
EnhancedVolcano(
  res,
  lab = res$GeneID,  # Gene name labels
  x = 'log2FoldChange',
  y = 'padj',
  title = plot_title,
  subtitle = plot_subtitle,
  pCutoff = 0.05,
  FCcutoff = 0,
  pointSize = 2,
  labSize = 3.0,
  boxedLabels = FALSE,
  ylab = bquote("-Log"[10] * "(adjusted p-value)"),
  col = c("#052049", "#178CCB", "#A238BA", "#C42882"),
  labCol = "#2E2872",
  cutoffLineCol = "#052049",
  hlineCol = "#052049",
  vlineCol = "#052049",
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  borderColour = "#052049"
)
dev.off()
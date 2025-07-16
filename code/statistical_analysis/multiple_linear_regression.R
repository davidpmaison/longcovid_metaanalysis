# ============================================================
# run_age_covariate_regressions_and_plot.R
# ============================================================
# Author: Maison et al., 2025
# ============================================================
# This script:
# 1. Loads normalized counts and metadata (containing age and condition)
# 2. Subsets to a comprehensive gene list and overlapping samples
# 3. Transposes counts and merges with metadata
# 4. Runs a multiple linear model per gene to test the effect of age
#    while controlling for condition (expression ~ age + condition)
# 5. Extracts statistics for the 'age' term from each model
# 6. Adjusts P-values (Benjamini–Hochberg FDR)
# 7. Saves results and exports a PDF plot of the age estimates ± SE,
#    coloring bars by significance and outlining mitochondrial genes.
# ============================================================

# ------------------------------------------------------------
# 0) LOAD LIBRARIES
# ------------------------------------------------------------
library(data.table)
library(purrr)
library(ggplot2)

# ------------------------------------------------------------
# 1) SETUP: directory & gene lists
# ------------------------------------------------------------
# IMPORTANT: Set this to your working directory
dir <- "/Users/davidmaison/Documents/Maison_MetaAnalysis_SourceComparisons/Age_Final/"
setwd(dir)

# Gene IDs to be analyzed
genes_to_test <- c(
  # A comprehensive list of Ensembl gene IDs
  "ENSG00000002726", "ENSG00000137819", "ENSG00000146205",
  "ENSG00000147647", "ENSG00000156049", "ENSG00000164695",
  "ENSG00000174944", "ENSG00000198075", "ENSG00000211626",
  "ENSG00000223551", "ENSG00000229344", "ENSG00000253957",
  "ENSG00000288349", "ENSG00000054938", "ENSG00000070729",
  "ENSG00000102239", "ENSG00000104332", "ENSG00000151224",
  "ENSG00000171403", "ENSG00000171747", "ENSG00000196634",
  "ENSG00000198153", "ENSG00000198727", "ENSG00000198938",
  "ENSG00000200227", "ENSG00000200520", "ENSG00000200959",
  "ENSG00000210077", "ENSG00000214097", "ENSG00000214321",
  "ENSG00000215874", "ENSG00000219163", "ENSG00000222932",
  "ENSG00000225284", "ENSG00000225349", "ENSG00000225630",
  "ENSG00000226348", "ENSG00000226629", "ENSG00000232117",
  "ENSG00000238420", "ENSG00000243295", "ENSG00000243777",
  "ENSG00000243883", "ENSG00000248449", "ENSG00000250403",
  "ENSG00000250669", "ENSG00000252072", "ENSG00000253317",
  "ENSG00000253961", "ENSG00000254132", "ENSG00000254331",
  "ENSG00000255231", "ENSG00000255307", "ENSG00000256148",
  "ENSG00000259502", "ENSG00000259686", "ENSG00000259866",
  "ENSG00000260737", "ENSG00000262902", "ENSG00000266304",
  "ENSG00000270225", "ENSG00000272986", "ENSG00000273481",
  "ENSG00000273882", "ENSG00000274452", "ENSG00000277769",
  "ENSG00000278599", "ENSG00000278992", "ENSG00000279174",
  "ENSG00000279263", "ENSG00000279423", "ENSG00000279993",
  "ENSG00000285648", "ENSG00000285972", "ENSG00000287189",
  "ENSG00000289209", "ENSG00000290597", "ENSG00000290866"
)

# Mitochondrial Ensembl IDs to highlight
mito_genes <- c("ENSG00000198727", "ENSG00000198938", "ENSG00000198153")

# ------------------------------------------------------------
# 2) LOAD DATA
# ------------------------------------------------------------
full_counts_dt <- fread("normalized_counts.csv")
metadata_dt    <- fread("metadata_subset_with_age_CR_LC.csv")

# Ensure metadata has required columns for the multiple regression
if (!all(c("sample_id", "age", "condition") %in% colnames(metadata_dt))) {
  stop("Metadata must contain 'sample_id', 'age', and 'condition' columns.")
}

# ------------------------------------------------------------
# 3) CHECK SAMPLE OVERLAP
# ------------------------------------------------------------
sample_overlap <- intersect(metadata_dt$sample_id, colnames(full_counts_dt))
if (length(sample_overlap) == 0) stop("No overlapping samples found.")
cat("Found", length(sample_overlap), "overlapping samples.\n")

# ------------------------------------------------------------
# 4) SUBSET COUNTS
# ------------------------------------------------------------
counts_subset_dt <- full_counts_dt[
  GeneID %in% genes_to_test,
  c("GeneID", sample_overlap),
  with = FALSE
]

# ------------------------------------------------------------
# 5) TRANSPOSE
# ------------------------------------------------------------
counts_long_dt <- melt(
  counts_subset_dt,
  id.vars       = "GeneID",
  variable.name = "sample_id",
  value.name    = "expr"
)
counts_wide_dt <- dcast(
  counts_long_dt,
  sample_id ~ GeneID,
  value.var = "expr"
)

# ------------------------------------------------------------
# 6) MERGE METADATA
# ------------------------------------------------------------
regression_data <- merge(
  metadata_dt[, .(sample_id, age, condition)],
  counts_wide_dt,
  by = "sample_id",
  all = FALSE
)
if (nrow(regression_data) == 0) stop("No data after merging metadata and counts.")

# Ensure age is numeric and condition is a factor for modeling
regression_data[, age := as.numeric(age)]
regression_data[, condition := as.factor(condition)]

# ------------------------------------------------------------
# 7) RUN MULTIPLE REGRESSION MODELS (expression ~ age + condition)
# ------------------------------------------------------------
gene_cols <- intersect(genes_to_test, colnames(regression_data))

all_results <- map_dfr(gene_cols, function(gene) {
  # Create a subset for the model, removing any rows with missing data
  sub_data <- regression_data[
    !is.na(age) & !is.na(condition) & !is.na(get(gene)),
    .(age, condition, expr = get(gene))
  ]
  
  # Check for sufficient data: at least 3 samples and more than one condition level
  if (nrow(sub_data) < 3 || length(unique(sub_data$condition)) < 2) {
    return(NULL)
  }
  
  # Fit the multiple linear regression model
  fit <- lm(expr ~ age + condition, data = sub_data)
  s   <- summary(fit)
  coef_table <- s$coefficients
  
  # Ensure the 'age' term is present in the model output before proceeding
  if (!"age" %in% rownames(coef_table)) {
    return(NULL)
  }
  
  # Extract statistics specifically for the 'age' coefficient
  data.frame(
    GeneID          = gene,
    N               = nrow(sub_data),
    Estimate_Age    = coef_table["age", "Estimate"],
    Std.Error_Age   = coef_table["age", "Std. Error"],
    T.Value_Age     = coef_table["age", "t value"],
    P.Value_Age     = coef_table["age", "Pr(>|t|)"],
    Model.R.Squared = s$r.squared,
    stringsAsFactors = FALSE
  )
})

# ------------------------------------------------------------
# 8) MULTIPLE TEST CORRECTION
# ------------------------------------------------------------
# Adjust p-values from the 'age' term in the models
all_results$adj.P.Value <- p.adjust(all_results$P.Value_Age, method = "BH")
all_results$Signif      <- all_results$adj.P.Value < 0.05
all_results$IsMito      <- all_results$GeneID %in% mito_genes

# Reorder genes for plotting based on their adjusted P-value for age
all_results$GeneID <- factor(
  all_results$GeneID,
  levels = all_results$GeneID[order(all_results$adj.P.Value)]
)

# ------------------------------------------------------------
# 9) SAVE RESULTS
# ------------------------------------------------------------
fwrite(
  all_results[order(all_results$adj.P.Value), ],
  file.path(dir, "age_as_covariate_results.csv")
)

# ------------------------------------------------------------
# 10) EXPORT PLOT TO PDF
# ------------------------------------------------------------
pdf(file = file.path(dir, "age_as_covariate_plot.pdf"),
    width = 12, height = 6)
print(
  ggplot(all_results, aes(x = GeneID, y = Estimate_Age, fill = Signif, color = IsMito)) +
    geom_col(size = 0.5) +
    geom_errorbar(aes(
      ymin = Estimate_Age - Std.Error_Age,
      ymax = Estimate_Age + Std.Error_Age
    ), width = 0.2, size = 0.5) +
    scale_fill_manual(
      values = c(`TRUE` = "steelblue", `FALSE` = "grey80"),
      name   = "adj.P < 0.05 (Age Term)"
    ) +
    scale_color_manual(
      values = c(`TRUE` = "firebrick", `FALSE` = NA),
      guide  = guide_legend(override.aes = list(fill = NA)),
      name   = "Mitochondrial"
    ) +
    theme_bw() +
    theme(
      axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top"
    ) +
    labs(
      x     = "Ensembl Gene ID",
      y     = "Slope for Age (Controlling for Condition)",
      title = "Effect of Age on Gene Expression, Controlling for Condition"
    )
)
dev.off()


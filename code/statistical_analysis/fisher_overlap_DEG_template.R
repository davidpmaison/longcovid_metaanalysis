# ---------------------------
# Template: Fisher's Exact and Hypergeometric Overlap Test for DEGs
# ---------------------------
# Author: Maison et al., 2025

# 1. Set working directory (update to your path)
setwd("/your/path/to/analysis_directory")

# 2. Load DEG results (update filenames)
mono_all <- read.csv("mononuclear_DEG_results.csv", row.names = 1, stringsAsFactors = FALSE)
pbmc_all <- read.csv("pbmc_DEG_results.csv", row.names = 1, stringsAsFactors = FALSE)

# 3. Filter to significant genes (adjust threshold as needed)
mono_sig <- subset(mono_all, !is.na(padj) & padj < 0.05)
pbmc_sig <- subset(pbmc_all, !is.na(padj) & padj < 0.05)

# 4. Count DEGs and extract gene sets
n_mono_DEG <- nrow(mono_sig)
n_pbmc_DEG <- nrow(pbmc_sig)

genes_mono <- rownames(mono_sig)
genes_pbmc <- rownames(pbmc_sig)

# 5. Determine overlaps
common_genes  <- intersect(genes_mono, genes_pbmc)
n_common      <- length(common_genes)
n_only_mono   <- length(setdiff(genes_mono, genes_pbmc))
n_only_pbmc   <- length(setdiff(genes_pbmc, genes_mono))

# 6. Define background gene universe
all_genes <- union(rownames(mono_all), rownames(pbmc_all))
N <- length(all_genes)

# 7. Hypergeometric test
p_hyper <- phyper(
  q = n_common - 1,
  m = n_mono_DEG,
  n = N - n_mono_DEG,
  k = n_pbmc_DEG,
  lower.tail = FALSE
)

# 8. Fisher's exact test
contingency <- matrix(
  c(n_common,
    n_pbmc_DEG - n_common,
    n_mono_DEG - n_common,
    N - n_common - (n_pbmc_DEG - n_common) - (n_mono_DEG - n_common)),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    PBMC = c("DE", "not_DE"),
    Mononuclear = c("DE", "not_DE")
  )
)

fisher_res <- fisher.test(contingency, alternative = "greater")

# 9. Print summary
cat("DE genes (padj<0.05):\n",
    "  Mononuclear:", n_mono_DEG, "\n",
    "  PBMC:       ", n_pbmc_DEG, "\n\n",
    "Overlap:\n",
    "  Shared:          ", n_common, "\n",
    "  Only Mononuclear:", n_only_mono, "\n",
    "  Only PBMC:       ", n_only_pbmc, "\n\n",
    "Background genes:", N, "\n\n",
    "Hypergeometric p-value: ", signif(p_hyper, 3), "\n",
    "Fisher’s test:\n"
)
print(fisher_res)

# 10. Save overlap classification
write.csv(
  data.frame(
    gene = all_genes,
    in_mononuc = all_genes %in% genes_mono,
    in_pbmc = all_genes %in% genes_pbmc
  ),
  file = "gene_overlap_classification.csv",
  row.names = FALSE,
  quote = FALSE
)
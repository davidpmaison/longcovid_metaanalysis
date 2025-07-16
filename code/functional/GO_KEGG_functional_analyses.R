# ===================================================================
# Master Script for Automated Downstream RNA-seq Analysis
# Author: Maison et al., 2025
#
# This script automates the following for multiple DESeq2 result files:
# 1. Filters for significant DEGs (padj < 0.05).
# 2. Annotation of ENSEMBL IDs to Gene Symbols for significant DEGs.
# 3. Generation of a Volcano Plot and a GO Enrichment Plot.
# 4. Combines plots into a single summary figure.
# 5. KEGG Pathway Enrichment Analysis.
# ===================================================================


# --- 1. Load All Required Libraries ---
# Ensure all these packages are installed before running
# install.packages(c("biomaRt", "readr", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "pathview", "ggplot2", "patchwork"))
library(biomaRt)
library(readr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)
library(ggplot2)
library(patchwork) # For combining plots


# --- 2. Configuration ---
# Set the main working directory
main_dir <- "/Users/davidmaison/Documents/Maison_MetaAnalysis_SourceComparisons/Figures/mono_new"
setwd(main_dir)

# Define sample sizes for Mononuclear cells based on the provided image
sample_sizes <- list(
  AC = 466,
  NC = 360,
  CR = 81,
  LC = 23
)

# Define full names for plot titles
condition_full_names <- list(
  AC = "Acute COVID (AC)",
  NC = "No known COVID (NC)",
  CR = "COVID Recovered (CR)",
  LC = "Long COVID (LC)"
)


# --- 3. Find and Order DEG Result Files ---
# Find all CSV files that match the DESeq2 output pattern
all_deg_files <- list.files(path = main_dir, pattern = "^deg_results_.*\\.csv$", full.names = TRUE)

# Define the desired processing order to match your filenames
desired_order <- c(
  "deg_results_AC_vs_NC.csv",
  "deg_results_CR_vs_NC.csv",
  "deg_results_LC_vs_NC.csv",
  "deg_results_CR_vs_AC.csv",
  "deg_results_LC_vs_AC.csv",
  "deg_results_LC_vs_CR.csv"
)

# Match the found files to the desired order
# This ensures the loop runs in the specified sequence
deg_files <- all_deg_files[match(desired_order, basename(all_deg_files))]
deg_files <- deg_files[!is.na(deg_files)] # Remove any not found

cat("Found", length(deg_files), "DEG result files to process in the specified order.\n\n")


# --- 4. Main Processing Loop ---
# Loop through each DEG file and perform the full analysis pipeline
for (input_file in deg_files) {
  
  # --- A. SETUP FOR CURRENT COMPARISON ---
  
  # Extract comparison name from filename (e.g., "AC_vs_NC")
  comparison_name <- sub("^deg_results_", "", basename(input_file))
  comparison_name <- sub("\\.csv$", "", comparison_name)
  
  groups <- strsplit(comparison_name, "_vs_")[[1]]
  group1_short <- groups[1] # This is the Treatment
  group2_short <- groups[2] # This is the Control
  
  cat("=========================================================\n")
  cat("Processing comparison:", comparison_name, "\n")
  cat("=========================================================\n\n")
  
  # Create a dedicated output directory for this comparison
  output_dir <- file.path(main_dir, comparison_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # --- B. PRE-FILTERING FOR SIGNIFICANT GENES ---
  
  cat("Step 1: Filtering for significant DEGs (padj < 0.05)...\n")
  
  # Read DEG results
  deg_results_full <- read.csv(input_file)
  colnames(deg_results_full)[1] <- "ENSEMBL"
  
  # Filter for significant genes first
  sig_genes_df <- subset(deg_results_full, !is.na(padj) & padj < 0.05)
  
  if (nrow(sig_genes_df) == 0) {
    cat("  -> No significant genes found for this comparison. Skipping.\n\n")
    next # Skip to the next file
  }
  cat("  -> Found", nrow(sig_genes_df), "significant genes.\n\n")
  
  # --- C. ANNOTATION ---
  
  cat("Step 2: Annotating significant results...\n")
  
  # Annotate with biomaRt
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
    filters = "ensembl_gene_id",
    values = sig_genes_df$ENSEMBL,
    mart = mart
  )
  
  # Merge annotations
  annotated_results <- merge(sig_genes_df, gene_info, by.x = "ENSEMBL", by.y = "ensembl_gene_id", all.x = TRUE)
  # Use gene biotype if gene symbol is missing
  annotated_results$external_gene_name[is.na(annotated_results$external_gene_name)] <- annotated_results$gene_biotype[is.na(annotated_results$external_gene_name)]
  
  # Save annotated results
  annotated_file_path <- file.path(output_dir, paste0("mononuclear_", comparison_name, "_annotated_significant_results.csv"))
  write.csv(annotated_results, annotated_file_path, row.names = FALSE)
  cat("  -> Saved annotated significant results to:", basename(annotated_file_path), "\n\n")
  
  # --- D. VOLCANO PLOT GENERATION (CORRECTED) ---
  
  cat("Step 3: Generating Volcano Plot object...\n")
  
  # Generate plot title with full names, using the order from the filename
  plot_title <- paste(condition_full_names[[group1_short]], "vs", condition_full_names[[group2_short]])
  plot_subtitle <- paste0("Mononuclear: ", group1_short, " n=", sample_sizes[[group1_short]], "; ", group2_short, " n=", sample_sizes[[group2_short]])
  
  # The data is already in the correct orientation (Treatment vs Control), so no inversion is needed.
  plot_data <- annotated_results
  cat("  -> Plotting contrast as", group1_short, "vs", group2_short, "\n")
  
  # Generate and store the volcano plot as an object
  p_volcano <- EnhancedVolcano(
    plot_data, # Use the original data
    lab = plot_data$external_gene_name,
    x = 'log2FoldChange',
    y = 'padj',
    title = plot_title,
    subtitle = plot_subtitle,
    pCutoff = 0.05,
    FCcutoff = 0,
    pointSize = 2.0,
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
  
  
  # --- E. FUNCTIONAL ENRICHMENT (GO & KEGG) ---
  
  cat("Step 4: Performing Functional Enrichment (GO & KEGG)...\n")
  
  # Prepare gene list for clusterProfiler from the already filtered data
  gene_symbols <- annotated_results$external_gene_name
  entrez_map <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Merge with original data to keep logFC values
  enrich_data <- merge(entrez_map, annotated_results, by.x = "SYMBOL", by.y = "external_gene_name")
  enrich_data <- enrich_data[!duplicated(enrich_data$ENTREZID), ]
  
  gene_list_fc <- enrich_data$log2FoldChange
  names(gene_list_fc) <- enrich_data$ENTREZID
  
  # --- GO ENRICHMENT ---
  
  go_ontologies <- c("BP", "MF", "CC")
  go_results_list <- list()
  for (ont in go_ontologies) {
    cat("  -> Running GO analysis for:", ont, "...\n")
    ego <- enrichGO(gene = names(gene_list_fc), OrgDb = org.Hs.eg.db, ont = ont, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    if (!is.null(ego) && nrow(ego) > 0) {
      go_results_list[[ont]] <- ego
      go_results_path <- file.path(output_dir, paste0("mononuclear_", comparison_name, "_", tolower(ont), "_go_results.csv"))
      write.csv(as.data.frame(ego), go_results_path, row.names = FALSE)
      cat("    -> Saved", ont, "results to:", basename(go_results_path), "\n")
    } else {
      cat("    -> No significant", ont, "terms found.\n")
    }
  }
  cat("\n")
  
  # --- COMBINED GO PLOT ---
  p_go <- NULL # Initialize GO plot as NULL
  if (length(go_results_list) > 0) {
    plot_df_list <- list()
    if (!is.null(go_results_list$BP) && nrow(go_results_list$BP) > 0) plot_df_list$BP <- head(as.data.frame(go_results_list$BP), 5)
    if (!is.null(go_results_list$MF) && nrow(go_results_list$MF) > 0) plot_df_list$MF <- head(as.data.frame(go_results_list$MF), 5)
    if (!is.null(go_results_list$CC) && nrow(go_results_list$CC) > 0) plot_df_list$CC <- head(as.data.frame(go_results_list$CC), 5)
    
    if (length(plot_df_list) > 0) {
      combined_df <- do.call(rbind, lapply(names(plot_df_list), function(name) {
        data.frame(Category = name, Term = plot_df_list[[name]]$Description, EnrichmentScore = -log10(plot_df_list[[name]]$p.adjust))
      }))
      
      p_go <- ggplot(combined_df, aes(x = reorder(Term, EnrichmentScore), y = EnrichmentScore, fill = Category)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        facet_wrap(~Category, scales = "free_y", ncol = 1) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.y = element_text(size = 10), 
          strip.text = element_text(size = 12, face = "bold"), 
          legend.position = "none",
          plot.title = element_blank() # Remove plot title
        ) +
        labs(
          y = "Enrichment Score (-log10 adjusted p-value)", 
          x = "GO Term"
        ) +
        scale_fill_manual(values = c("BP" = "#052049", "MF" = "#B8E6FA", "CC" = "#178CCB"))
    }
  }
  
  # --- F. SAVE ALL PLOTS ---
  cat("Step 5: Saving individual and combined plots...\n")
  
  # Save the individual volcano plot
  volcano_path <- file.path(output_dir, paste0("mononuclear_", comparison_name, "_volcano_plot_significant.pdf"))
  ggsave(volcano_path, plot = p_volcano, width = 7, height = 7, device = "pdf")
  cat("  -> Saved individual volcano plot to:", basename(volcano_path), "\n")
  
  # If a GO plot was created, save it individually and then create the combined figure
  if (!is.null(p_go)) {
    # Save the individual GO plot
    go_plot_path <- file.path(output_dir, paste0("mononuclear_", comparison_name, "_Top_GO_Terms_Plot.pdf"))
    ggsave(go_plot_path, plot = p_go, width = 12, height = 8, device = "pdf")
    cat("  -> Saved individual GO plot to:", basename(go_plot_path), "\n")
    
    # Combine plots using patchwork
    combined_figure <- p_volcano + p_go + plot_layout(widths = c(1, 1.2)) # Give GO plot slightly more space
    
    # Save the combined figure
    summary_figure_path <- file.path(output_dir, paste0("mononuclear_", comparison_name, "_summary_figure.pdf"))
    ggsave(summary_figure_path, plot = combined_figure, width = 18, height = 9, device = "pdf", limitsize = FALSE)
    cat("  -> Saved combined summary figure to:", basename(summary_figure_path), "\n\n")
  } else {
    cat("  -> No significant GO terms found, so no GO or summary plot was generated.\n\n")
  }
  
  
  # --- G. KEGG ENRICHMENT ---
  
  cat("Step 6: Running KEGG analysis...\n")
  kegg <- enrichKEGG(gene = names(gene_list_fc), organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)
  if (!is.null(kegg) && nrow(kegg) > 0) {
    kegg_results_path <- file.path(output_dir, paste0("mononuclear_", comparison_name, "_KEGG_results.csv"))
    write.csv(as.data.frame(kegg), kegg_results_path, row.names = FALSE)
    cat("    -> Saved KEGG results to:", basename(kegg_results_path), "\n")
    
    kegg_plot_path <- file.path(output_dir, paste0("mononuclear_", comparison_name, "_KEGG_barplot.pdf"))
    pdf(file = kegg_plot_path, width = 10, height = 8)
    print(barplot(kegg, showCategory = 15) + labs(title = paste("KEGG Pathways:", comparison_name)))
    dev.off()
    cat("    -> Saved KEGG bar plot to:", basename(kegg_plot_path), "\n")
    
  } else {
    cat("    -> No significant KEGG pathways found.\n")
  }
  
  cat("\nFinished processing", comparison_name, ".\n\n")
}

cat("=========================================================\n")
cat("All comparisons have been processed successfully.\n")
cat("=========================================================\n")


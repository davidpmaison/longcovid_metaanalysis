FUNCTIONAL ANALYSIS README (GO and KEGG Enrichment Combined)
Author: Maison et al., 2025

This folder contains a unified R script for functional enrichment analysis of differentially expressed genes (DEGs), combining:

1. Gene Ontology (GO) enrichment analysis
2. KEGG pathway enrichment analysis

================================================================================
FILE CONTENT:
- GO_KEGG_functional_analyses.R
  A single script that performs both GO and KEGG enrichment analyses, outputs result tables, and generates visualizations.

================================================================================
INPUT FILE FORMAT:
The input CSV file must include the following columns:
GeneID	log2FoldChange	padj

Example:
GeneID	log2FoldChange	padj  
H4C3	6.528554176	6.03E-89  
ZNHIT2	3.355050009	6.43E-55  
NICOL1	4.136140501	1.56E-54  

- `GeneID`: gene symbols  
- `log2FoldChange`: log fold change from DESeq2  
- `padj`: adjusted p-value from DESeq2  

================================================================================
SCRIPT FUNCTIONALITY:
- Converts `GeneID` to Entrez IDs using `org.Hs.eg.db`
- Filters for significant genes (`padj < 0.05`)
- Performs:
  - GO enrichment (Biological Process, Molecular Function, Cellular Component)
  - KEGG pathway enrichment
- Generates:
  - CSV tables of enrichment results for each ontology (GO) and KEGG
  - PDF bar plots for top GO and KEGG terms
  - Optional KEGG pathway diagrams using `pathview` (based on log2FoldChange)

================================================================================
OUTPUT:
- `GO_enrichment_BP.csv`, `GO_enrichment_MF.csv`, `GO_enrichment_CC.csv`
- `KEGG_enrichment_results.csv`
- Bar plots: `GO_enrichment_barplot.pdf`, `KEGG_enrichment_barplot.pdf`
- Optional: KEGG pathway diagrams (e.g., `hsa04110.pathview.png`)

================================================================================
REQUIRED R PACKAGES:
Ensure the following libraries are installed:
- `clusterProfiler`
- `org.Hs.eg.db`
- `enrichplot`
- `ggplot2`
- `pathview` (for pathway diagrams)

Install with:
```R
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "pathview"))
install.packages("ggplot2")
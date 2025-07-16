# ---------------------------------------------
# Viral RNA Analysis Template
# ---------------------------------------------
# USER: update the following file paths and parameters as needed
metadata_file   <- "path/to/metadata_with_classification.csv"  # Metadata CSV file
count_dir       <- "path/to/counts/"                          # Directory containing per-sample count files
output_prefix   <- "Fig"                                       # Prefix for output files

# Load required packages
library(dplyr)
library(ggstatsplot)
library(ggplot2)
library(patchwork)

# Set a consistent base theme and font sizes
theme_set(
  theme_bw(base_size = 16) +
    theme(
      plot.title     = element_text(size = 18, face = "bold"),
      axis.title     = element_text(size = 16),
      axis.text      = element_text(size = 14),
      legend.title   = element_text(size = 16),
      legend.text    = element_text(size = 14),
      plot.tag       = element_text(size = 18, face = "bold")  # Tags 2 pts larger
    )
)

# 1. Read the augmented metadata CSV
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# 2. Compute total and SARS-CoV-2 read sums for each sample
sums_df <- bind_rows(lapply(metadata$file_name, function(f) {
  # Construct full path to count file
  file_path <- file.path(count_dir, f)
  
  # Read STAR 'ReadsPerGene.out.tab' format: Gene, Counts, (ignored)
  counts <- read.table(
    file_path,
    col.names = c("Gene", "Counts", "Ignore1", "Ignore2"),
    stringsAsFactors = FALSE
  )
  # sum_total: sum of all counts
  sum_total <- sum(counts$Counts)
  # sum_sars: sum of SARS-CoV-2-specific rows (customize indices or gene IDs as needed)
  # Default here: rows 5-16; adjust per your file structure
  sum_sars  <- sum(counts$Counts[5:16])
  
  data.frame(sum_total = sum_total, sum_sars = sum_sars)
}))

# 3. Combine with metadata and compute normalized SARS ratio
df <- bind_cols(metadata, sums_df) %>%
  mutate(
    SARS = sum_sars / sum_total
  )

# 4A. Figure A: all cells by cell_classification
p1 <- ggbetweenstats(
  data                 = df,
  x                    = cell_classification,
  y                    = SARS,
  xlab                 = "Cell Classification",
  ylab                 = "SARS-CoV-2 RNA / Total RNA",
  title                = "SARS-CoV-2 / Total RNA by Cell Classification",
  pairwise.comparisons = TRUE,
  p.adjust.method      = "holm"
) +
  scale_color_manual(
    values = c("Granulocyte" = "#178CCB", "Mononuclear" = "#C42882"),
    guide  = FALSE
  )

# 4B. Figure B: mononuclear only, by condition
p2 <- df %>%
  filter(cell_classification == "Mononuclear") %>%
  ggbetweenstats(
    data                 = ., 
    x                    = condition,
    y                    = SARS,
    xlab                 = "Condition",
    ylab                 = "SARS-CoV-2 RNA / Total RNA",
    title                = "SARS-CoV-2 / Total RNA by Condition\n(Mononuclear)",
    pairwise.comparisons = TRUE,
    p.adjust.method      = "holm"
  ) +
  scale_color_manual(
    values = c("NC" = "#178CCB", "AC" = "#C42882", 
               "CR" = "#A238BA", "LC" = "#16A0AC"),
    guide  = FALSE
  )

# 4C. Figure C: mononuclear only, by cell_type
p3 <- df %>%
  filter(cell_classification == "Mononuclear") %>%
  ggbetweenstats(
    data                 = ., 
    x                    = cell_type,
    y                    = SARS,
    xlab                 = "Cell Type",
    ylab                 = "SARS-CoV-2 RNA / Total RNA",
    title                = "SARS-CoV-2 / Total RNA by Cell Type\n(Mononuclear)",
    pairwise.comparisons = TRUE,
    p.adjust.method      = "holm"
  ) +
  scale_color_manual(
    values = c(
      "Leukocytes"   = "#178CCB",
      "Monocytes"    = "#C42882",
      "PBMC"         = "#16A0AC",
      "G-MDSC"       = "#A238BA",
      "Treg"         = "#00BFC4",
      "memory CD4+"  = "#ECAD00",
      "CD4,CD8,B"    = "#F8766D",
      "CD4+ / CD8+"  = "#7CAE00"
    ),
    guide = FALSE
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

# 5. Display plots
print(p1)
print(p2)
print(p3)

# 6. Save each figure as a separate PDF
ggsave(paste0(output_prefix, "_CellClassification.pdf"), p1, width = 7, height = 6)
ggsave(paste0(output_prefix, "_Condition_Mononuclear.pdf"), p2, width = 7, height = 6)
ggsave(paste0(output_prefix, "_CellType_Mononuclear.pdf"), p3, width = 7, height = 10)

# 7. Combine all three into one tall, single-column PDF with tags A), B), C)
combined <- (p1 / p2 / p3) +
  plot_layout(ncol = 1, heights = c(1, 1, 1.5), tag_level = 'new') +
  plot_annotation(
    tag_levels = 'A',    # Start tags at A
    tag_prefix = '',     # No prefix
    tag_suffix = ')',    # Add closing parenthesis
    theme      = theme(plot.tag = element_text(size = 18, face = "bold"))
  ) &
  theme(plot.margin = margin(10, 10, 10, 10))

ggsave(paste0(output_prefix, "_Combined.pdf"), combined, width = 10, height = 24)

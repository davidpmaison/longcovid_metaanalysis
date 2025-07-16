# DESeq2 Differential Expression Analysis Template
# Author: Maison et al., 2025

# ------------------ Load Required Packages ------------------
library(DESeq2)
library(data.table)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ------------------ 1. SETUP: Define Paths ------------------
# Set working directory
dir <- "<INSERT/PROJECT/DIRECTORY>"
setwd(dir)

# Load metadata
sampleData <- fread(file.path(dir, "<INSERT_METADATA_FILENAME.csv>"))
rownames(sampleData) <- sampleData$sample_id

# Define condition factor
sampleData$condition <- as.factor(sampleData$condition)
sampleData$condition <- relevel(sampleData$condition, ref = "<INSERT_REFERENCE_GROUP>")

# ------------------ 2. LOAD: Import Count Files ------------------
files <- list.files(dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)

# Extract counts (column 2: unstranded)
countList <- lapply(files, function(file) {
  dat <- fread(file, header = FALSE)
  dat <- dat[, .(V1, V2)]
  setnames(dat, c("gene", basename(file)))
  return(dat)
})

# Merge all count files
mergedCounts <- Reduce(function(...) merge(..., by = "gene"), countList)
rownames(mergedCounts) <- mergedCounts$gene
countData <- as.matrix(mergedCounts[, -1, with = FALSE])

# Ensure column names match sampleData
colnames(countData) <- gsub(".*\\/|\\.ReadsPerGene.out.tab", "", colnames(countData))

# ------------------ 3. RUN DESeq2 ------------------
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleData,
                              design = ~ condition)
dds <- DESeq(dds)

# ------------------ 4. EXTRACT RESULTS ------------------
# Replace contrast groups as needed
res <- results(dds, contrast = c("condition", "<GROUP1>", "<GROUP2>"))
res <- lfcShrink(dds, coef = 2, res = res, type = "apeglm")

# Annotate with gene names
res$ENSEMBL <- rownames(res)
res$GeneID <- mapIds(org.Hs.eg.db,
                     keys = res$ENSEMBL,
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
res <- as.data.table(res)

# ------------------ 5. SAVE RESULTS ------------------
fwrite(res, file = file.path(dir, "DESeq2_results_<GROUP1>_vs_<GROUP2>.csv"))

# ------------------ 6. OPTIONAL: Volcano Plot ------------------
EnhancedVolcano(res,
                lab = res$GeneID,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = '<GROUP1> vs <GROUP2>',
                pCutoff = 0.05,
                FCcutoff = 1.0)

# ------------------ END ------------------
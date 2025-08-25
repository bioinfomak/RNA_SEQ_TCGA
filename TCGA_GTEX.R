library("TCGAbiolinks")
library("recount3")
library("survminer")
library("survival")
library("tidyverse")
library("DESeq2")
library("tidyverse")
library("maftools")
library("pheatmap")
library("DT")
library("SummarizedExperiment")
library("data.table")

# STEP 2: Create project directory structure
# Create main project folder
project_dir <- "~/TCGA_GTEx_Analysis"
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)

# Create subfolders
dir.create(file.path(project_dir, "data"), showWarnings = FALSE)
dir.create(file.path(project_dir, "results"), showWarnings = FALSE)
dir.create(file.path(project_dir, "plots"), showWarnings = FALSE)

# Set working directory
setwd(project_dir)

# STEP 3: Explore available TCGA datasets
#Get all TCGA cancer projects
all_projects <- getGDCprojects()

# Filter for TCGA projects (they start with "TCGA-")
tcga_projects <- all_projects[grepl("^TCGA-", all_projects$project_id), ]

# Show available cancer types
print(tcga_projects[, c("project_id", "name")])


# STEP 4: Check available samples for TCGA-BRCA
# Query what data is available for breast cancer
brca_query <- GDCquery(
  project = "TCGA-BRCA",                    # Breast cancer project
  data.category = "Transcriptome Profiling", # RNA-seq data
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"           # Count data (not FPKM)
)

# Get sample information
sample_info <- getResults(brca_query)

# Check sample types available
sample_types <- table(sample_info$sample_type)
print(sample_types)


# STEP 5: Download TCGA breast cancer data
# We want both tumor and normal samples
target_samples <- c("Primary Tumor", "Solid Tissue Normal")

# Create focused query
brca_query_filtered <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = target_samples
)

# Check how many samples we'll download
filtered_samples <- getResults(brca_query_filtered)
print(table(filtered_samples$sample_type))

# Download the data (this takes several minutes)
GDCdownload(brca_query_filtered)

# Process the downloaded data into R format
tcga_query_data <- GDCprepare(brca_query_filtered, summarizedExperiment = TRUE)
cat(paste("Downloaded", ncol(tcga_query_data), "samples with", nrow(tcga_query_data), "genes\n"))   

# Downloaded 1224 samples with 60660 genes

# STEP 6: Clean and examine TCGA data
# Get the expression matrix (genes x samples)
tcga_expression <- assay(tcga_query_data, "unstranded")  # Raw counts
tcga_samples <- colData(tcga_query_data)                 # Sample information
tcga_genes <- rowData(tcga_query_data)                   # Gene information

# Check sample distribution
sample_dist <- table(tcga_samples$sample_type)
print(sample_dist)

# Keep only protein-coding genes (reduces noise)
if ("gene_type" %in% colnames(tcga_genes)) {
  protein_coding <- tcga_genes$gene_type == "protein_coding"
  tcga_data_clean <- tcga_query_data[protein_coding, ]
} else {
  tcga_data_clean <- tcga_query_data
}

# STEP 7: Save TCGA data in multiple formats
# Save processed data
saveRDS(tcga_data_clean, "data/tcga_brca_processed.rds")

# Also save as text files for use in other programs
tcga_expression_clean <- assay(tcga_data_clean, "unstranded")
tcga_samples_clean <- as.data.frame(colData(tcga_data_clean))
tcga_genes_clean <- as.data.frame(rowData(tcga_data_clean))

# Save expression matrix with gene names
expression_table <- data.table(
  gene_id = rownames(tcga_expression_clean),
  gene_name = tcga_genes_clean$gene_name,
  tcga_expression_clean
)
fwrite(expression_table, "data/tcga_brca_expression.tsv", sep = "\t")

# Save sample information
tcga_samples_clean$sample_id <- rownames(tcga_samples_clean)

# Convert list columns to character (more concise approach)
tcga_samples_clean[, sapply(tcga_samples_clean, is.list)] <- lapply(tcga_samples_clean[, sapply(tcga_samples_clean, is.list)], as.character)

fwrite(tcga_samples_clean, "data/tcga_brca_samples.tsv", sep = "\t")



# STEP 8: Explore available GTEx datasets
# Check what projects are available in recount3
available_data <- available_projects()

# Look for GTEx data
gtex_data <- available_data[available_data$file_source == "gtex", ]
print(gtex_data[, c("project", "organism", "file_source")])

# STEP 9: Get GTEx sample information
# Create GTEx data object (this is lightweight - just gets metadata)
gtex_rse <- create_rse_manual(
  project = "BREAST",
  project_home = "data_sources/gtex", 
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

# Get sample information
gtex_samples <- colData(gtex_rse)

# See what tissue types are available
table(gtex_samples$gtex.smtsd)

# Breast - Mammary Tissue   482
# STEP 10: Download GTEx expression data
# Get the actual expression data (this triggers download)
gtex_expression <- assay(gtex_rse, "raw_counts")
gtex_samples_filtered <- colData(gtex_rse)
gtex_genes <- rowData(gtex_rse)

# STEP 11: Process GTEx data
# Filter for protein-coding genes (to match TCGA)
gtex_gene_df <- as.data.frame(gtex_genes)
if ("gene_type" %in% colnames(gtex_gene_df)) {
  protein_coding <- gtex_gene_df$gene_type == "protein_coding"
  gtex_clean <- gtex_rse[protein_coding, ]
  cat(paste("Kept", sum(protein_coding), "protein-coding genes\n"))
} else {
  gtex_clean <- gtex_rse
}

# Focus on breast tissue samples for comparison with TCGA
breast_mask <- gtex_clean$gtex.smtsd == "Breast - Mammary Tissue"
gtex_breast <- gtex_clean[, breast_mask]

# Save processed GTEx data
saveRDS(gtex_clean, "data/gtex_processed.rds")
saveRDS(gtex_breast, "data/gtex_breast.rds")

# STEP 12: Save GTEx data
# Get clean data for saving
gtex_expression_clean <- assay(gtex_clean, "raw_counts")  # Use raw_counts assay
gtex_samples_clean <- as.data.frame(colData(gtex_clean))
gtex_genes_clean <- as.data.frame(rowData(gtex_clean))

# Save expression matrix with gene information
gtex_table <- data.table(
  gene_id = rownames(gtex_expression_clean),
  gene_name = gtex_genes_clean$gene_name,
  gtex_expression_clean
)
fwrite(gtex_table, "data/gtex_breast_expression.tsv", sep = "\t")

# Save sample information (handle list columns)
gtex_samples_clean$sample_id <- rownames(gtex_samples_clean)
gtex_samples_clean[, sapply(gtex_samples_clean, is.list)] <- 
  lapply(gtex_samples_clean[, sapply(gtex_samples_clean, is.list)], as.character)
fwrite(gtex_samples_clean, "data/gtex_breast_samples.tsv", sep = "\t")


# STEP 13: Prepare data for differential expression analysis
# Load our processed data
tcga_data <- readRDS("data/tcga_brca_processed.rds")
gtex_breast <- readRDS("data/gtex_breast.rds")

# Get expression matrices
tcga_expression <- assay(tcga_data, "unstranded")  # TCGA counts
gtex_expression <- assay(gtex_breast, "raw_counts")    # GTEx counts

# Remove genes that are not expressed in both datasets
tcga_expression <- tcga_expression[rowSums(tcga_expression > 0) >= ncol(tcga_expression), ]

gtex_expression <- gtex_expression[rowSums(gtex_expression > 0) >= ncol(gtex_expression), ]

# Get sample information
tcga_samples <- colData(tcga_data)
gtex_samples <- colData(gtex_breast)

# Get gene information
tcga_gene <- rowData(tcga_data)[rownames(tcga_expression), ]
gtex_gene <- rowData(gtex_breast)[rownames(gtex_expression), ]

rownames(tcga_gene) <- gsub("\\.\\d+$", "", rownames(tcga_gene))
rownames(gtex_gene) <- gsub("\\.\\d+$", "", rownames(gtex_gene))

# Find common genes between datasets
common_genes <- intersect(rownames(tcga_gene), rownames(gtex_gene))
head(common_genes)
tcga_gene_common <- tcga_gene[common_genes, ]
gtex_gene_common <- gtex_gene[common_genes, ]

# Keep only common genes
rownames(tcga_expression) <- gsub("\\.\\d+$", "", rownames(tcga_expression))
rownames(gtex_expression) <- gsub("\\.\\d+$", "", rownames(gtex_expression))

tcga_expression <- tcga_expression[rownames(tcga_gene_common), ]
gtex_expression <- gtex_expression[rownames(gtex_gene_common), ]

# Focus on TCGA tumor samples only
tumor_samples <- tcga_samples$sample_type == "Primary Tumor"
tcga_tumor_expression <- tcga_expression[, tumor_samples]


# STEP 14: Combine datasets and address batch effects
# Combine expression matrices
combined_expression <- cbind(tcga_tumor_expression, gtex_expression)

# Create combined sample information
# Use data.table for efficient data manipulation
tcga_info <- data.table(
  sample_id = colnames(tcga_tumor_expression),
  dataset = "TCGA",
  condition = "Tumor"
)

gtex_info <- data.table(
  sample_id = colnames(gtex_expression),
  dataset = "GTEx", 
  condition = "Normal"
)

combined_samples <- rbind(tcga_info, gtex_info)

print(table(combined_samples$dataset, combined_samples$condition))

# Filter out very low expressed genes
# Keep genes with at least 10 counts in at least 10% of samples
min_samples <- round(0.1 * ncol(combined_expression))
keep_genes <- rowSums(combined_expression > 10) >= min_samples
filtered_expression <- combined_expression[keep_genes, ]


# STEP 15: Differential expression analysis with limma

# Load required packages
library("limma")
library("edgeR")

# Create DGEList object for normalization
dge <- DGEList(counts = filtered_expression)

# Add gene symbols to the DGEList object
dge$genes <- data.frame(
  gene_id = rownames(filtered_expression),
  gene_name = tcga_gene_common$gene_name,
  gene_type = tcga_gene_common$gene_type
)

# Add sample group information
dge$samples$condition <- factor(combined_samples$condition)
dge$samples$dataset <- factor(combined_samples$dataset)

# Calculate normalization factors (accounts for library size differences)
dge <- calcNormFactors(dge, method = "TMM")

# Create design matrix
# Include both condition (tumor vs normal) and dataset (TCGA vs GTEx)
design <- model.matrix(~ condition, data = combined_samples)
colnames(design) <- c("Intercept", "Tumor_vs_Normal")

# Apply voom transformation (converts counts to log2-cpm with weights)
v <- voom(dge, design, plot = FALSE)

# Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)  # Apply empirical Bayes smoothing

# Extract results for tumor vs normal comparison
results <- topTable(
  fit,
  coef = "Tumor_vs_Normal",  # This coefficient represents tumor vs normal
  number = Inf,              # Get all genes
  adjust.method = "BH"       # Multiple testing correction
)

fwrite(results, "DEGs_Tumor_vs_Normal.tsv", sep = "\t")

# STEP 16: Create visualization plots

# Load plotting library
library(ggplot2)

# 1. Volcano plot
volcano_data <- data.table(
  gene_name = results$gene_name,
  logFC = results$logFC,
  neg_log10_pval = -log10(results$adj.P.Val),
  significant = results$adj.P.Val < 0.05 & abs(results$logFC) > 1
)

ggplot(volcano_data, aes(x = logFC, y = neg_log10_pval)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Volcano Plot: Tumor vs Normal Breast Tissue",
    x = "Log2 Fold Change (Tumor vs Normal)",
    y = "-Log10(Adjusted P-value)",
    color = "Significant"
  ) +
  theme_minimal()

ggsave("plots/volcano_plot.png", volcano_plot, width = 8, height = 6, dpi = 300)

# 2. MA plot (shows relationship between expression level and fold change)
ggplot(results, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed")) +
  labs(
    title = "MA Plot: Average Expression vs Fold Change",
    x = "Average Expression (log2 CPM)",
    y = "Log2 Fold Change",
    color = "FDR < 0.05"
  ) +
  theme_minimal()

ggsave("plots/ma_plot.png", ma_plot, width = 8, height = 6, dpi = 300)


# STEP 17: Assess potential batch effects

# Let's check some "housekeeping" genes that shouldn't change much
housekeeping_genes <- c("ACTB", "GAPDH", "RPL13A", "SDHA")

# Find these genes in our results
hk_results <- results[results$gene_name %in% housekeeping_genes, ]

cat("Housekeeping gene results (should have small changes):\n")
print(hk_results[, c("gene_name", "logFC", "adj.P.Val")])

# If housekeeping genes show large changes, this suggests batch effects
large_hk_changes <- abs(hk_results$logFC) > 0.5
if (any(large_hk_changes)) {
  cat("<img draggable="false" role="img" class="emoji" alt="⚠️" src="https://s.w.org/images/core/emoji/16.0.1/svg/26a0.svg">  Warning: Some housekeeping genes show large changes\n")
  cat("This suggests possible batch effects between TCGA and GTEx\n")
} else {
  cat("✓ Housekeeping genes look stable\n")
}


# Focus on genes with very large effect sizes and high significance
high_confidence <- results[
  abs(results$logFC) > 2 &           # Large effect size
    results$adj.P.Val < 0.001 &       # Very significant
    results$AveExpr > 5,              # Well-expressed genes
]

cat(paste("High-confidence DE genes:", nrow(high_confidence), "\n"))

# Dichotomize risk
high_risk <- risk_score > median(risk_score)

# KM plot
survfit_obj <- survfit(Surv(time = surv_data$time, event = surv_data$status) ~ high_risk)
plot(survfit_obj, col = c("blue", "red"), lwd = 2)
legend("bottomleft", legend = c("Low Risk", "High Risk"), col = c("blue", "red"), lwd = 2)


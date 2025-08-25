library("TCGAbiolinks")
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

#get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-BRCA')


query_exp <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling', # parameter enforced by GDCquery
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts', # Count data (not FPKM)
                       data.type = "Gene Expression Quantification",
                       access = 'open')

GDCdownload(query_exp)
data_exp <- GDCprepare(query_exp)

# Get sample information
sample_info <- getResults(query_exp)

# Check sample types available
sample_types <- table(sample_info$sample_type)
print(sample_types)

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

# Get the expression matrix (genes x samples)
tcga_expression <- assay(tcga_query_data, "unstranded")  # Raw counts
tcga_samples <- colData(tcga_query_data)                 # Sample information
tcga_genes <- rowData(tcga_query_data)                   # Gene information

# Check sample distribution
sample_dist <- table(tcga_samples$sample_type)
print(sample_dist)

#       Primary Tumor 1111
# Solid Tissue Normal  113

# Keep only protein-coding genes (reduces noise)
if ("gene_type" %in% colnames(tcga_genes)) {
  protein_coding <- tcga_genes$gene_type == "protein_coding"
  tcga_data_clean <- tcga_query_data[protein_coding, ]
} else {
  tcga_data_clean <- tcga_query_data
}

# Save processed data
saveRDS(tcga_data_clean, "data/tcga_brca_processed.rds")

# Save the processed data as RDS (R's native format)
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

dds <- DESeqDataSetFromMatrix(countData = tcga_expression_clean,
                              colData = tcga_samples_clean,
                              design = ~ 1)


# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10


dds <- dds[keep,]
dds <- DESeq(dds)

res1 <- results( dds )
res1

idx <- which.min(res1$pvalue)
counts(dds)[idx, ]

# Assuming 'dds' is your DESeqDataSet object
# Apply the variance stabilizing transformation
# vsd now contains the variance-stabilized dataset
# vst to vsd
vsd <- vst(dds, blind=FALSE)
tcga_matrix_vst <- assay(vsd)
tcga_matrix_vst[1:10,1:10]


idx3 <- tcga_expression_clean %>% 
  assay %>% 
  rowVars(var = 'gene_name') %>% 
  order(decreasing = TRUE) %>% 
  head(5)

# plot
pheatmap(assay(vsd)[idx1,])

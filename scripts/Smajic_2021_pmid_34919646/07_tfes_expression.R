# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(DESeq2)
library(edgeR)

# Set working directory
setwd("/scratch/ben/rnaseq/")

# Define analysis parameters ----------------------------------------------

# Define model formula 
# model_formula <- "~ 1"
model_formula <- "~ disease_status:disease_status_member + cell_type + disease_status"

# Define if protein-coding genes only
protein_coding_only <- T

# Define FDR threshold
fdr <- 0.05

# Load metadata -----------------------------------------------------------

# Load gene metadata
gene_metadata <- readRDS("ref_seq/feature_metadata/gene_metadata.rds")

# Load and tidy sample metadata
sample_metadata <- read_csv("billy/data/input/Smajic_2021_pmid_34919646/sample_metadata_pb.rds") %>% 
  mutate(disease_status_member = str_replace(donor_id, "[A-Z]+", "m"))

# Construct count matrix --------------------------------------------------

# Load count matrix
count_matrix <- read.delim("billy/data/input/Smajic_2021_pmid_34919646/raw_count_matrix_pb.txt",
                     check.names = F)

# Subset to protein-coding genes only
if(protein_coding_only) {
  protein_coding_genes <- gene_metadata[gene_metadata$gene_type == "protein_coding", "ensembl_gene_id_version"][[1]]
  count_matrix <- count_matrix[rownames(count_matrix) %in% protein_coding_genes,]
} 

# Set up design matrix ----------------------------------------------------

# Construct design matrix
design <- model.matrix(object = formula(model_formula), 
                       data = sample_metadata)

# Make the matrix full-rank
# => Remove columns with all 0 from design matrix
design <- design[, which(colSums(design) != 0), drop = F]

# DESeq2 workflow ---------------------------------------------------------

# Construct DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_metadata,
                              design = design)

# dds <- estimateSizeFactors(dds)
# 
# norm_counts2 <- counts(dds, normalized = TRUE)

# Perform DGE
# => Using Wald test
dds_wald <- DESeq(object = dds, 
                  test = "Wald",
                  parallel = T,
                  BPPARAM = BiocParallel::MulticoreParam(workers = 8))

# => Using LRT
dds_lrt <- DESeq(object = DESeqDataSetFromMatrix(countData = count_matrix,
                                                 colData = sample_metadata,
                                                 design = design), 
                 test = "LRT",
                 reduced = design[, which(colnames(design) != "disease_statusPD")],
                 parallel = T,
                 BPPARAM = BiocParallel::MulticoreParam(workers = 8))

# Extract results table for each test
for(i in c("wald", "lrt")) {
  
  # Construct res
  assign(x = paste0("res_", i),
         value = results(object = eval(sym(paste0("dds_", i))),
                         alpha = fdr,
                         name = "disease_statusPD"))
  
  # Construct res_shrunken_lfc
  assign(x = paste0("res_", i, "_shrunken_lfc"),
         value = lfcShrink(dds = eval(sym(paste0("dds_", i))),
                           coef = "disease_statusPD",
                           type = "apeglm",
                           parallel = T,
                           BPPARAM = BiocParallel::MulticoreParam(workers = 8)))
  
  # Construct res_df
  assign(x = paste0("res_", i, "_df"),
         value = as_tibble(eval(sym(paste0("res_", i))),
                           rownames = "ensembl_gene_id_version") %>% 
           full_join(as_tibble(eval(sym(paste0("res_", i, "_shrunken_lfc"))),
                               rownames = "ensembl_gene_id_version") %>%
                       dplyr::rename(shrunken_log2FoldChange = log2FoldChange,
                                     shrunken_lfcSE = lfcSE) %>% 
                       select(ensembl_gene_id_version, shrunken_log2FoldChange, shrunken_lfcSE),
                     by = join_by(ensembl_gene_id_version)) %>% 
           mutate(sign_log2fc_times_minus_log10pvalue = sign(log2FoldChange) * -log(x = pvalue, base = 10),
                  log2fc = log2FoldChange) %>% 
           left_join(gene_metadata, by = join_by(ensembl_gene_id_version)))
  
  # Construct ddeg
  assign(x = paste0("ddeg_", i),
         value = eval(sym(paste0("res_", i, "_df"))) %>% 
           filter(padj < fdr) %>% 
           arrange(desc(log2FoldChange)))
  
}

# Extract log2 of normalized counts + 1
log_norm_counts <- as_tibble(log2(counts(dds_wald, normalized = TRUE) + 1), rownames = NA)
norm_counts <- counts(dds_wald, normalized = TRUE)

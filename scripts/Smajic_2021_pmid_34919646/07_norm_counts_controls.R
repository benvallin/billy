# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(DESeq2)
library(edgeR)

# Set working directory
setwd("/scratch/ben/rnaseq/")

# Define analysis parameters ----------------------------------------------

# Define model formula 
model_formula <- "~ 1"

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

dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)

norm_counts <- as_tibble(counts(dds, normalized = TRUE), rownames = "ensembl_gene_id_version") %>% 
  left_join(unique(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id","gene_name")]), 
            by = join_by(ensembl_gene_id_version)) %>% 
  pivot_longer(cols = unique(sample_metadata$donor_id_cell_type),
               names_to = "donor_id_cell_type",
               values_to = "norm_count") %>% 
  left_join(sample_metadata, by = join_by(donor_id_cell_type))

ggplot(data = norm_counts %>% 
         filter(gene_name %in% c("TFEB", "TFE3")),
       mapping = aes(x = cell_type, y = norm_count, fill = disease_status_member, shape = disease_status)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  scale_shape_manual(values = c(22, 21)) +
  facet_wrap(~ gene_name) +
  ggpubr::theme_pubr(legend = "right") +
  ggpubr::rotate_x_text(angle = 45) 

mean_norm_counts <- as_tibble(counts(dds, normalized = TRUE), rownames = "ensembl_gene_id_version") %>% 
  left_join(unique(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id","gene_name")]), 
            by = join_by(ensembl_gene_id_version))  %>% 
  mutate(mean_norm_count = rowMeans(.[, sample_metadata$sample_name]),
         sd_norm_count = rowSds(as.matrix(.[, sample_metadata$sample_name])),
         sem_norm_count = sd_norm_count / length(sample_metadata$sample_name)) 

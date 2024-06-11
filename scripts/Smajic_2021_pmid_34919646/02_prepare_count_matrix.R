# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(Matrix)
library(DropletUtils)

# Set working directory
setwd("/scratch/ben/rnaseq/")

# Construct count matrix --------------------------------------------------

# Load raw UMI count matrix
matrix <- as.matrix(fread("public_datasets/Smajic_2021_pmid_34919646/supplementary_files/IPDCO_hg_midbrain_UMI.tsv"))

# Set row names to ensembl gene IDs 
rownames(matrix) <- read_tsv("public_datasets/Smajic_2021_pmid_34919646/supplementary_files/IPDCO_hg_midbrain_genes.tsv") %>% 
  arrange(row) %>% 
  pull(gene)

# Convert to sparse matrix
matrix <- Matrix(matrix, sparse = TRUE)

# Write count matrix to file
write10xCounts(path = "billy/data/input/Smajic_2021_pmid_34919646/counts.h5",
               x = matrix,
               barcodes = colnames(matrix),
               gene.id = rownames(matrix),
               overwrite = T)


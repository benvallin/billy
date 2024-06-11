# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd("/scratch/ben/rnaseq/")


# Combinations of donor ID and cell type ----------------------------------

# Construct sample metadata
sample_metadata <- read_csv("billy/data/input/Smajic_2021_pmid_34919646/sample_metadata_pb.csv") %>% 
  mutate(donor_id_cell_type = paste0(donor_id, "_", cell_type),
         batch = "Smajic") %>% 
  select(batch, donor_id_cell_type, broad_cell_type, cell_type, donor_id:sn_neuron_loss, psbulk_n_cells, psbulk_counts)

write_csv(sample_metadata, "billy/data/input/Smajic_2021_pmid_34919646/sample_metadata_pb.rds")

# Construct count matrix
count_matrix <- read_csv(paste0("billy/data/input/Smajic_2021_pmid_34919646/raw_count_matrix_pb.csv"))%>% 
  column_to_rownames(var = "donor_id_cell_type")

count_matrix <- round(t(count_matrix))

identical(colnames(count_matrix), sample_metadata[["donor_id_cell_type"]])

write.table(count_matrix, "billy/data/input/Smajic_2021_pmid_34919646/raw_count_matrix_pb.txt", sep = "\t")


# Combinations of donor ID and broad cell type ----------------------------

# Construct sample metadata
sample_metadata <- read_csv("billy/data/input/Smajic_2021_pmid_34919646/sample_metadata_pb_broad_cell_type.csv") %>% 
  mutate(donor_id_broad_cell_type = paste0(donor_id, "_", broad_cell_type),
         batch = "Smajic") %>% 
  select(batch, donor_id_broad_cell_type, broad_cell_type, donor_id:sn_neuron_loss, psbulk_n_cells, psbulk_counts)

write_csv(sample_metadata, "billy/data/input/Smajic_2021_pmid_34919646/sample_metadata_pb_broad_cell_type.rds")

# Construct count matrix
count_matrix <- read_csv(paste0("billy/data/input/Smajic_2021_pmid_34919646/raw_count_matrix_pb_broad_cell_type.csv"))%>% 
  column_to_rownames(var = "donor_id_broad_cell_type")

count_matrix <- round(t(count_matrix))

identical(colnames(count_matrix), sample_metadata[["donor_id_broad_cell_type"]])

write.table(count_matrix, "billy/data/input/Smajic_2021_pmid_34919646/raw_count_matrix_pb_broad_cell_type.txt", sep = "\t")



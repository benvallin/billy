# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd("/scratch/ben/rnaseq/")

# Construct sample metadata -----------------------------------------------

# Load donor metadata
donor_metadata <- read_csv("public_datasets/Smajic_2021_pmid_34919646/supplementary_files/donor_metadata.csv")

# Construct sample metadata
sample_metadata <- read_tsv("public_datasets/Smajic_2021_pmid_34919646/supplementary_files/IPDCO_hg_midbrain_cell.tsv") %>% 
  select(barcode, 
         cell_type = cell_ontology,
         donor_id = patient) %>% 
  mutate(donor_id = str_replace(donor_id, "C", "CTRL")) %>% 
  left_join(donor_metadata, by = join_by(donor_id))

write_csv(sample_metadata, "billy/data/input/Smajic_2021_pmid_34919646/sample_metadata.csv")

# Construct tissue metadata -----------------------------------------------

# tissue_metadata <- read_delim(file = "public_datasets/Smajic_2021_pmid_34919646/SraRunTable.txt",
#                               delim = ",") %>%
#   select(run_id = Run,
#          sample_id = `Sample Name`,
#          disease_status = disease_state,
#          donor_sex = sex,
#          age_at_death,
#          tissue,
#          molecule_subtype,
#          LibraryLayout,
#          AvgSpotLen,
#          Instrument) %>%
#   mutate(donor_sex = ifelse(donor_sex == "female", "F", "M"),
#          disease_status = ifelse(disease_status == "Healthy control", "CTRL", "IPD"),
#          tissue = tolower(tissue))
# 
# write_csv(tissue_metadata, "billy/data/input/Smajic_2021_pmid_34919646/tissue_metadata.csv")
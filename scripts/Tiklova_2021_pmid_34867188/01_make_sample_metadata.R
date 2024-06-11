# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(GEOquery)

# Set working directory
setwd("/scratch/ben/rnaseq/")

# Construct sample metadata -----------------------------------------------

# Load and tidy donor metadata
donor_metadata <- read_csv("public_datasets/Tiklova_2021_pmid_34867188/supplementary_files/donor_metadata.csv") %>% 
  mutate(disease_status = case_when(str_detect(donor_id, "^C-.*$") ~ "CTRL",
                                    str_detect(donor_id, "^PD-.*$") ~ "PD",
                                    str_detect(donor_id, "^ILBD-.*$") ~ "ILBD",
                                    TRUE ~ NA_character_))

# Retrieve sample ID - sample name pairs
GSE182622_series_matrix <- GEOquery::getGEO(filename = "public_datasets/Tiklova_2021_pmid_34867188/supplementary_files/GSE182622_series_matrix.txt.gz")

sample_ids <- GSE182622_series_matrix@phenoData@data %>% 
  select(sample_id = geo_accession,
         sample_name = title) %>% 
  mutate(sample_name = str_replace(sample_name, "-I-", "-ILBD-"),
         donor_id = ifelse(str_detect(sample_name, "DA-"),
                           str_remove(sample_name, "DA-") %>% 
                             str_remove("-[1-9]{1,}-[A-Z]{1}.*"),
                           sample_name))

# Construct sample metadata
sample_metadata <- read_delim(file = "public_datasets/Tiklova_2021_pmid_34867188/SraRunTable.txt",
                              delim = ",") %>%
  select(run_id = Run,
         sample_id = `Sample Name`,
         species = Organism,
         tissue = source_name,
         cell_type,
         LibraryLayout,
         AvgSpotLen,
         Instrument) %>% 
  left_join(sample_ids, by = join_by(sample_id)) %>% 
  left_join(donor_metadata, by = join_by(donor_id)) %>% 
  mutate(batch = "Tiklova",
         files = setNames(
           object = map_chr(.x = run_id,
                            .f = ~ list.files(path = paste0("public_datasets/Tiklova_2021_pmid_34867188/03_salmon_selective_alignment_trim_galore/", .x),
                                              pattern = "quant.sf",
                                              full.names = TRUE)),
           nm = run_id)
  ) %>%
  select(files, batch, donor_id, disease_status, age_at_death:n_samples, sample_name, everything())

write_rds(sample_metadata, "billy/data/input/Tiklova_2021_pmid_34867188/sample_metadata.rds")



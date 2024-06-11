# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd("/scratch/ben/rnaseq/")

# Construct sample metadata -----------------------------------------------

sample_metadata <- read_delim(file = "public_datasets/Monzon_Sandoval_2020_pmid_32167609/SraRunTable.txt",
                              delim = ",") %>%
  select(run_id = Run,
         sample_id = `Sample Name`,
         donor_id = Origin,
         donor_sex = sex,
         donor_age = AGE,
         tier,
         rin = RIN,
         LibraryLayout,
         AvgSpotLen,
         Instrument) %>%
  mutate(batch = "Monzon-Sandoval",
         species = "human",
         disease_status = "CTRL",
         donor_sex = ifelse(donor_sex == "female", "F", "M"),
         tissue = "SNpc",
         precise_tissue = paste0(tier, " tier ", tissue),
         sample_name = paste0(donor_id, " ", tier, " ", tissue, " DaNs"),
         files = setNames(
           object = map_chr(.x = run_id,
                            .f = ~ list.files(path = paste0("public_datasets/Monzon_Sandoval_2020_pmid_32167609/03_salmon_selective_alignment_trim_galore/", .x),
                                              pattern = "quant.sf",
                                              full.names = TRUE)), 
           nm = run_id)
  ) %>% 
  select(files, batch, run_id:donor_id, sample_name, 
         species, disease_status, donor_sex, donor_age, 
         tissue, tier, precise_tissue, rin, LibraryLayout:Instrument) 

write_rds(sample_metadata, "billy/data/input/Monzon_Sandoval_2020_pmid_32167609/sample_metadata.rds")

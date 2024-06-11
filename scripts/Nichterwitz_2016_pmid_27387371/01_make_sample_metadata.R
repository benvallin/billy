# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd("/scratch/ben/rnaseq/")

# Construct sample metadata -----------------------------------------------

sample_metadata <- read_delim(file = "public_datasets/Nichterwitz_2016_pmid_27387371/SraRunTable.txt",
                              delim = ",") %>%
  select(run_id = Run,
         sample_id = `Sample Name`,
         species = Organism,
         cell_number = Cell_number,
         source = source_name,
         tissue,
         LibraryLayout,
         AvgSpotLen,
         Instrument) %>%
  mutate(batch = "Nichterwitz",
         species = ifelse(species == "Homo sapiens", "human", "mouse"),
         disease_status = "CTRL",
         source = ifelse(str_detect(source, "Hb9"), 
                         source, 
                         paste0("LCM ", source)),
         files = setNames(
           object = map2_chr(.x = species,
                             .y = run_id,
                             .f = ~ list.files(path = paste0("public_datasets/Nichterwitz_2016_pmid_27387371/",
                                                             ifelse(.x == "human", 
                                                                    "03_salmon_selective_alignment_trim_galore/",
                                                                    "04_salmon_selective_alignment_trim_galore_mouse/"), 
                                                             .y),
                                               pattern = "quant.sf",
                                               full.names = TRUE)), 
           nm = run_id)
  ) %>% 
  select(files, batch, run_id:species, disease_status, cell_number:tissue, LibraryLayout:Instrument) 

write_rds(sample_metadata, "billy/data/input/Nichterwitz_2016_pmid_27387371/sample_metadata.rds")

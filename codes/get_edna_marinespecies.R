#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
######################## Download species list and save ########################

# Load packages ----
library(data.table)
library(dplyr)


# Download species list ----
download.file("https://obis-edna-results.s3.amazonaws.com/output.zip", "species_list.zip",
              method = "wget")

outdir <- "data/species_lists"

fs::dir_create(outdir)

unzip("species_list.zip", exdir = outdir)

file.remove("species_list.zip")


# Load species list and get unique ----
fl <- list.files("data/species_lists/output/", full.names = T, pattern = "Occurrence")

for (f in 1:length(fl)) {
  d <- fread(fl[f])
  
  d_filt <- d %>%
    filter(taxonRank == "species") %>%
    select(scientificName) %>%
    distinct()
  
  if (f == 1) {
    d_tot <- d_filt
  } else {
    d_tot <- bind_rows(d_tot, d_filt) %>%
      distinct()
  }
}

# Load list of marine species
marine_taxa <- readRDS("data/taxa.rds")

marine_taxa <- marine_taxa %>% 
  select(scientificName, scientificName_gbif = input, AphiaID)

d_tot_marine <- d_tot %>%
  left_join(marine_taxa) %>%
  filter(!is.na(AphiaID))

# Save
saveRDS(d_tot_marine, file = "data/marine_species.rds")

### END
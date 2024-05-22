#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
############# Retrieve SST summaries for species (full version) ################

# Load packages ----
library(DBI)
library(furrr)
library(storr)
library(dplyr)
library(tidyr)
library(arrow)
library(terra)
library(polars)
fs::dir_create("results")

# Load files/settings ----
# Database
database <- "~/Research/mpa_europe/mpaeu_shared/occurrence_h3_db.sqlite"
con <- dbConnect(RSQLite::SQLite(), database)

# Settings 
occurrence_table <- "occurrence"

env_vars <- c(
  "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmean_mean",
  "thetao_baseline_depthsurf_mean",
  "thetao_baseline_depthmax_min",
  "thetao_baseline_depthmean_min",
  "thetao_baseline_depthsurf_min",
  "thetao_baseline_depthmax_max",
  "thetao_baseline_depthmean_max",
  "thetao_baseline_depthsurf_max"
)

# Load species list
species_list <- bind_rows(
  lapply(list.files("data/species_list_full", full.names = T), read.csv)
)
species_list <- species_list %>%
  select(species, AphiaID, source_obis, source_gbif, source_dna, group) %>%
  distinct(species, .keep_all = T)

sel_species <- species_list$species
sel_species <- paste0("'", sel_species, "'", collapse = ", ")


# Query database ----
db_query <- glue::glue("
      select species, h3, records
      from {occurrence_table}
      where species in ({sel_species})
                         ")

db_res <- dbSendQuery(con, db_query)
db_species <- dbFetch(db_res)


# Get environmental data ----
# Load environmental layers
env <- rast(paste0("~/Research/mpa_europe/mpaeu_sdm/data/env/current/", env_vars, ".tif"))
names(env) <- env_vars

# Get H3 ids from the environmental layer
biooracle_h3 <- open_dataset("data/biooracle_h3id.parquet")

# Get unique h3 ids
un_h3s <- unique(db_species$h3)

# Get cells for those available
bio_h3_presel <- biooracle_h3 %>%
  filter(h3 %in% un_h3s) %>%
  select(cell, h3) %>%
  collect()

# Get not available
not_available <- un_h3s[!un_h3s %in% unique(bio_h3_presel$h3)]

# Get disks
# Create a temporary folder to hold disks
fs::dir_create("t_disks")

plan(multisession, workers = 8)
not_available_h3_b <- split(not_available, as.integer((seq_along(not_available) - 1) / 2000))

un_h3s_disk <- future_map(not_available_h3_b, function(x){
  
  bio_h3 <- open_dataset("data/biooracle_h3id.parquet")
  
  disk_cells <- lapply(x, function(z){
    dc <- unlist(h3jsr::get_disk(z))
    data.frame(h3_original = z, h3 = dc)
  })
  disk_cells <- bind_rows(disk_cells)
  
  pc <- bio_h3 %>%
    filter(h3 %in% unique(disk_cells$h3)) %>%
    select(cell, h3) %>%
    collect()
  
  disk_cells <- left_join(disk_cells, pc) %>%
    filter(!is.na(cell))
  
  write_parquet(disk_cells,
                paste0("t_disks/", x[1], ".parquet"))
  NULL
}, .progress = T)

# Aggregate expanded disks data
disk_f <- open_dataset("t_disks/")
disk_f <- disk_f %>%
  collect()

# Merge with previous one
bio_h3_presel <- bio_h3_presel %>%
  mutate(h3_original = h3) %>%
  select(h3_original, h3, cell)

bio_h3_sel <- bind_rows(disk_f, bio_h3_presel)


# Extract env data
env_ext <- terra::extract(env, y = as.vector(bio_h3_sel$cell)) %>%
  bind_cols(bio_h3_sel) %>%
  select(-cell, -h3) %>%
  rename(h3 = h3_original) %>%
  group_by(h3) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# Join to species table
all_result <- left_join(db_species, env_ext)


# Get summaries ----
rec_by_sp <- all_result %>%
  group_by(species) %>%
  summarise(records = sum(records))

rec_by_h3 <- all_result %>%
  group_by(h3) %>%
  summarise(records = sum(records))

species_by_h3 <- all_result %>%
  group_by(h3) %>%
  distinct(species, .keep_all = T) %>%
  summarise(n_species = n())

quantile_df <- function(x, probs = c(0, 0.01, 0.05, 0.1, 0.25,
                                     0.5,
                                     0.75, 0.9, 0.95, 0.99, 1)) {
  res <- tibble(metric = c(paste0("q_", probs), "mean", "sd", "no_data", "with_data"), 
                value = c(quantile(x, probs, na.rm = T),
                          mean(x, na.rm = T),
                          sd(x, na.rm = T),
                          sum(is.na(x)),
                          sum(!is.na(x))))
  res
}

final_result <- all_result %>%
  select(-h3, -records) %>%
  pivot_longer(2:length(.), names_to = "variable", values_to = "values") %>%
  group_by(variable, species) %>%
  reframe(quantile_df(values)) %>%
  separate_wider_delim(cols = "variable",
                       names = c("variable", "baseline", "depth", "variant"),
                       delim = "_") %>%
  select(-baseline, -variable)

write_parquet(final_result, "results/species_tsummaries.parquet")
write_parquet(rec_by_sp, "results/records_by_sp.parquet")
write_parquet(rec_by_h3, "results/records_by_h3.parquet")
write_parquet(species_by_h3, "results/species_by_h3.parquet")
write_parquet(species_list, "results/species_list.parquet")

dbDisconnect(con)

rm(list = ls());gc()

### END
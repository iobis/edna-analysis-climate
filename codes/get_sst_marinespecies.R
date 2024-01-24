#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
############# Retrieve SST summaries for species (full version) ################

# Load packages ----
library(DBI)
library(dplyr)
library(tidyr)
library(parquet)
fs::dir_create("results")

# Load files/settings ----
# Database
database <- "~/Downloads/protectedseas/database.sqlite"
con <- dbConnect(RSQLite::SQLite(), database)

# Settings 
gbif_occurrence_table <- "gbif_occurrence"
obis_occurrence_table <- "obis_occurrence"

env_vars <- c(
  "thetao_baseline_depthmax_max", "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_min", "thetao_baseline_depthmean_max",
  "thetao_baseline_depthmean_mean", "thetao_baseline_depthmean_min",
  "thetao_baseline_depthsurf_max", "thetao_baseline_depthsurf_mean",
  "thetao_baseline_depthsurf_min"
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

obis_query <- glue::glue("
      with filtered_data as (
        select species, h3, records
        from {obis_occurrence_table}
        where {obis_occurrence_table}.species in ({sel_species})
        )
      select species, filtered_data.h3, filtered_data.records, {paste(env_vars, collapse = ', ')}
      from filtered_data
      inner join env_current on filtered_data.h3 = env_current.h3;
                         ")

obis_res <- dbSendQuery(con, obis_query)
obis_species <- dbFetch(obis_res)

gbif_query <- glue::glue("
      with filtered_data as (
        select species, h3, records
        from {gbif_occurrence_table}
        where {gbif_occurrence_table}.species in ({sel_species})
        )
      select species, filtered_data.h3, filtered_data.records, {paste(env_vars, collapse = ', ')}
      from filtered_data
      inner join env_current on filtered_data.h3 = env_current.h3;
                         ")

gbif_res <- dbSendQuery(con, gbif_query)
gbif_species <- dbFetch(gbif_res)

all_result <- bind_rows(obis_species, gbif_species)

# Get summaries
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
  res <- tibble(metric = c(paste0("q_", probs), "mean", "sd"), 
                value = c(quantile(x, probs), mean(x), sd(x)))
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
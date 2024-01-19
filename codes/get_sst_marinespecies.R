#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
###################### Retrieve SST summaries for species ######################

# Load packages ----
library(DBI)
library(data.table)
library(dplyr)
library(tidyr)
library(storr)


# Load files ----
species_list <- readRDS("data/marine_species.rds")

database <- "~/Downloads/protectedseas/database.sqlite"

gbif_occurrence_table <- "gbif_occurrence"
obis_occurrence_table <- "obis_occurrence"


# Select environmental variables ----
env_vars <- c(
  "thetao_baseline_depthmax_max", "thetao_baseline_depthmax_mean",
  "thetao_baseline_depthmax_min", "thetao_baseline_depthmean_max",
  "thetao_baseline_depthmean_mean", "thetao_baseline_depthmean_min",
  "thetao_baseline_depthsurf_max", "thetao_baseline_depthsurf_mean",
  "thetao_baseline_depthsurf_min"
)


st <- storr_rds("data/speciestemp_storr")


# Get temperature for each species  ----
job::job(
  {
    con <- dbConnect(RSQLite::SQLite(), database)
    
    for (sp in 1:nrow(species_list)) {
      tsp <- species_list[sp,]
      cli::cat_line("Running species ", cli::col_cyan(tsp$scientificName))
      
      cli::cli_progress_step("Querying OBIS")
      
      sel_species <- tsp$scientificName
      
      obis_query <- glue::glue("
      with filtered_data as (
        select species, h3
        from {obis_occurrence_table}
        where {obis_occurrence_table}.species = '{sel_species}'
        )
      select species, filtered_data.h3, {paste(env_vars, collapse = ', ')}
      from filtered_data
      inner join env_current on filtered_data.h3 = env_current.h3;
                         ")
      
      obis_res <- dbSendQuery(con, obis_query)
      obis_species <- dbFetch(obis_res)
      
      cli::cli_progress_step("Querying GBIF")
      sel_species <- ifelse(!is.na(tsp$scientificName_gbif),
                            ifelse(
                              tsp$scientificName != tsp$scientificName_gbif,
                              tsp$scientificName_gbif,
                              tsp$scientificName
                            ), tsp$scientificName)
      
      gbif_query <- glue::glue("
      with filtered_data as (
        select species, h3
        from {gbif_occurrence_table}
        where {gbif_occurrence_table}.species = '{sel_species}'
        )
      select species, filtered_data.h3, {paste(env_vars, collapse = ', ')}
      from filtered_data
      inner join env_current on filtered_data.h3 = env_current.h3;
                         ")
      
      gbif_res <- dbSendQuery(con, gbif_query)
      gbif_species <- dbFetch(gbif_res)
      
      all_result <- bind_rows(obis_species, gbif_species)
      
      if (nrow(all_result) == 0) {
        all_result <- NULL
      } else {
        q25_fun <- function(x) quantile(x, .25)
        q75_fun <- function(x) quantile(x, .75)
        all_result <- all_result %>%
          select(-h3, -species) %>%
          summarise(across(starts_with("thetao"),
                           list(max = max, min = min,
                                mean = mean, sd = sd, median = median,
                                q25 = q25_fun, q75 = q75_fun),
                           .names = "{.col}${.fn}")) %>%
          pivot_longer(1:length(.), names_to = "variable", values_to = "value") %>%
          separate_wider_delim(cols = "variable", names = c("variable", "variant"),
                               delim = "$") %>%
          pivot_wider(names_from = variant, values_from = value) %>%
          mutate(scientificName = tsp$scientificName,
                 AphiaID = tsp$AphiaID)
        
        
      }
      
      st$set(sp, all_result)
      
      cli::cli_progress_done()
      
    }
    
    dbDisconnect(con)
    job::export("none")
  }
)


# Get ecological information of species
job::job({
  obissdm::mp_get_ecoinfo(species_list = species_list$AphiaID)
  job::export("none")
})

### END
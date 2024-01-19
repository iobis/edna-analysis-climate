#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
######################### Add temperature on database ##########################

# Load packages ----
library(dplyr)
library(future)
library(h3jsr)
library(storr)
library(terra)
library(furrr)
library(ggplot2)
library(DBI)

# Define settings ----
h3_res <- 7

# Get H3 indexes for environmental raster ----

# Environmental layers are downloaded from Bio-ORACLE v3 
# (through the package biooracler)

# Get a base raster
f <- list.files("~/Research/mpa_europe/mpaeu_sdm/data/env/current/",
                full.names = T)
f <- f[!grepl("json", f)]

f_st <- f[grepl("thetao_baseline_depthsurf_mean", f)]

env <- rast(f_st)

# Convert to data frame (cells and xy)
env_df <- as.data.frame(env, xy = T, cells = T)
env_df <- env_df[,1:3]

# Open a storr to hold results
st <- storr_rds("data/envh3_storr")

# Divide in batches to run
env_batches <- split(env_df, as.integer((seq_along(1:nrow(env_df)) - 1) / 300))
ids <- 1:length(env_batches)

rm(env_df)

# Create a function to assign h3 index
to_h3 <- function(obj, idx) {
  h3 <- point_to_cell(obj[,c("x", "y")], res = h3_res)
  h3_df <- cbind(obj, h3 = h3)
  st$set(idx, h3_df)
  return(1)
}

# Map in parallel
plan(multisession, workers = 7)
h3_mapping <- future_map2(env_batches, ids, to_h3, .progress = TRUE) 

# Verify it worked
sum(unlist(h3_mapping))

# Remove unnescessary objects
rm(env_batches, h3_mapping)

# Combine
env_h3 <- st$mget(as.character(1:length(st$list())))

env_h3 <- bind_rows(env_h3)


# Aggregate environmental information ----

# Create a storr to hold results
env_st <- storr_rds("data/envagg_storr/")

# Divide in batches to run
to_split <- unique(env_h3$h3)
to_split <- split(to_split,
                  as.integer((seq_along(to_split) - 1) / 200))
to_split <- lapply(to_split, function(x) data.frame(h3 = x))
names(to_split) <- 1:length(to_split)
to_split <- bind_rows(to_split, .id = "id")

env_h3 <- left_join(env_h3, to_split)

env_batches <- split(env_h3, env_h3$id)
ids <- 1:length(env_batches)

rm(to_split, env_h3)
gc()

#all_rast <- rast(f)

# Create a function to assign h3 index
agg_h3 <- function(obj, idx) {
  
  all_rast <- rast(f)
  
  r_vals <- all_rast[obj$cell]
  
  cols <- gsub("\\.tif", "", basename(f))
  
  names(r_vals) <- cols
  
  agg_vals <- r_vals %>%
    mutate(h3 = obj$h3) %>%
    group_by(h3) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
  
  env_st$set(idx, agg_vals)
  return(1)
}

# Run
plan(multisession, workers = 7)
h3_agg <- future_map2(env_batches, ids, agg_h3, .progress = TRUE)

# Verify it worked
sum(unlist(h3_agg))



# Add to database ----
# Information on how to obtain the database are available on
# https://github.com/iobis/protectedseas-statistics

sqlite_file <- "~/Downloads/protectedseas/database.sqlite"
con <- dbConnect(RSQLite::SQLite(), sqlite_file)

purrr::map(1:length(ids), function(x) {
  env_st$get(as.character(x)) %>%
    dbWriteTable(con, "env_current", ., append = TRUE)
}, .progress = TRUE)

dbSendQuery(con, glue::glue("create index env_current_h3 on env_current(h3)"))

dbDisconnect(con)


# Test query ----
con <- dbConnect(RSQLite::SQLite(), sqlite_file)

site_cells_table <- "site_cells"
gbif_occurrence_table <- "gbif_occurrence"
obis_occurrence_table <- "obis_occurrence"
env_table <- "env_current"
env_vars <- c("thetao_baseline_depthmax_min",
              "thetao_baseline_depthmax_mean",
              "thetao_baseline_depthmax_max")

sel_species <- "Platichthys stellatus"

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

dbDisconnect(con)

head(obis_species)

species_points <- sf::st_as_sf(h3jsr::cell_to_point(obis_species$h3))
species_points <- bind_cols(species_points, obis_species)

base <- rnaturalearth::ne_countries(returnclass = "sf")

ggplot() +
  geom_sf(data = base, fill = "grey90") +
  geom_sf(data = species_points, aes(color = thetao_baseline_depthmax_mean)) +
  scale_color_distiller(name = "thetao mean") +
  coord_sf() +
  theme_minimal()

#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
################### Get H3 indices for environmental layers ####################

# Load packages ----
library(dplyr)
library(future)
library(h3jsr)
library(storr)
library(terra)
library(furrr)

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

# Write as parquet
arrow::write_parquet(env_h3, "data/biooracle_h3id.parquet")

### END
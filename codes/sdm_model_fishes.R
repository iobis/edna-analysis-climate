############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################## Run multiple species distribution models ####################


# Load packages ----
library(obissdm)
library(furrr)
library(progressr)
library(storr)
library(polars)
library(arrow)
library(dplyr)
library(DBI)
library(terra)
source("../../obis/research/marineheritage_sst/functions/sdm_model_function.R")
set.seed(2023)
handlers("cli")
options("progressr.enable" = TRUE)
setwd("~/Research/mpa_europe/mpaeu_sdm")
fs::dir_create("../../obis/research/marineheritage_sst/results/sdms")

# Define settings ----

# General
# The output folder
outfolder <- "../../obis/research/marineheritage_sst/results/sdms"
# An unique code for identifying this model run
outacro <- "mhs"
# Run in parallel? 
run_parallel <- FALSE
# Number of cores for parallel processing
n_cores <- 4

# Modelling
# Algorithms to be used
algos <- c("maxent")
# Should areas be masked by the species depth?
limit_by_depth <- TRUE
# A buffer to be applied to the depth limitation
depth_buffer <- 500
# Assess spatial bias?
assess_bias <- FALSE
# Try to correct spatial bias?
correct_bias <- FALSE

# Create storr to hold results
st <- storr_rds(paste0(outacro, "_storr"))
# If does not exist, add start date
if (!st$exists("startdate")) {
  st$set("startdate", format(Sys.Date(), "%Y%m%d"))
}
# Should the storr object be destructed at the end if the model succeeded?
destroy_storr <- FALSE

# Get list of species
splist <- read_parquet("~/Research/obis/research/marineheritage_sst/results/species_list.parquet")

fishes <- splist[splist$group == "fish", ]

head(fishes)

fishes$where <- apply(fishes[,3:5], 1, function(x){
  if (x[1] | x[2]) {
    "databases"
  } else if (!x[1] & !x[2]) {
    "edna"
  }
})

fishes$sdm_group <- "fishes"



# Fit models ----
# Create a function to save results in a storr object for better control
pmod <- function(sp, gp, sdat, outf, outac, alg, lmd, lmd_buf, assb, corb, p) {
  
  p()
  
  if (st$exists(sp)) {
    if (st$get(as.character(sp))[[1]] %in% c("done", "succeeded")) {
      to_do <- FALSE
    } else {
      to_do <- TRUE
    }
  } else {
    to_do <- TRUE
  }
  
  if (to_do) {
    st$set(sp, "running")
    
    sp_name <- fishes$species[fishes$AphiaID == sp][1]
    
    sel_species <- paste0("'", sp_name, "'")
    
    con <- dbConnect(RSQLite::SQLite(), "~/Research/mpa_europe/mpaeu_shared/occurrence_h3_db.sqlite")
    
    query <- glue::glue('
      select h3, species, records
      from occurrence
      where species in ({sel_species})
    ')
    
    # Get results
    res <- dbSendQuery(con, query)
    records <- dbFetch(res)
    
    dbDisconnect(con)
    
    # Get occurrence coordinates
    rec_coords <- h3jsr::cell_to_point(records$h3)
    rec_coords <- sf::st_coordinates(rec_coords)
    colnames(rec_coords) <- c("decimalLongitude", "decimalLatitude")
    rec_coords <- as.data.frame(rec_coords)
    rec_coords$species <- sp_name
    rec_coords$data_type <- "fit_points"
    
    #print(head(rec_coords))
    
    fit_result <- try(model_species(species = sp,
                                    group = gp,
                                    species_data = rec_coords,
                                    outfolder = outf,
                                    outacro = outac,
                                    algorithms = alg,
                                    limit_by_depth = lmd,
                                    depth_buffer = lmd_buf,
                                    assess_bias = assb,
                                    correct_bias = corb, verbose = T),
                      silent = F)
    
    if (!inherits(fit_result, "try-error")) {
      st$set(sp, fit_result)
    } else {
      st$set(sp, list(status = "failed",
                      error = fit_result))
    }
  }
  
  return(invisible(NULL))
}

# Run models according to the strategy
if (run_parallel) {
  plan(multisession, workers = n_cores)
  
  with_progress({
    p <- progressor(steps = nrow(species_list))
    result <- future_map2(species_list$taxonID, species_list$sdm_group, pmod,
                          outf = outfolder,
                          outac = outacro, alg = algos, lmd = limit_by_depth,
                          lmd_buf = depth_buffer, assb = assess_bias,
                          corb = correct_bias,
                          p = p, .options = furrr_options(seed = T))
  })
} else {
  with_progress({
    p <- progressor(steps = nrow(fishes[1:2,]))
    result <- purrr::map2(fishes$AphiaID[1:2], fishes$sdm_group[1:2], pmod,
                          outf = outfolder,
                          outac = outacro, alg = algos, lmd = limit_by_depth,
                          lmd_buf = depth_buffer, assb = assess_bias,
                          corb = correct_bias,
                          p = p)
  })
}



# Check results ----
# Check if everything was processed
cli::cli_alert_warning("{.val {length(st$list())}} out of {.val {nrow(species_list)}} model{?s} processed.")

# And if so, destroy storr object
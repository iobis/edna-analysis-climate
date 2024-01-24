#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
######################## Sea temperature from Copernicus #######################

# Load packages ----
library(dplyr)
library(stringr)
library(purrr)
library(parallel)
library(terra)
library(glue)
library(arrow)
fs::dir_create("data/agg_temperature")


# Get list of temperature files ----
outfolder <- "data/temperature"

fl <- list.files(outfolder, full.name = T, pattern = "\\.nc")


# Load occurrence info ----
# Bind all occurrence files
occurrence_files <- list.files("data/species_lists/output", "*Occurrence*", full.names = TRUE)

occurrence <- map(occurrence_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
  bind_rows() %>%
  mutate_if(is.character, na_if, "") %>%
  mutate(
    species = ifelse(taxonRank == "species", scientificName, NA),
    aphiaid = as.numeric(str_replace(scientificNameID, "urn:lsid:marinespecies.org:taxname:", ""))
  )

# Get unique sites / get the cell of site ----
base <- rast(res = 1/12)

sites <- occurrence %>%
  group_by(higherGeography, locality, decimalLongitude, decimalLatitude) %>%
  distinct(decimalLongitude, decimalLatitude) %>%
  ungroup() %>%
  filter(!is.na(decimalLongitude)) %>%
  mutate(unsid = 1:nrow(.)) %>%
  rowwise() %>%
  mutate(cell = cellFromXY(base, data.frame(decimalLongitude, decimalLatitude)))
  

sites_unique <- sites %>%
  ungroup() %>%
  distinct(cell, .keep_all = T)


# Load MHS shapefile ----
mhs <- vect("data/shapefiles/marine_world_heritage.gpkg")
# Add an index - we put from 0 to keep consistency with the JupytherHub approaches
mhs$code <- 0:(length(mhs)-1)

# Remove special strings
mhs$name <- stringi::stri_trans_general(str = mhs$name, id = "Latin-ASCII")
sites$higherGeography <- stringi::stri_trans_general(str = sites$higherGeography, id = "Latin-ASCII")
sites_unique$higherGeography <- stringi::stri_trans_general(str = sites_unique$higherGeography, id = "Latin-ASCII")


# Extract info for each site ----
# Define selected depths
depths <- c(0, 25, 50, 75, 100)

for (i in 1:nrow(sites_unique)) {
  
  cli::cli_h1(glue("Running site {i} out of {nrow(sites_unique)}"))
  # Get site code
  site_code <- mhs$code[mhs$name == sites_unique$higherGeography[i]]
  
  tgt <- as.data.frame(sites_unique[i, c("decimalLongitude", "decimalLatitude")])
  
  outfile <- glue("data/agg_temperature/site={site_code[1]}_locality={sites_unique$locality[i]}_unsid={sites_unique$unsid[i]}.parquet")
  
  if (!file.exists(outfile)) {
    
    # Get for the selected depths
    for (d in depths) {
      cli::cli_progress_message(d)
      # Load raster
      # We select the first because for some there is with and without buffer
      # But for this purpose they are the same
      sst <- rast(glue("data/temperature/var=thetao_site={site_code}_depth={d}_product=glorys.nc")[1])
      
      if (d == 0) {
        # Test if is NA
        trast <- sst[[1]]
        tval <- trast[cellFromXY(trast, tgt)]
        
        if (is.na(tval[1,1])) {
          st <- 1
          while (all(is.na(tval[,1])) & st < 5) {
            mat <- c(3, 5, 7, 9)[st]
            adj_cells <- adjacent(trast, cellFromXY(trast, tgt), 
                                  directions = matrix(1, nrow = mat, ncol = mat))
            tval <- trast[as.vector(adj_cells)]
            st <- st+1
          }
          tgc <- as.vector(adj_cells)[which(!is.na(tval[,1]))][1]
        } else {
          tgc <- cellFromXY(trast, tgt)
        }
      }
      
      values <- sst[tgc]
      
      values <- tidyr::pivot_longer(values, 1:length(values),
                                    names_to = "time", values_to = "temperature")
      
      values$time <- gsub(" UTC", "", time(sst))
      values$site <- sites_unique$higherGeography[i]
      values$depth <- d
      
      if (d == 0) {
        all_values <- values
      } else {
        if (nrow(values) > 0 & !is.na(values$temperature[1])) {
          all_values <- bind_rows(all_values, values)
        }
      }
      
      cli::cli_progress_done()
    }
    
    write_parquet(all_values, outfile)
    
    
  } else {
    cli::cat_line("Already done, skipping.")
  }
  
}

write.csv(sites, "data/agg_temperature/all_sites.csv", row.names = F)
write.csv(sites_unique, "data/agg_temperature/unique_sites.csv", row.names = F)

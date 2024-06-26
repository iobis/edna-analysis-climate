#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
######################## Sea temperature from Copernicus #######################

# Load packages ----
library(terra)
library(reticulate)
library(arrow)
library(parallel)
library(dplyr)
library(tidyr)
cm <- import("copernicus_marine_client")


# Settings ----
.user <- rstudioapi::askForPassword("Enter your user")
.pwd <- rstudioapi::askForPassword("Enter your password")

outfolder <- "data/temperature"
fs::dir_create(outfolder)
filename <- "var=thetao"


# Load areas ----
mhs <- vect("data/shapefiles/marine_world_heritage.gpkg")
# Add an index - we put from 0 to keep consistency with the JupytherHub approaches
mhs$code <- 0:(length(mhs)-1)


# Define target dataset, time and depths
dataset <- "cmems_mod_glo_phy_my_0.083deg_P1D-m"
product <- "glorys"
years <- 1992:2021
depths <- c(0, 25, 50, 75, 100, 150, 200, 250, 500, 1000, 2000)
variables <- list("thetao")
lon_lat_buffer <- 0.5 # in degrees
failed <- c() # To see if any failed

for (site_idx in mhs$code) {
  
  sel_site <- mhs[mhs$code == site_idx, ]
  
  long_range <- ext(sel_site)[1:2] + c(-lon_lat_buffer, lon_lat_buffer)
  lat_range <- ext(sel_site)[3:4] + c(-lon_lat_buffer, lon_lat_buffer)
  
  for (depth in depths) {
    
    outfile <- paste0(filename, "_site=", site_idx, "_depth=", depth, "_product=", product, ".nc")
    
    if (!file.exists(paste0(outfolder, "/", outfile))) {
      success <- try(cm$subset(
        dataset_id = dataset,
        variables = variables,
        username = .user,
        password = .pwd,
        minimum_longitude = long_range[1],
        maximum_longitude = long_range[2],
        minimum_latitude = lat_range[1],
        maximum_latitude = lat_range[2],
        start_datetime = paste0(min(years), "-01-01T00:00:00"),
        end_datetime = paste0(max(years), "-12-31T23:59:59"),#"2022-01-31T23:59:59",
        minimum_depth = depth,
        maximum_depth = depth+0.5,
        output_filename = outfile,
        output_directory = outfolder,
        force_download = TRUE
      ), silent = T)
      
      # It will rarely fail, but can happen due to server connection problems.
      # In those cases, sleep and retry
      if (inherits(success, "try-error")) {
        cat("Retrying... \n")
        Sys.sleep(1)
        success <- try(cm$subset(
          dataset_id = dataset,
          variables = variables,
          username = .user,
          password = .pwd,
          minimum_longitude = long_range[1],
          maximum_longitude = long_range[2],
          minimum_latitude = lat_range[1],
          maximum_latitude = lat_range[2],
          start_datetime = paste0(min(years), "-01-01T00:00:00"),
          end_datetime = paste0(max(years), "-12-31T23:59:59"),#"2022-01-31T23:59:59",
          minimum_depth = depth,
          maximum_depth = depth+0.5,
          output_filename = outfile,
          output_directory = outfolder,
          force_download = TRUE
        ), silent = T)
        if (inherits(success, "try-error")) {failed <- c(failed, outfile)}
      }
    } else {
      cat(glue::glue("File for site {site_idx} depth {depth} already exists. Skipping.\n"))
    }
    
  }
  
}

# Use the code below to explore a dataset
# 
# sst_l3s = cm$open_dataset(
#   dataset_id = dataset,
#   #variables = variables,
#   username = .user,
#   password = .pwd,
#   minimum_longitude = -10,
#   maximum_longitude = 10,
#   minimum_latitude = -10,
#   maximum_latitude = 10
# )
# 
# # Print loaded dataset information
# print(sst_l3s)
# 
# 
# # Convert to parquet
# proc <- job::job({
#   lapply(list.files(outfolder, full.name = T, pattern = "\\.nc"), 
#          function(fname, redo = F) {
#            
#            outfile <- gsub("\\.nc", ".parquet", fname)
#            
#            if (file.exists(outfile) & !redo) {
#              cat("File already processed, skipping.\n")
#            } else {
#              cat("Processing", outfile, "\n")
#              r <- rast(fname)
#              
#              if (file.size(fname)/1000/1000 > 1000) {
#                
#                cat("Splitting in pieces (large file) \n")
#                
#                fs::dir_create("temp")
#                
#                grid <- sf::st_make_grid(sf::st_bbox(r), cellsize = 1)
#                
#                mclapply(1:length(grid), function(x){
#                  r_b <- r
#                  r_b <- crop(r, grid[x])
#                  r_b <- as.data.frame(r_b, xy = TRUE, time = TRUE)
#                  r_b <- r_b %>%
#                    pivot_longer(3:length(.), names_to = "time", values_to = "values") %>%
#                    mutate(layer = gsub("\\..*", "", names(r)[1]))
#                  arrow::write_parquet(r_b, sprintf("temp/temp_%06d.parquet", x),
#                                       compression = "gzip")
#                }, mc.cores = 3)
#                
#                gen_files <- open_dataset("temp")
#                
#                gen_files %>%
#                  write_parquet(., outfile, compression = "gzip")
#                
#                fs::dir_delete("temp")
#              } else {
#                
#                r_dat <- as.data.frame(r, xy = TRUE, time = TRUE)
#                
#                r_dat <- r_dat %>%
#                  pivot_longer(3:length(.), names_to = "time", values_to = "values") %>%
#                  mutate(layer = gsub("\\..*", "", names(r)[1]))
#                
#                arrow::write_parquet(r_dat, outfile, compression = "gzip")
#              }
#              
#              
#            }
#            
#            return(invisible(NULL))
#          })
#   
#   lf <- list.files(outfolder, full.name = T, pattern = "\\.parquet")
#   base <- list.files(outfolder, full.name = T, pattern = "\\.nc")[1]
#   
#   zip("temperature_dayly.zip", c(lf, base))
#   
#   job::export("none")
# })




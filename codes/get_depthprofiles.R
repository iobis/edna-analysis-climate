#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
######################## Define depth profiles for sites #######################

# Load packages ----
library(terra)
library(biooracler)
library(reticulate)
cm <- import("copernicus_marine_client")



# Download bathymetry ----
fs::dir_create("data/terrain")
bath <- download_layers("terrain_characteristics", "bathymetry_mean",
                        list(latitude = c(-89.975, 89.975),
                             longitude = c(-179.975, 179.975)), fmt = "raster",
                        directory = "data/terrain",
                        verbose = TRUE)

# Load sites shapefiles ----
mhs <- vect("data/shapefiles/marine_world_heritage.gpkg")

# Get depths ----
max_depth <- rep(NA, length(mhs))
max_depth_within <- rep(NA, length(mhs))
lon_lat_buffer <- 0.5 # in degrees

for (site_idx in 1:length(mhs)) {
  
  sel_site <- mhs[site_idx, ]
  
  site_ext <- ext(sel_site)
  
  site_ext[1:2] <- site_ext[1:2] + c(-lon_lat_buffer, lon_lat_buffer)
  site_ext[3:4] <- site_ext[3:4] + c(-lon_lat_buffer, lon_lat_buffer)
  
  sel_depth <- crop(bath, site_ext)
  
  max_depth[site_idx] <- as.numeric(global(sel_depth, min, na.rm = T))
  
  max_depth_within[site_idx] <- as.numeric(extract(sel_depth, sel_site, fun = "min", na.rm = T)[,2])
}

min(max_depth)
min(max_depth_within, na.rm = T)

max_depth <- data.frame(max_depth = max_depth, site_id = 1:length(mhs))
max_depth_within <- data.frame(max_depth_within = max_depth_within, site_id = 1:length(mhs))

# Set the depths to be used ----
depths <- c(0, 25, 50, 75, 100, 150, 200, 250, 500, 1000, 2000)

# See which depths are being retrieved for each product ----
temp_depths <- sapply(depths, function(x){
  temp_dataset <- cm$open_dataset(
    dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
    username = .user,
    password = .pwd,
    minimum_longitude = -10,
    maximum_longitude = 10,
    minimum_latitude = -10,
    maximum_latitude = 10,
    minimum_depth = x,
    maximum_depth = x+0.5
  )
  temp_dataset$depth$to_dataframe()[,1]
})

o2_depths <- sapply(depths, function(x){
  temp_dataset <- cm$open_dataset(
    dataset_id = "cmems_mod_glo_bgc_my_0.25_P1D-m",
    username = .user,
    password = .pwd,
    minimum_longitude = -10,
    maximum_longitude = 10,
    minimum_latitude = -10,
    maximum_latitude = 10,
    minimum_depth = x,
    maximum_depth = x+0.5
  )
  temp_dataset$depth$to_dataframe()[,1]
})

# Save object
saveRDS(list(
  max_depth_sites = max_depth_within,
  max_depth_sites_buff = max_depth,
  depths_used = depths,
  temperature_depths = temp_depths,
  o2_depths = o2_depths
), file = "data/depth_profile.rds")

#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
##################### Retrieve SST summaries for the sites #####################

# Load packages ----
library(terra)
library(dplyr)
# Folder where environmental files are located, from another project
ffold <- "../../../mpa_europe/mpaeu_sdm/" 

# Load sites list
sites <- read.csv("data/agg_temperature/unique_sites.csv")

# Load a base raster to check which cell is valid
base <- rast(paste0(ffold, "data/env/current/thetao_baseline_depthsurf_mean.tif"))

cells <- terra::extract(base, sites[,c("decimalLongitude", "decimalLatitude")], xy = T, cells = T)

valid_cells <- lapply(1:nrow(cells), function(id){
  if (is.na(cells[id, 2])) {
    
    tval <- cells[id, 2]
    
    starter <- cells$cell[id]
    
    st <- 1
    while (all(is.na(tval)) & st < 5) {
      mat <- c(3, 5, 7, 9)[st]
      adj_cells <- adjacent(base, starter, 
                            directions = matrix(1, nrow = mat, ncol = mat))
      tval <- unname(unlist(base[as.vector(adj_cells)]))
      st <- st+1
    }
    
    tgc <- as.vector(adj_cells)[which(!is.na(tval))][1]
    
    return(tgc)
    
  } else {
    return(cells$cell[id])
  }
})

cells$valid <- unlist(valid_cells)

# Do last check
any(is.na(base[cells$valid]))

cells_xy <- xyFromCell(base, cells$valid)


# Get polygons for the h3 sites
sites_pol <- h3jsr::cell_to_polygon(h3jsr::point_to_cell(cells_xy, res = 7))
sites_pol <- vect(sites_pol)

# Get temperature summaries for each period
all_vars <- list.files(paste0(ffold, "data/env"), recursive = T, full.names = T)
all_vars <- all_vars[!grepl("aux", all_vars)]
all_vars <- all_vars[grepl("thetao", all_vars)]
all_vars <- all_vars[!grepl("range", all_vars)]

# List scenarios and decades
scenarios <- c("current", paste0("ssp", c(126, 245, 370, 460, 585)))
dec <- c("dec50", "dec100")

# Get summaries
summaries <- lapply(scenarios, function(sc) {
  
  layers <- all_vars[grepl(sc, all_vars)]
  
  if (sc == "current") {
    r <- rast(layers)
    names(r) <- gsub("thetao_baseline_", "", gsub(".tif", "", basename(layers)))
    
    r_vals <- terra::extract(r, sites_pol, fun = mean, na.rm = T)
    
    r_vals <- r_vals %>%
      select(-ID) %>%
      mutate(unsid = sites$unsid) %>%
      left_join(sites %>% select(higherGeography, locality, unsid)) %>%
      mutate(period = sc) %>%
      pivot_longer(1:9, names_to = "variable", values_to = "temperature")
    
  } else {
    
    r <- rast(layers[grepl(dec[1], layers)])
    names(r) <- gsub("thetao_", "", gsub(".tif", "", basename(layers[grepl(dec[1], layers)])))
    
    r_vals_a <- terra::extract(r, sites_pol, fun = mean, na.rm = T)
    
    r_vals_a <- r_vals_a %>%
      select(-ID) %>%
      mutate(unsid = sites$unsid) %>%
      left_join(sites %>% select(higherGeography, locality, unsid)) %>%
      mutate(period = paste0(sc, "_", dec[1])) %>%
      pivot_longer(1:9, names_to = "variable", values_to = "temperature")
    
    rm(r)
    
    r <- rast(layers[grepl(dec[2], layers)])
    names(r) <- gsub("thetao_", "", gsub(".tif", "", basename(layers[grepl(dec[2], layers)])))
    
    r_vals_b <- terra::extract(r, sites_pol, fun = mean, na.rm = T)
    
    r_vals_b <- r_vals_b %>%
      select(-ID) %>%
      mutate(unsid = sites$unsid) %>%
      left_join(sites %>% select(higherGeography, locality, unsid)) %>%
      mutate(period = paste0(sc, "_", dec[2])) %>%
      pivot_longer(1:9, names_to = "variable", values_to = "temperature")
    
    r_vals <- bind_rows(r_vals_a, r_vals_b)
    
  }
  
  return(r_vals)
  
})

final_list <- bind_rows(summaries)

final_list <- final_list %>%
  mutate(depth = ifelse(grepl("depthmax", variable), "depthmax",
                        ifelse(grepl("depthsurf", variable), "depthsurf", "depthmean"))) %>%
  mutate(variant = ifelse(grepl("_max", variable), "max",
                      ifelse(grepl("_mean", variable), "mean", "min")))

arrow::write_parquet(final_list, "results/sites_tsummaries.parquet")

### END
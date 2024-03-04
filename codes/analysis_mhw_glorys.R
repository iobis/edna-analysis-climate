#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Mike Burrows, Silas Principe
# Contact: michael.burrows@sams.ac.uk, s.principe@unesco.org
#
######### Analyze Marine Heat Waves on MWHS sites for the eDNA project #########

# Load packages ----
library(terra)
library(heatwaveR)
library(arrow)
library(purrr)
library(dplyr)
library(rjson)
library(future)
library(furrr)

# Settings
outdir <- "results/mhw/"
fs::dir_create(outdir)
set.seed(2024)


# Load data ----
all_sites <- fromJSON(file = "https://raw.githubusercontent.com/iobis/edna-tracker-data/data/generated.json")
all_sites <- all_sites$samples

all_sites_simp <- lapply(all_sites, function(x) {
  if (!x$blank) {
    data.frame(materialSampleID = x$name,
               station = x$station,
               area = x$parent_area_name,
               decimalLongitude = ifelse(is.null(x$area_longitude), NA, x$area_longitude),
               decimalLatitude = ifelse(is.null(x$area_latitude), NA, x$area_latitude))
  } else {
    NULL
  }
})

all_sites_simp <- do.call("rbind", all_sites_simp)

# Load data from used sites
used_sites <- read_parquet("results/sites_tsummaries.parquet")
used_sites <- used_sites$higherGeography
used_sites <- unique(used_sites)

all_sites_filt <- all_sites_simp %>%
  filter(area %in% used_sites) %>%
  distinct(station, .keep_all = T)

# Remove blanks
all_sites_filt <- all_sites_filt[!is.na(all_sites_filt$decimalLongitude),]

# Get sites indexes
# Data for the sites was downloaded based on the shapefile
# We need to traceback the IDs to the sites names
mhs <- vect("data/shapefiles/marine_world_heritage.gpkg")
# Add an index - we put from 0 to keep consistency with the JupytherHub approache
mhs$code <- 0:(length(mhs)-1)

mhs_info <- data.frame(area = mhs$name, code = mhs$code, buffer = mhs$buffer)


# Run analysis for each site in parallel ----
plan(multisession, workers = 4)

retrieve_mhw <- function(id) {
  
  depth_val <- 0 # To try other depths change here
  
  mhs_sel <- mhs_info[mhs_info$area == id,]
  
  if (nrow(mhs_sel) > 0) {
    
    if (nrow(mhs_sel) > 1) {
      
      # Load rasters
      to_merge <- lapply(mhs_sel$code, function(x){
        r <- rast(paste0("data/temperature/var=thetao_site=", x, "_depth=", depth_val, "_product=glorys.nc"))
      })
      
      tvar <- merge(sprc(to_merge))
      
    } else {
      tvar <- rast(paste0("data/temperature/var=thetao_site=", mhs_sel$code[1], "_depth=", depth_val, "_product=glorys.nc"))
    }
    
    dates <- terra::time(tvar)
    dates <- as.Date(dates)
    
    valid_cells <- as.data.frame(tvar[[1]], cell = T, xy = T)[,1:3]
    
    if (nrow(valid_cells) > 2000) {
      sa <- sample(1:nrow(valid_cells), 2000)
      valid_cells <- valid_cells[sa,]
    }
    
    for (i in 1:nrow(valid_cells)) {
      wdf <- unname(unlist(tvar[valid_cells$cell[i]]))
      wdf <- data.frame(t = dates, temp = wdf)
      ts <- ts2clm(wdf, climatologyPeriod = c("1993-01-01", "2012-12-31"))
      mhw <- detect_event(ts)
      mhw <- mhw$event
      mhw$rast_cell <- valid_cells$cell[i]
      mhw$decimalLongitude <- valid_cells$x[i]
      mhw$decimalLatitude <- valid_cells$y[i]
      if (i == 1) {
        mhw_result <- mhw
      } else {
        mhw_result <- bind_rows(mhw_result, mhw)
      }
    }
    
    write_parquet(mhw_result, paste0(outdir, "mhw_area=",
                                     gsub("\\s", "-", gsub("[[:punct:]]", "", tolower(id))),
                                     "_depth=", depth_val, ".parquet"))
  }
  
  return(invisible(NULL))
  
}

future_map(unique(all_sites_filt$area), retrieve_mhw, .progress = T, .options = furrr_options(seed = T))

### END
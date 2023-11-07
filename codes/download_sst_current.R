############# eDNA Expeditions - UNESCO World Heritage Marine Sites ############
################# Data analysis of the project - Thermal niche #################
# August of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################ Retrieve SST satellite data for current period ################

# Load packages and define settings ----
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(arrow)
# Path to save files
save_path <- "data/sst/current/"


# Load Python functions ----
reticulate::source_python("functions/inspect_cop.py")
# reticulate::py_help(inspect_cop) # To see help run this line
reticulate::source_python("functions/retrieve_sst_fun.py")
# reticulate::py_help(retrieve_cop) # To see help run this line


# Get credentials for Copernicus Marine service ----
cop_user <- rstudioapi::askForPassword("Copernicus username")
cop_pass <- rstudioapi::askForPassword("Copernicus password")

# Load marine heritage sites shapefile ----
mhs_path <- "data/shapefiles/WorldMarineHeritageSites_v2.shp"
mhs <- st_read(mhs_path)

# Limit to those that are confirmed on the website
# (ask for a csv list)
mhs <- mhs[grepl(paste0(
  c(
    "wadden", "shark", "noronha", "rocas", "french austral", "lord howe",
    "sundarbans", "coiba", "revillagigedo", "caledonia", "calanche",
    "sanganeb", "everglades", "aldabra", "belize", "tubbataha",
    "simangaliso", "banc", "ningaloo", "socotra", "Península Valdés"
  ), collapse = "|"
), mhs$FULL_NAME, ignore.case = T),]

mhs_sites <- unique(mhs$MRGID)

# Oceans (optional) #
oceans <- mregions::mr_shp(key = "MarineRegions:goas", read = TRUE, maxFeatures = 1000)
oceans <- st_as_sf(oceans)

inter <- st_intersects(mhs, st_make_valid(oceans), sparse = F)

mhs$OCEAN <- oceans$name[apply(inter, 1, which.max)]


# Open dataset ----
# Dataset name
dataset <- "METOFFICE-GLO-SST-L4-REP-OBS-SST"
# Variable of interest
variable <- "analysed_sst"
# Time range
time_window <- 1992:2021

# Inspect and open dataset
dataset_info <- inspect_cop(dataset, cop_user, cop_pass, variable = variable, plot = T)


# Get mean, maximum, minimum and sd for each site and save ----
# Because all will have the same configurations, we first create a function
get_metric <- function(site, metric, time_window) {
  retrieve_cop(c(dataset, dataset_info), # Supplying dataset name and dataset object
               variable, # Variable to open
               cop_user, cop_pass, # Copernicus access information
               shape = mhs_path, # Shapefile path
               shape_var = "MRGID", # Shapefile variable to filter
               shape_filter = site, # Filter value
               time_range = time_window, # Time range
               ret_type = "df", # Returning type
               metric = metric, # Metric of summary
               group_by_lonlat = TRUE, # Group by grouper and summarise by whole raster
               res_by_month = TRUE, # Aggregate by month
               grouper = "time", # Get for each time step
               mask_by_shape = FALSE, # Mask by the shapefile
               k_to_celsius = TRUE, # Convert from degree to celsius
               plot = F) # Plot during execution
}

# We run in loop, that way if there is any problem we can try again
sites_metrics <- list()

for (z in 16:length(mhs_sites)) {
  # We load the dataset again at each start to avoid the server disconnecting
  dataset_info <- inspect_cop(dataset, cop_user, cop_pass, variable = variable, plot = F)
  for (k in 1:length(time_window)) {
    tr <- paste0(time_window[k], c("-01-01", "-12-31"))
    temp_data <- get_metric(mhs_sites[z],
                            metric = c("mean", "std", "max", "min", "median"),
                            time_window = tr)
    if (k == 1) {
      sites_metrics[[z]] <- temp_data
    } else {
      sites_metrics[[z]] <- rbind(sites_metrics[[z]],
                                  temp_data)
    }
  }
}

# 
# sites_metrics <- lapply(mhs_sites, get_metric,
#                         metric = c("mean", "std", "max", "min", "median"))

# Convert to a single data frame containing sites info
sites_metrics <- lapply(sites_metrics, function(x){x$time <- rownames(x);x})
names(sites_metrics) <- mhs_sites

sites_metrics <- bind_rows(sites_metrics, .id = "MRGID")

sites_info <- st_drop_geometry(mhs[,c("FULL_NAME", "COUNTRY", "MRGID")])
sites_info <- distinct(sites_info, MRGID, .keep_all = T)

sites_metrics <- left_join(sites_metrics, sites_info, by = "MRGID")


# Save results ----
write_parquet(sites_metrics, paste0(save_path, "mhs_sst_current.parquet"))


# Plot results ----
sites_metrics$time <- lubridate::as_date(sites_metrics$time)

sites_metrics <- sites_metrics %>%
  group_by(MRGID) %>%
  mutate(general_mean = mean(mean)) %>%
  mutate(detrend = mean - general_mean) %>%
  select(-general_mean)

ggplot(sites_metrics, aes(x = time, y = mean))+
  geom_line() +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .3)+
  theme_light() +
  facet_wrap(~ MRGID, scales = "free_y") 

sites_metrics$state <- ifelse(sites_metrics$detrend > 0, "Higher", "Lower")

ggplot(sites_metrics, aes(x = time, y = detrend))+
  geom_hline(yintercept = 0) +
  #geom_line(aes(color = state)) +
  geom_area(aes(x=time, y=ifelse(detrend<0, detrend, 0)), fill="#1093C8") +
  geom_area(aes(x=time, y=ifelse(detrend>0, detrend, 0)), fill="#C72B10") +
  theme_bw() +
  facet_wrap(~ MRGID, scales = "free_y")







####
# Animated grap
library(dygraphs)
library(xts)

wide_metrics <- sites_metrics %>%
  select(time, MRGID, mean, sd) %>%
  filter(MRGID %in% c(26836, 64215)) %>%
  mutate(upr = mean+sd, lwr = mean-sd) %>%
  select(-sd) %>%
  tidyr::pivot_wider(names_from = MRGID, values_from = c(mean, upr, lwr))

smetric <- xts(x = wide_metrics, order.by = wide_metrics$time)

# Finally the plot
(p <- dygraph(smetric) %>%
    dySeries(c("lwr_26836", "mean_26836", "upr_26836"), label = "26836") %>%
    dySeries(c("lwr_64215", "mean_64215", "upr_64215"), label = "64215") %>%
    dyOptions(labelsUTC = TRUE, fillGraph=FALSE, fillAlpha=0.1, drawGrid = FALSE,
              colors = RColorBrewer::brewer.pal(3, "Set2")) %>%
    dyAxis("y", label = "Temperature (°C)") %>%
    dyRangeSelector() %>%
    dyCrosshair(direction = "vertical") %>%
    dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE,
                highlightSeriesOpts = list(strokeWidth = 2))  #%>%
  #dyRoller(rollPeriod = 1)
)


smetric_a <- xts(x = wide_metrics[,c("time", "mean_26836", "lwr_26836", "upr_26836")], order.by = wide_metrics$time)
smetric_b <- xts(x = wide_metrics[,c("time", "mean_64215", "lwr_64215", "upr_64215")], order.by = wide_metrics$time)


(p1 <- dygraph(smetric_a, group = "sst") %>%
    dySeries(c("lwr_26836", "mean_26836", "upr_26836"), label = "26836") %>%
    dyOptions(labelsUTC = TRUE, fillGraph=FALSE, fillAlpha=0.1, drawGrid = FALSE,
              colors = RColorBrewer::brewer.pal(3, "Set2")) %>%
    dyAxis("y", label = "Temperature (°C)") %>%
    dyRangeSelector() %>%
    dyCrosshair(direction = "vertical") %>%
    dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE,
                highlightSeriesOpts = list(strokeWidth = 2))  #%>%
  #dyRoller(rollPeriod = 1)
)

(p2 <- dygraph(smetric_b, group = "sst") %>%
    dySeries(c("lwr_64215", "mean_64215", "upr_64215"), label = "64215") %>%
    dyOptions(labelsUTC = TRUE, fillGraph=FALSE, fillAlpha=0.1, drawGrid = FALSE,
              colors = RColorBrewer::brewer.pal(3, "Set2")) %>%
    dyAxis("y", label = "Temperature (°C)") %>%
    dyRangeSelector() %>%
    dyCrosshair(direction = "vertical") %>%
    dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE,
                highlightSeriesOpts = list(strokeWidth = 2))  #%>%
  #dyRoller(rollPeriod = 1)
)
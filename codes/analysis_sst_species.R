#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data analysis #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
############################# Species thermal limits ###########################

# Load packages ----
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


# Load files ----
# Temperature per species dataset
temp_sp <- open_dataset("results/species_tsummaries.parquet")

temp_sp_filt <- temp_sp %>%
  filter(depth == "depthsurf") %>%
  filter(variant == "mean") %>%
  collect()

# Temperature on sites
temp_sites <- read_parquet("results/sites_tsummaries.parquet")

# Get sites list
sites_f <- list.files("data/species_list_full", full.names = T)

sites_list <- lapply(sites_f, read.csv)
names(sites_list) <- gsub("_", " ", gsub("\\.csv", "", basename(sites_f)))

sites_list <- bind_rows(
  sites_list, .id = "higherGeography"
)

sites_list_distinct <- sites_list %>%
  select(species, AphiaID, source_obis, source_gbif, source_dna, higherGeography) %>%
  group_by(higherGeography) %>%
  distinct(species, .keep_all = T)

# Get number of records 
species_rec <- read_parquet("results/records_by_sp.parquet")


# Bind information ----
# Bind species temperature
temp_sp_wide <- temp_sp_filt %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  select(-depth, -variant) %>%
  relocate(no_data, .after = with_data)

sites_list_binded <- sites_list_distinct %>%
  left_join(temp_sp_wide) %>%
  filter(!is.na(q_0.25))

# Bind temperature on sites
temp_sites_filt <- temp_sites %>%
  filter(depth == "depthsurf") %>%
  filter(variant == "mean") %>%
  group_by(higherGeography, period) %>%
  summarise(temperature = mean(temperature)) %>%
  pivot_wider(names_from = period, values_from = temperature) %>%
  mutate(higherGeography = tolower(higherGeography)) %>%
  mutate(higherGeography = gsub("[']", " ", higherGeography)) %>%
  mutate(higherGeography = gsub("[:,]", "", higherGeography))

colnames(temp_sites_filt)[2:length(temp_sites_filt)] <- paste0("site_", colnames(temp_sites_filt)[2:length(temp_sites_filt)])

sites_list_binded <- sites_list_binded %>%
  left_join(temp_sites_filt) %>%
  filter(!is.na(site_current))

# Bind number of records 
sites_list_binded <- sites_list_binded %>%
  left_join(species_rec)

# Create tag for from where it comes
sites_list_binded <- sites_list_binded %>%
  mutate(source_obis = ifelse(as.numeric(source_obis) == 1, 1, 0)) %>%
  mutate(source_gbif = ifelse(as.numeric(source_gbif) == 1, 2, 0)) %>%
  mutate(source_dna = ifelse(as.numeric(source_dna) == 1, 4, 0)) %>%
  mutate(where = source_obis + source_gbif + source_dna) %>%
  mutate(where = case_when(
    where == 4 ~ "eDNA",
    where %in% c(3, 2, 1) ~ "OBIS/GBIF",
    where %in% c(6, 5, 7) ~ "Both"
  ))


# Save file
write_parquet(sites_list_binded, "results/tsummaries_aggregated.parquet")


# Get thermal ranges ----
species_thermal_range <- sites_list_binded %>%
  filter(records >= 10) %>%
  mutate(loc_range = (site_current - q_0.01)/(q_0.99 - q_0.01))

# Plot to verify
ggplot(species_thermal_range) +
  geom_boxplot(aes(x = loc_range))


# Create a function to plot for just one site
plot_site <- function(site) {
  
  if (is.numeric(site)) {
    site <- unique(sites_list_binded$higherGeography)[site]
  }
  
  sel_site <- sites_list_binded %>%
    filter(higherGeography == site) %>%
    mutate(where = ifelse(source_obis & source_gbif & !source_dna, "OBIS/GBIF", 
                          ifelse(source_dna & !source_obis & !source_gbif, "eDNA", "Both"))) 
  
  plot_a <- sel_site %>%
    mutate(status = ifelse(q_0.5 > min(sel_site$site_current), "Above", "Below")) %>%
    ggplot() +
    geom_jitter(aes(y = higherGeography, x = q_0.5, color = status), alpha = .2, height = 0.2) +
    geom_vline(xintercept = min(sel_site$site_current), linetype = 2) +
    theme_minimal() +
    ylab(NULL) + xlab (NULL) + ggtitle("All species", site) +
    theme(legend.title = element_blank(),
          axis.text.y = element_blank())
  
  plot_b <- sel_site %>%
    mutate(status = ifelse(q_0.5 > min(sel_site$site_current), "Above", "Below")) %>%
    ggplot() +
    geom_jitter(aes(y = where, x = q_0.5, color = status), alpha = .2, height = 0.2) +
    geom_vline(xintercept = min(sel_site$site_current), linetype = 2) +
    theme_minimal() +
    ylab(NULL) + ggtitle("By type") +
    theme(legend.title = element_blank())
  
  plot_a / plot_b + plot_layout(guides = "collect", heights = c(1,2)) & theme(legend.position = "bottom")
  
}

plot_site(3)

### END
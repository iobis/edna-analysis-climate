#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas Principe, Mike Burrows
# Contact: s.principe@unesco.org; michael.burrows@sams.ac.uk
#
############## Plot eDNA data from UNESCO World Heritage Sites #################

# Load packages ----
library(dplyr)
library(terra)
library(tidyr)
library(stringr)
library(purrr)
library(arrow)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggdist)
library(see)
fs::dir_create("figures")

# Load data ----
# Load temperature summaries by site
speciesthermsite <- read_parquet("data/species_sites_risk.parquet")
colnames(speciesthermsite)[1] <- "species"

# Process the data ----
# Convert to data.table
speciesthermsitedt <- data.table(speciesthermsite)

# Load species list
specieslist <- open_dataset("data/output/occurrence.parquet")
groups <- read.csv("data/groups.csv")
group_lookup <- setNames(groups$group, paste(groups$rank, groups$taxon, sep = "_"))

# Add a new column for group
specieslist <- specieslist %>%
    filter(taxonRank == "species") %>%
    collect() %>%
    distinct(scientificName, .keep_all = T) %>%
    mutate(group = coalesce(
        group_lookup[paste("class", class, sep = "_")],
        group_lookup[paste("order", order, sep = "_")],
        group_lookup[paste("phylum", phylum, sep = "_")]
    )) %>%
    select(species = scientificName, group)


# Merge group information
speciesthermsitedt_gp <- merge(speciesthermsitedt, specieslist[,c("species", "group")], by = "species")

# Filter for fish only
fishthermsitedt <- speciesthermsitedt_gp[group %in% c("fish", "sharks"), , ]
fishthermsitedt <- fishthermsitedt[!is.na(fishthermsitedt$baseline_depthsurf_mean),]

# Get functional traits data
funcdiv <- rfishbase::species(unique(fishthermsitedt$species))
funcdiv$to_include <- NA
funcdiv$to_include[funcdiv$DepthRangeDeep <= 200] <- TRUE

table(funcdiv$DemersPelag[is.na(funcdiv$to_include)])

funcdiv$to_include[is.na(funcdiv$to_include) & 
    funcdiv$DemersPelag %in% c("reef-associated", "pelagic-neritic")] <- TRUE

table(funcdiv$DemersPelag[is.na(funcdiv$to_include)])

sel_species <- funcdiv %>%
    filter(to_include) %>%
    select(species = Species, DepthRangeDeep, habitat = DemersPelag)

# Convert to longer format
fishthermsitedt_long <- fishthermsitedt %>%
  select(species, higherGeography, baseline_depthsurf_mean, contains("depthsurf_dec50"), limit_dsurf) %>%
  filter(species %in% unique(sel_species$species)) %>% # Remove this line to not filter by the traits
  pivot_longer(cols = contains("depthsurf"), names_to = "scenario",
               values_to = "sst") %>%
  mutate(scenario = case_when(
    scenario == "baseline_depthsurf_mean" ~ "site_current",
    scenario == "ssp126_depthsurf_dec50_mean" ~ "site_ssp126_dec50",
    scenario == "ssp245_depthsurf_dec50_mean" ~ "site_ssp245_dec50",
    scenario == "ssp370_depthsurf_dec50_mean" ~ "site_ssp370_dec50",
    scenario == "ssp585_depthsurf_dec50_mean" ~ "site_ssp585_dec50",
  ))


# Create a function to see at which scenario the species become at risk
scenario_order <- c("site_current", "site_ssp126_dec50", 
                    "site_ssp585_dec50")
get_class <- function(x) {
  if(all(x < 0)) {
    return("not_risk")
  } else {
    return(scenario_order[min(which(x >= 0))])
  }
}

# Retrieve sites names and latitudes
sites <- vect("data/shapefiles/marine_world_heritage.gpkg")
sites$latitude <- unlist(sapply(1:length(sites), function(x){
  unname(terra::geom(centroids(sites[x,]))[,4])
}))
sites_org <- data.frame(sites = sites$name, latitude = sites$latitude)

# Change names of sites
sites_org <- sites_org %>%
  rename(higherGeography = sites) %>%
  mutate(higherGeography = case_when(
    higherGeography == "Wadden Sea" ~ "Wadden Sea",
    higherGeography == "Archipiélago de Revillagigedo" ~ "Revillagigedo",
    higherGeography == "iSimangaliso Wetland Park" ~ "iSimangaliso",
    higherGeography == "Ningaloo Coast" ~ "Ningaloo coast",
    higherGeography == "Socotra Archipelago" ~ "Socotra Archipelago",
    higherGeography == "Belize Barrier Reef Reserve System" ~ "Belize Barrier Reef",
    higherGeography == "Cocos Island National Park" ~ "Cocos Island",
    higherGeography == "Lord Howe Island Group" ~ "Lord Howe Island",            
    higherGeography == "Sundarbans National Park" ~ "Sundarbans",                    
    higherGeography == "Península Valdés" ~ "Peninsula Valdès",
    higherGeography == "Gulf of Porto: Calanche of Piana, Gulf of Girolata, Scandola Reserve" ~ "Scandola",         
    higherGeography == "Coiba National Park and its Special Zone of Marine Protection" ~ "Coiba",             
    higherGeography == "Lagoons of New Caledonia: Reef Diversity and Associated Ecosystems" ~ "Lagoons of New Caledonia",         
    higherGeography == "Shark Bay, Western Australia" ~ "Shark Bay",                                               
    higherGeography == "Tubbataha Reefs Natural Park" ~ "Tubbataha Reefs",                                              
    higherGeography == "Brazilian Atlantic Islands: Fernando de Noronha and Atol das Rocas Reserves" ~ "Fernando de Noronha",
    higherGeography == "Everglades National Park" ~ "Everglades",                                                  
    higherGeography == "Aldabra Atoll" ~ "Aldabra Atoll",                                                             
    higherGeography == "Banc d'Arguin National Park" ~ "Banc d'Arguin",                                               
    higherGeography == "French Austral Lands and Seas" ~ "French Austral Lands",
  )) %>%
  filter(!is.na(higherGeography)) %>%
  arrange(desc(latitude)) %>%
  distinct(higherGeography, .keep_all = T)

# Update names on the longer object
fishthermsitedt_long <- fishthermsitedt_long %>%
  mutate(higherGeography = case_when(
    higherGeography == "Wadden Sea" ~ "Wadden Sea",
    higherGeography == "Archipiélago de Revillagigedo" ~ "Revillagigedo",
    higherGeography == "iSimangaliso Wetland Park" ~ "iSimangaliso",
    higherGeography == "Ningaloo Coast" ~ "Ningaloo coast",
    higherGeography == "Socotra Archipelago" ~ "Socotra Archipelago",
    higherGeography == "Belize Barrier Reef Reserve System" ~ "Belize Barrier Reef",
    higherGeography == "Cocos Island National Park" ~ "Cocos Island",
    higherGeography == "Lord Howe Island Group" ~ "Lord Howe Island",            
    higherGeography == "The Sundarbans" ~ "Sundarbans",                    
    higherGeography == "Peninsula Valdès" ~ "Peninsula Valdès",
    higherGeography == "Gulf of Porto: Calanche of Piana, Gulf of Girolata, Scandola Reserve" ~ "Scandola",         
    higherGeography == "Coiba National Park and its Special Zone of Marine Protection" ~ "Coiba",             
    higherGeography == "Lagoons of New Caledonia: Reef Diversity and Associated Ecosystems" ~ "Lagoons of New Caledonia",         
    higherGeography == "Shark Bay, Western Australia" ~ "Shark Bay",                                               
    higherGeography == "Tubbataha Reefs Natural Park" ~ "Tubbataha Reefs",                                              
    higherGeography == "Brazilian Atlantic Islands: Fernando de Noronha and Atol das Rocas Reserves" ~ "Fernando de Noronha",
    higherGeography == "Everglades National Park" ~ "Everglades",                                                  
    higherGeography == "Aldabra Atoll" ~ "Aldabra Atoll",                                                             
    higherGeography == "Banc d'Arguin National Park" ~ "Banc d'Arguin",                                               
    higherGeography == "French Austral Lands and Seas" ~ "French Austral Lands",
  ))

# Add means
species_lims_full <- read_parquet(file.path("data", "species_thermal_lims.parquet"))
species_lims_full <- species_lims_full %>%
    filter(metric == "q_0.95") %>%
    filter(depth == "depthsurf") %>%
    select(species, q_0.95 = value)

fishthermsitedt_long <- left_join(fishthermsitedt_long, species_lims_full)

# Get the number of species in each site
sites_species <- fishthermsitedt_long %>%
  group_by(higherGeography, species) %>%
  summarise(n = n()) %>%
  group_by(higherGeography) %>%
  summarise(n = n())

# Join with the sites object
sites_org <- left_join(sites_org, sites_species)

# Create the plot object with all information
plot_obj <- fishthermsitedt_long %>%
  filter(scenario %in% c("site_current", "site_ssp126_dec50", "site_ssp585_dec50")) %>%
  select(species, higherGeography, scenario, differences = sst, q_0.95) %>%
  filter(!is.na(differences)) %>%
  group_by(species, higherGeography) %>%
  mutate(point_class = get_class(differences)) %>%
  mutate(higherGeography = factor(higherGeography, levels = rev(sites_org$higherGeography),
                                  labels = rev(paste(sites_org$higherGeography, paste0("(", sites_org$n, ")")))))
plot_obj$point_class <- factor(plot_obj$point_class, levels = c("not_risk", "site_current", "site_ssp126_dec50", "site_ssp585_dec50"))

# Extract sites SST
sites_sst <- read_parquet("data/sites_sst_all.parquet")

sites_sst <- sites_sst %>%
    mutate(higherGeography = case_when(
    higherGeography == "Wadden Sea" ~ "Wadden Sea",
    higherGeography == "Archipiélago de Revillagigedo" ~ "Revillagigedo",
    higherGeography == "iSimangaliso Wetland Park" ~ "iSimangaliso",
    higherGeography == "Ningaloo Coast" ~ "Ningaloo coast",
    higherGeography == "Socotra Archipelago" ~ "Socotra Archipelago",
    higherGeography == "Belize Barrier Reef Reserve System" ~ "Belize Barrier Reef",
    higherGeography == "Cocos Island National Park" ~ "Cocos Island",
    higherGeography == "Lord Howe Island Group" ~ "Lord Howe Island",            
    higherGeography == "The Sundarbans" ~ "Sundarbans",                    
    higherGeography == "Peninsula Valdès" ~ "Peninsula Valdès",
    higherGeography == "Gulf of Porto: Calanche of Piana, Gulf of Girolata, Scandola Reserve" ~ "Scandola",         
    higherGeography == "Coiba National Park and its Special Zone of Marine Protection" ~ "Coiba",             
    higherGeography == "Lagoons of New Caledonia: Reef Diversity and Associated Ecosystems" ~ "Lagoons of New Caledonia",         
    higherGeography == "Shark Bay, Western Australia" ~ "Shark Bay",                                               
    higherGeography == "Tubbataha Reefs Natural Park" ~ "Tubbataha Reefs",                                              
    higherGeography == "Brazilian Atlantic Islands: Fernando de Noronha and Atol das Rocas Reserves" ~ "Fernando de Noronha",
    higherGeography == "Everglades National Park" ~ "Everglades",                                                  
    higherGeography == "Aldabra Atoll" ~ "Aldabra Atoll",                                                             
    higherGeography == "Banc d'Arguin National Park" ~ "Banc d'Arguin",                                               
    higherGeography == "French Austral Lands and Seas" ~ "French Austral Lands",
  )) %>%
  filter(!is.na(higherGeography))

sites_sst <- sites_sst %>%
  select(higherGeography, baseline_depthsurf_mean, contains("depthsurf_dec50")) %>%
  pivot_longer(cols = contains("depthsurf"), names_to = "scenario",
               values_to = "sst") %>%
  mutate(scenario = case_when(
    scenario == "baseline_depthsurf_mean" ~ "site_current",
    scenario == "ssp126_depthsurf_dec50_mean" ~ "site_ssp126_dec50",
    scenario == "ssp245_depthsurf_dec50_mean" ~ "site_ssp245_dec50",
    scenario == "ssp370_depthsurf_dec50_mean" ~ "site_ssp370_dec50",
    scenario == "ssp585_depthsurf_dec50_mean" ~ "site_ssp585_dec50",
  )) %>%
  filter(scenario %in% c("site_current", "site_ssp126_dec50", "site_ssp585_dec50")) %>%
  filter(higherGeography %in% unique(fishthermsitedt_long$higherGeography)) %>%
  mutate(point_class = scenario) %>%
  left_join(sites_org[,c("higherGeography", "n")]) %>%
  mutate(higherGeography = paste(higherGeography, paste0("(", n, ")")))
sites_sst$point_class <- factor(sites_sst$point_class, levels = c("not_risk", "site_current", "site_ssp126_dec50", "site_ssp585_dec50"))


# Make plots ----
# Left plot - dots
plot_a <- plot_obj %>%
  group_by(higherGeography) %>%
  distinct(species, .keep_all = T) %>%
  ggplot() +
  geom_violinhalf(aes(x = higherGeography, y = q_0.95),
                  position = position_nudge(x = 0.3)) +
  geom_jitter(aes(x = higherGeography, y = q_0.95, color = point_class), alpha = .3,
              width = 0.12, shape = 16, show.legend = F) +
  geom_linerange(data = sites_sst, aes(ymin = sst-0.1, ymax = sst+0.1, x = higherGeography, color = point_class),
                 linewidth = 7, show.legend = T) +
  scale_color_manual(values = c("grey90", "#d4b9da",  "#df65b0",  "#980043"),
                     labels = c("Not at risk", "Current", "SSP1", "SSP5")) +
  coord_flip() +
  ylab("Temperature (°C)") + xlab(NULL) +
  theme_light() +
  scale_y_continuous(breaks = c(0, 10, 20, 30), limits = c(0, 32)) +
  theme(panel.border = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0), "cm"));plot_a

# Right plot - bars
plot_b <- plot_obj %>%
  group_by(higherGeography, scenario) %>%
  summarise(percentage = sum(differences > 0)/length(differences)) %>%
  ggplot() +
  geom_bar(aes(x = higherGeography, y = percentage, fill = scenario), stat = "identity", position="dodge", width = .6) +
  scale_fill_manual(values = c("#d4b9da",  "#df65b0",  "#980043")) +
  scale_y_continuous(labels = c("0%", "", "50%", "", "100%")) +
  coord_flip() +
  ylab("Potentially Affected Fraction") + xlab(NULL) +
  theme_light() +
  theme(panel.border = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        plot.margin = unit(c(0,0,0,0), "cm"));plot_b

# Put together
plot_a + plot_b + plot_layout(widths = c(4, 1))

# Save
ggsave("figures/paf_fishes_filtered_v2_nl.png", width = 12, height = 7)
ggsave("figures/paf_fishes_filtered_v2_nl.svg", width = 12, height = 7)

plot_c <- plot_obj %>%
  group_by(higherGeography, scenario) %>%
  summarise(percentage = sum(differences > 0)/length(differences)) %>%
  ggplot() +
  geom_bar(aes(x = higherGeography, y = percentage, fill = scenario), stat = "identity", position="dodge", width = .6) +
  scale_fill_manual(values = c("#d4b9da",  "#df65b0",  "#980043")) +
  scale_y_continuous(labels = c("0%", "", "50%", "", "100%")) +
  coord_flip() +
  ylab("Potentially Affected Fraction") + xlab(NULL) +
  theme_light() +
  theme(panel.border = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        plot.margin = unit(c(0,0,0,0), "cm"));plot_c

ggsave("figures/paf_fishes_filtered_onlybars_v2_nl.png", width = 8, height = 7)
ggsave("figures/paf_fishes_filtered_onlybars_v2_nl.svg", width = 8, height = 7)

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
library(fishualize)
library(ggdist)
library(see)
fs::dir_create("figures")

# Load data ----
# Load functional traits data
funcdiv <- read.csv("data/funcdiversity.csv", sep = "\t")
table(funcdiv$Habitat)
summary(funcdiv$DepthRangeDeep)

# Select shallow water species
summary(funcdiv$DepthRangeDeep) # There are 1348 species with NA values

# Limit to those with maximum depth 200m
sel_species <- funcdiv %>%
  filter(!is.na(DepthRangeDeep)) %>%
  filter(DepthRangeDeep <= 200)



# Load temperature summaries by site
speciesthermsite <- read_parquet("results/tsummaries_aggregated.parquet")

# Process the data ----
# Convert to data.table
speciesthermsitedt <- data.table(speciesthermsite)

# Change `where` column IDs
speciesthermsitedt$where <- ifelse(
  speciesthermsitedt$where == "Both" | speciesthermsitedt$where == "OBIS/GBIF",
  "Databases",
  "eDNA"
)

# Load species list
specieslist <- read_parquet("results/species_list.parquet")

# Merge group information
speciesthermsitedt_gp <- merge(speciesthermsitedt, specieslist[,c("species", "group")], by = "species")

# Filter for fish only
fishthermsitedt <- speciesthermsitedt_gp[group == "fish", , ]

# Convert to longer format
fishthermsitedt_long <- fishthermsitedt %>%
  filter(species %in% unique(sel_species$X)) %>% # Remove this line to not filter by the traits
  pivot_longer(cols = starts_with("site"), names_to = "scenario",
               values_to = "sst") %>%
  filter(scenario %in% c("site_current", "site_ssp126_dec50", "site_ssp245_dec50",
                         "site_ssp370_dec50", "site_ssp585_dec50"))

# Create a function to see at which scenario the species become at risk
scenario_order <- c("site_current", "site_ssp126_dec50", 
                    "site_ssp585_dec50")
get_class <- function(x) {
  if(all(x < 0)) {
    return("not_risk")
  } else {
    return(scenario_order[min(which(x > 0))])
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
    higherGeography == "Península Valdés" ~ "Peninsula Valdes",
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
    higherGeography == "wadden sea" ~ "Wadden Sea",
    higherGeography == "archipielago de revillagigedo" ~ "Revillagigedo",
    higherGeography == "isimangaliso wetland park" ~ "iSimangaliso",
    higherGeography == "ningaloo coast" ~ "Ningaloo coast",
    higherGeography == "socotra archipelago" ~ "Socotra Archipelago",
    higherGeography == "belize barrier reef reserve system" ~ "Belize Barrier Reef",
    higherGeography == "cocos island national park" ~ "Cocos Island",
    higherGeography == "lord howe island group" ~ "Lord Howe Island",            
    higherGeography == "the sundarbans" ~ "Sundarbans",                    
    higherGeography == "peninsula valdes" ~ "Peninsula Valdes",
    higherGeography == "gulf of porto calanche of piana gulf of girolata scandola reserve" ~ "Scandola",         
    higherGeography == "coiba national park and its special zone of marine protection" ~ "Coiba",             
    higherGeography == "lagoons of new caledonia reef diversity and associated ecosystems" ~ "Lagoons of New Caledonia",         
    higherGeography == "shark bay western australia" ~ "Shark Bay",                                               
    higherGeography == "tubbataha reefs natural park" ~ "Tubbataha Reefs",                                              
    higherGeography == "brazilian atlantic islands fernando de noronha and atol das rocas reserves" ~ "Fernando de Noronha",
    higherGeography == "everglades national park" ~ "Everglades",                                                  
    higherGeography == "aldabra atoll" ~ "Aldabra Atoll",                                                             
    higherGeography == "banc d arguin national park" ~ "Banc d'Arguin",                                               
    higherGeography == "french austral lands and seas" ~ "French Austral Lands",
  ))

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
  select(species, higherGeography, q_0.95, scenario, sst) %>%
  mutate(differences = sst - q_0.95) %>%
  group_by(higherGeography, species) %>%
  mutate(q_0.95 = mean(q_0.95), point_class = get_class(differences)) %>%
  mutate(higherGeography = factor(higherGeography, levels = rev(sites_org$higherGeography),
                                  labels = rev(paste(sites_org$higherGeography, paste0("(", sites_org$n, ")")))))
plot_obj$point_class <- factor(plot_obj$point_class, levels = c("not_risk", "site_current", "site_ssp126_dec50", "site_ssp585_dec50"))

# Extract sites SST
sites_sst <- fishthermsitedt_long %>%
  filter(scenario %in% c("site_current", "site_ssp126_dec50", "site_ssp585_dec50")) %>%
  group_by(higherGeography, scenario) %>%
  summarise(sst = mean(sst)) %>%
  mutate(point_class = scenario) %>%
  left_join(sites_org[,c("higherGeography", "n")]) %>%
  mutate(higherGeography = paste(higherGeography, paste0("(", n, ")")))
sites_sst$point_class <- factor(sites_sst$point_class, levels = c("not_risk", "site_current", "site_ssp126_dec50", "site_ssp585_dec50"))


# Make plots ----
# Left plot - dots
plot_a <- plot_obj %>%
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
ggsave("figures/paf_fishes_filtered_v1.png", width = 12, height = 7)
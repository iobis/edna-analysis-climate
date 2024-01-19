#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
######################## Plot SST summaries for species ########################

# Load packages ----
library(storr)
library(tidyverse)


# Load storr
st <- storr_rds("data/speciestemp_storr")

# Select species
sp <- "Gnatholepis anjerensis"

dat <- st$get(st$list()[2])

# Select mode of life/depth
depth <- "depthmean"

# Plot
sel_dat <- dat %>%
  filter(grepl(depth, dat$variable)) %>%
  separate_wider_delim(cols = variable, delim = "_",
                       names = c("variable", "base", "depth", "variant")) %>%
  select(-variable, -base, -depth)

ggplot(sel_dat) +
  geom_linerange(aes(y = variant, xmin = q25, xmax = q75, color = variant),
                 linewidth = 2) +
  geom_pointrange(aes(y = variant, x = median, xmin = min, xmax = max, color = variant),
                  size = 1) +
  geom_vline(xintercept = min(sel_dat$min), linetype = 2,
             color = "#0B4F71") +
  geom_text(aes(y = "min", x = min(min), label = paste("Lower limit: ", round(min(min), 2))),
            nudge_y = 0.3, nudge_x = 4,
            color = "#0B4F71") +
  geom_vline(xintercept = max(sel_dat$max), linetype = 2,
             color = "#860D35") +
  geom_text(aes(y = "max", x = max(max), label = paste("Upper limit: ", round(max(max), 2))),
            nudge_y = 0.3, nudge_x = -4,
            color = "#860D35") +
  scale_color_manual(values = c("#390099", "#9e0059", "#ff0054")) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank()) +
  xlab("Temperature") + ylab(NULL)
  

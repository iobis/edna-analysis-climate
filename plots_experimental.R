#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas Principe, Mike Burrows
# Contact: s.principe@unesco.org, michael.burrows@sams.ac.uk
#
########## Read and analyse eDNA data from UNESCO World Heritage Sites #########

# Load packages ----
library(dplyr)
library(stringr)
library(purrr)
library(arrow)
library(data.table)
library(ggplot2)
fs::dir_create("figures")
set.seed(2023)

# Load files and convert to data.table ----
# Species summaries
speciestherm <- read_parquet("results/species_tsummaries.parquet")
speciesthermdt <- data.table(speciestherm)

# Site summaries
speciesthermsite <- read_parquet("results/tsummaries_aggregated.parquet")
speciesthermsitedt <- data.table(speciesthermsite)

# Species list
specieslist <- read_parquet("results/species_list.parquet")


# Select relevant depth/variant ----
speciesthermsst <- speciesthermdt[depth == "depthsurf" & variant == "mean", , ]
speciesthermsstc <- dcast(speciesthermdt,
                          species ~ metric,
                          value.var = "value",
                          fun = mean)
  
# Get community thermal index
ctihighergeog <- speciesthermsitedt[, list(
  ctiavg = mean(q_0.5),
  sdcti = sd(q_0.5),
  str = mean(q_0.9 - q_0.1),
  sstavg = mean(site_current),
  nspp = .N
), by = c("higherGeography", "where")]

for (i in 1:100) {
  
  resampled <- speciesthermsitedt[, list(
    q_0.5 = sample(q_0.5, .N, replace = T),
    q_0.9 = sample(q_0.9, .N, replace = T),
    q_0.1 = sample(q_0.1, .N, replace = T),
    site_current = sample(site_current, .N, replace = T)
    #resample = 1:.N#sample(1:.N, .N, replace = T)
  ), by = c("higherGeography", "where")]
  
  ctihighergeog_temp <- resampled[, list(
    ctiavg = mean(q_0.5),
    sdcti = sd(q_0.5),
    str = mean(q_0.9 - q_0.1),
    sstavg = mean(site_current),
    nspp = .N
  ), by = c("higherGeography", "where")]
  
  if (i == 1) {
    ctihighergeog_boot <- ctihighergeog_temp
  } else {
    ctihighergeog_boot <- rbind(ctihighergeog_boot, ctihighergeog_temp)
  }
  
}

# Make plots ----
plot(data=ctihighergeog,ctiavg~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,
     xlab="Sea surface temperature (°C)", ylab="Community Temperature Index (°C)",main="Fish")

ggplot(ctihighergeog) +
  geom_point(aes(x = sstavg, y = ctiavg, color = where), shape = 16, size = 2) +
  xlab("Sea surface temperature (°C)") + 
  ylab("Community Temperature Index (°C)") +
  theme_light() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())


### Bootstrap version
ctihighergeog$higherGeography <- substr(ctihighergeog$higherGeography, 1, 10)
ctihighergeog$higherGeography <- factor(ctihighergeog$higherGeography,
                                        levels = unique(ctihighergeog$higherGeography[order(ctihighergeog$sstavg)]))

ctihighergeog_boot$higherGeography <- substr(ctihighergeog_boot$higherGeography, 1, 10)
ctihighergeog_boot$higherGeography <- factor(ctihighergeog_boot$higherGeography,
                                             levels = levels(ctihighergeog$higherGeography))

ctihighergeog_boot_sum <- ctihighergeog_boot[, list(
  ctiavg_low = quantile(ctiavg, 0.25),
  ctiavg_median = median(ctiavg),
  ctiavg_high = quantile(ctiavg, 0.75),
  sdcti_low = quantile(sdcti, 0.25),
  sdcti_median = median(sdcti),
  sdcti_high = quantile(sdcti, 0.75),
  str_low = quantile(str, 0.25),
  str_median = median(str),
  str_high = quantile(str, 0.75)
), by = c("higherGeography", "where")]

ctihighergeog_sst <- ctihighergeog_boot[,list(
  sstavg = mean(sstavg)
), by = c("higherGeography", "where")]
ctihighergeog_sst$what = "Average site SST"

ggplot() +
  geom_point(data = ctihighergeog_sst, aes(x = higherGeography, y = sstavg, fill = what),
             shape = "|", size = 6, alpha = .1) +
  geom_jitter(data = ctihighergeog_boot, aes(x = higherGeography, y = ctiavg, color = where), 
             shape = 16, size = 1, alpha = .05,
             position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.height = 0, jitter.width = 0.1)) +
  geom_linerange(data = ctihighergeog_boot_sum,
                 aes(x = higherGeography, ymin = ctiavg_low, ymax = ctiavg_high,
                     color = where), linewidth = 1,
                 position = position_dodge(width = 0.5)) +
  geom_point(data = ctihighergeog, aes(x = higherGeography, y = ctiavg, color = where), 
             shape = 16, size = 2.5, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("#00c3c9", "#5151db", "#f58000")) +
  scale_fill_manual(values = c("grey70")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  ylab("Temperature (°C)") + 
  xlab(NULL) +
  theme_light() +
  coord_flip() +
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3, linetype = 2),
        legend.title = element_blank())

ggsave("figures/test_newcti.png", width = 7, height = 10)

ctihighergeog_sst_b <- ctihighergeog_sst[,list(
  sstavg = mean(sstavg)
), by = c("higherGeography")]
ctihighergeog_sst_b$what <- ctihighergeog_sst$what[1]

ggplot() +
  geom_col(data = ctihighergeog_sst_b, aes(x = higherGeography, y = sstavg, fill = what),
             alpha = .3, width = 0.8) +
  geom_jitter(data = ctihighergeog_boot, aes(x = higherGeography, y = ctiavg, color = where), 
              shape = 16, size = 1, alpha = .05,
              position = position_jitterdodge(dodge.width = 0.5,
                                              jitter.height = 0, jitter.width = 0.1)) +
  geom_linerange(data = ctihighergeog_boot_sum,
                 aes(x = higherGeography, ymin = ctiavg_low, ymax = ctiavg_high,
                     color = where), linewidth = 1,
                 position = position_dodge(width = 0.5)) +
  geom_point(data = ctihighergeog, aes(x = higherGeography, y = ctiavg, color = where), 
             shape = 16, size = 2.5, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("#00c3c9", "#5151db", "#f58000")) +
  scale_fill_manual(values = c("grey70")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(expand = c(0,0,0,2)) +
  ylab("Temperature (°C)") + 
  xlab(NULL) +
  theme_light() +
  coord_flip() +
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3, linetype = 2),
        legend.title = element_blank())

ggsave("figures/test_newcti_opt2.png", width = 7, height = 10)

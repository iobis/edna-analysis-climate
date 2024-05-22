# Plot SST maps with differences for glossy report

# Load packages
library(terra)
library(ggplot2)
library(patchwork)

# Set directory
sh <- "../../../mpa_europe/mpaeu_sdm/data/env/"

# Load environmental layers
current <- rast(paste0(sh, "current/thetao_baseline_depthsurf_mean.tif"))
ssp1 <- rast(paste0(sh, "future/ssp126/thetao_ssp126_depthsurf_dec50_mean.tif"))
ssp5 <- rast(paste0(sh, "future/ssp585/thetao_ssp585_depthsurf_dec50_mean.tif"))

# Aggregate for faster plotting
current <- aggregate(current, 2)
ssp1 <- aggregate(ssp1, 2)
ssp5 <- aggregate(ssp5, 2)

# Get deltas
ssp1_d <- ssp1 - current
ssp5_d <- ssp5 - current

# Convert to data.frame
current_df <- as.data.frame(current, xy = T)
ssp1_d_df <- as.data.frame(ssp1_d, xy = T)
ssp5_d_df <- as.data.frame(ssp5_d, xy = T)

# Change column names
colnames(ssp1_d_df)[3] <- colnames(ssp5_d_df)[3] <- "delta"
colnames(current_df)[3] <- "val"

# Verify range for breaks
range(ssp1_d_df$delta)
range(ssp5_d_df$delta)
# See how many are above 2.5, as difference is small
length(ssp5_d_df$delta[ssp5_d_df$delta > 2.5])

# Add extremes on one small cell to overcome ggplot limits
ssp1_d_df[1,3] <- -1.5
ssp1_d_df[2,3] <- 2.5
ssp5_d_df[1,3] <- -1.5
ssp5_d_df[2,3] <- 2.5

# Change high values of SSP5 (i.e 2.53) to be within the scale
ssp5_d_df$delta[ssp5_d_df$delta > 2.5] <- 2.5

# Load world shapefile
world <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")

# Make SST plot
a <- ggplot() +
  geom_sf(data = world, color = "white", fill = "white") +
  geom_contour_filled(data = current_df, aes(x = x, y = y, z = val)#,
                      #breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)#,
                      #color = "grey80", linewidth = 0.2
  ) +
  scale_fill_brewer(palette = "BuPu", direction = 1,
                    guide = guide_coloursteps(show.limits = T),
                    name = "SST (°C)") +
  coord_sf(expand = F, label_axes = list(bottom = "", right = "")) +
  ggtitle("Current period") +
  theme_light() +
  theme(plot.background = element_blank(),
        #panel.border = element_blank(),
        axis.title = element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.1, "npc"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.justification = "center",
        legend.title = element_text(hjust = 0.5),
        panel.grid = element_blank());a

# Make SSP1 plot
b <- ggplot() +
  geom_sf(data = world, color = "white", fill = "white") +
  geom_contour_filled(data = ssp1_d_df, aes(x = x, y = y, z = delta),
                      breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5)#,
                      #color = "grey80", linewidth = 0.2
  ) +
  scale_fill_manual(values = rev(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3")),
                    guide = guide_coloursteps(show.limits = T),
                    name = "Difference in temperature (°C)") +
  coord_sf(expand = F,
           label_axes = list(bottom = "", right = "")) +
  ggtitle("SSP1 (2050)") +
  theme_light() +
  theme(plot.background = element_blank(),
        #panel.border = element_blank(),
        axis.title = element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.1, "npc"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.justification = "center",
        legend.title = element_text(hjust = 0.5),
        panel.grid = element_blank());b

# Make SSP5 plot  
c <- ggplot() +
  geom_sf(data = world, color = "white", fill = "white") +
  geom_contour_filled(data = ssp5_d_df, aes(x = x, y = y, z = delta),
                      breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5)#,
                      #color = "grey80", linewidth = 0.2
  ) +
  scale_fill_manual(values = rev(c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3")),
                    guide = guide_coloursteps(show.limits = T),
                    name = "Difference in temperature (°C)") +
  coord_sf(expand = F,
           label_axes = list(bottom = "", right = "")) +
  ggtitle("SSP5 (2050)") +
  theme_light() +
  theme(plot.background = element_blank(),
        #panel.border = element_blank(),
        axis.title = element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.1, "npc"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.justification = "center",
        legend.title = element_text(hjust = 0.5),
        panel.grid = element_blank());c

# Compose
pc <- ((plot_spacer() + a + plot_layout(guides = "keep")) | (b + c + plot_layout(guides = "collect"))) & theme(legend.position = "bottom",
                                                                                                               legend.key.height = unit(0.025, "npc"),
                                                                                                               legend.key.width = unit(0.04, "npc"),
                                                                                                               legend.title = element_text(hjust = 0.5, size = 12),
                                                                                                               legend.text = element_text(hjust = 0.5, size = 12))

# Save
ggsave(filename = "figures/sst_isomaps.png", plot = pc, width = 28, height = 8)

# END
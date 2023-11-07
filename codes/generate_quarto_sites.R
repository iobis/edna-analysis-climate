#### Create quarto page for each site
library(tidyverse)

mhs <- sf::read_sf("../data/shapefiles/WorldMarineHeritageSites_v2.shp")

mhs <- mhs[grepl(paste0(
  c(
    "wadden", "shark", "noronha", "rocas", "french austral", "lord howe",
    "sundarbans", "coiba", "revillagigedo", "caledonia", "calanche",
    "sanganeb", "everglades", "aldabra", "belize", "tubbataha",
    "simangaliso", "banc", "ningaloo", "socotra", "Península Valdés"
  ), collapse = "|"
), mhs$FULL_NAME, ignore.case = T),]

sites_info <- sf::st_drop_geometry(mhs[,c("FULL_NAME", "COUNTRY", "MRGID")])
sites_info <- distinct(sites_info, MRGID, .keep_all = T)

temperature <- arrow::read_parquet("../data/sst/current/mhs_sst_current.parquet")

temperature <- temperature %>%
  mutate(time = as_date(time)) %>%
  filter(year(time) == 2021) %>%
  group_by(MRGID) %>%
  summarise(mean = mean(mean))

for (z in 1:nrow(sites_info)) {
  rl <- readLines("sites/_sitemodel.qmd")
  rl <- gsub("<SITE_NAME>", sites_info$FULL_NAME[z], rl)
  rl <- gsub("<SITE_CODE>", sites_info$MRGID[z], rl)
  rl <- gsub("<SITE_MEAN_SST>",
             round(temperature$mean[temperature$MRGID == sites_info$MRGID[z]], 1),
             rl)
  writeLines(rl, paste0("sites/site", sites_info$MRGID[z], ".qmd"))
}
# eDNA expeditions
## Environmental information on Marine Heritage Sites

This repository documents the process to obtain and analyse environmental information (sea temperature and oxygen) on Marine Heritage Sites. Environmental information is downloaded from Copernicus, from the following data sources:

- SST: [Global Ocean Physics Reanalysis](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description)
- Oxygen: [Global Ocean Biogeochemistry Hindcast](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description)

You can see the live page at https://iobis.github.io/marineheritage-sst

## Codes

### Data download

There are several ways of obtaining Copernicus satellite data. We provide code for two possible pathways:

1. From WEkEO (`download_temperature_wekeo.ipynb`) - WEkEO is a service from Copernicus that provide a virtual environment (JupyterHub) for satellite data processing. Because all Copernicus data is acessible directly from the virtual environment, this is the fastest way of obtaining the data. WEkEO is free and an account can be created here. Once you have set up your virtual environment, open the Jupyter notebook provided here to download the data.
2. Using the new Copernicus API (`download_*_toolbox.R`) - Copernicus introduced major changes in its Marine Data store in 2023, including a new [toolbox](https://help.marine.copernicus.eu/en/articles/7949409-copernicus-marine-toolbox-introduction) for data access. Unfortunately, the solution is based only on Python and thus there is no R equivalent for it. Using from R relies on system interface to the CLI or link through the `reticulate` package (approach used here). The code `get_depth_profiles.R` is used to obtain the nearest available depth from the chosen depths.

To be in line with the future changes in the Copernicus services, we adopted the pathway 2 to obtain the data. It is necessary to have **Python** installed and the toolbox ([instructions here](https://help.marine.copernicus.eu/en/articles/7970514-copernicus-marine-toolbox-installation), we recommend using `pip`).

Oxygen data download proceeded the same way as temperature (use the code with the word "oxygen").

### Data processing

__Marine heat waves__

Data analysis for Marine Heat Waves was done using the package [`heatwaveR`](https://robwschlegel.github.io/heatwaveR/). All analysis is on the code `analysis_mhw.R` (an alternative version for GLORYS product is on the code `analysis_mhw.R`).

__Temperature by species__

Codes are in the order they are needed.

- `get_biooracle_h3.R`: obtain H3 index for Bio-ORACLE layers.
- `get_sst_marinespecies.R`: aggregate species list, get unique species, retrieve records from database and extract temperature from Bio-ORACLE layers based on the H3 index. After that, get summaries and save results.
- `get_sst_sites.R`: extract temperature information for each site (collection site, not _higherGeography_).
- `analysis_occurrence_data.R`: calculate CTI (Community Thermal Index) and make plots by higher geography.
- `analysis_sst_species.R`: visualize thermal limits by species/sites.

__Others__

- `get_depthprofiles.R`: see which is the actual depth that was used when obtaining the data from Copernicus
- `generate_quarto_sites.R`: generate quarto pages by site
- `plot_sst_marinespecies.R`: plot temperature information for a selected species

## Data analysis

After retrieving temperature information for each species occurrence from OBIS/GBIF (aggregated by H3 cells on a resolution of 7), we generate a table with the different quantiles for each species. This is done still on the `get_sst_marinespecies.R` script. The same data is obtained for each site.

We join both tables to analyse the thermal limits for each species and to assess how species in each site might be at risk by climate change.

----

<img src="https://obis.org/images/logo.png" width="200">

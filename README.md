# eDNA expeditions
## Environmental information on Marine Heritage Sites

This repository documents the process to obtain and analyse environmental information (sea temperature and oxygen) on Marine Heritage Sites. Environmental information is downloaded from Copernicus, from the following data sources:

- SST: [Global Ocean Physics Reanalysis](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description)
- Oxygen: [Global Ocean Biogeochemistry Hindcast](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/description)

You can see the live page at https://iobis.github.io/marineheritage_sst

## Codes

### Data download

There are several ways of obtaining Copernicus satellite data. We provide code for three possible pathways:

1. From WEkEO (`download_temperature_wekeo.ipynb`) - WEkEO is a service from Copernicus that provide a virtual environment (JupyterHub) for satellite data processing. Because all Copernicus data is acessible directly from the virtual environment, this is the fastest way of obtaining the data. WEkEO is free and an account can be created here. Once you have set up your virtual environment, open the Jupyter notebook provided here to download the data.
2. Using OpenDAP (`download_temperature_opendap.R`) - OpenDAP is a data access protocol that is widely used to access satellite data. OpenDAP is also a very quick way to acess Copernicus data and here we adapt code provided by Jorge Assis to subset and download the information we need. Note, however, that the OpenDAP support for Copernicus is being deprecated in mid 2024 in favor of the new Python API (see below).
3. Using the new Copernicus API (`download_*_toolbox.R`) - Copernicus introduced major changes in its Marine Data store in 2023, including a new [toolbox](https://help.marine.copernicus.eu/en/articles/7949409-copernicus-marine-toolbox-introduction) for data access. Unfortunately, the solution is based only on Python and thus there is no R equivalent for it. Using from R relies on system interface to the CLI or link through the `reticulate` package (approach used here). The code `get_depth_profiles.R` is used to obtain the nearest available depth from the chosen depths.

To be in line with the future changes in the Copernicus services, we adopted the pathway 3 to obtain the data. It is necessary to have **Python** installed and the toolbox ([instructions here](https://help.marine.copernicus.eu/en/articles/7970514-copernicus-marine-toolbox-installation), we recommend using `pip`).

Oxygen data download proceeded the same way as temperature (use the code with the word "oxygen").

### Data processing

Data analysis for Marine Heat Waves and Cold Spells was done using the package [`heatwaveR`](https://robwschlegel.github.io/heatwaveR/). More details to be added after the workshop.

----

<img src="https://obis.org/images/logo.png" width="200">

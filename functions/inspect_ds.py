############# eDNA Expeditions - UNESCO World Heritage Marine Sites ############
################# Data analysis of the project - Thermal niche #################
# August of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################# Inspect satellite dataset from other sources #################
def inspect_ds(dataset, username, password, plot = False, variable = None):
  
  """Inspect datasets stored in Zarr/Xarray format
  
  This function enables to inspect satellite data prior to retrieving
  the full dataset or part of that. For Copernicus datasets, use inspect_cop
  instead.
  
  Parameters
  ----------
  dataset : str
      Link to the dataset (e.g. https://mur-sst.s3.us-west-2.amazonaws.com/zarr-v1)
  plot : bool
      Should the function plot the maps during the process? Default is False.
  plot_bbox : str, optional
      If ploting, a bounding box to the plot. This is specially relevant for
      datasets that are very large (i.e. high-resolution). Should be list on the
      format [xmin, xmax, ymin, ymax].
  variable : str, optional
      Name of the variable (e.g. analyzed_sst). Should be supplied if plot = True.
      
  Returns
  -------
  xarray.Dataset
    This object can be supplied to the function retrieve_cop
      
  Typical usage example:
  -------
    inspect_ds(dataset = 'https://mur-sst.s3.us-west-2.amazonaws.com/zarr-v1',
               plot = True, plot_bbox = [-90, -30, -42.5, 42.5],
               variable = 'analysed_sst')
                  
  Depends on the following packages:
  -------                
  matplotlib, xarray, zarr
  You may need to install lxml, aiohttp, requests
  """
  
  import matplotlib.pyplot as plt
  import xarray as xr
  
  if plot and variable == None:
    raise ValueError('When plot=True, variable should be supplied. Check arguments.')
  
  print('Loading dataset...', flush = True) 
  DS = xr.open_dataset(dataset)
  
  print(DS, flush = True)
  
  if plot:
    sliced = DS[variable].isel(time=1)
    sliced = sliced.sel(
            lon=slice(plot_bbox[0], plot_bbox[1]),
            lat=slice(plot_bbox[2], plot_bbox[3]))
    if 'depth' in DS.dims:
      sliced = sliced.isel(depth=1)
    plt.figure()
    sliced.plot()
    plt.show()

  return DS

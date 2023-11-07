############# eDNA Expeditions - UNESCO World Heritage Marine Sites ############
################# Data analysis of the project - Thermal niche #################
# August of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############## Retrieve satellite data from Marine Copernicus (3D) #############
def retrieve_cop_3d(dataset, variable, username, password, shape, shape_var, shape_filter, time_range, depth_range,
                ret_type = 'raster', agg_depth = False, metric = None, group_by_lonlat = False, res_by_month = False,
                grouper = 'time', save_path = None, mask_by_shape = True, k_to_celsius = False, plot = False):
  
  """Retrieve information from Copernicus Marine data products with 3 dimensions
  
  Retrieve, select and optionally summarise Copernicus satellite data with a depth 
  component. This function also enable to retrieve data from other datasets opened
  with the function 'inspect_ds'. In that case, you should supply the result of 
  the function as the dataset argument, following the instructions on the 
  parameters section. Note: the 'depth' component should be named 'depth' in the
  dataset.
  
  Parameters
  ----------
  dataset : str or list
      Name of the dataset (e.g. dataset-armor-3d-rep-monthly) or a list containing
      the name of the dataset and the Xarray dataset object. The second approach is
      relevant when you will loop over several areas or make distinct metric operations.
      For getting the Xarray dataset, use first the function inspect_cop (or alternatively
      inspect_ds for non-Copernicus datasets). NOTE: the list SHOULD be supplied 
      in that order (i.e. name, dataset).
  variable : str
      Name of the variable (e.g. to)
  username : str
      Copernicus username
  password : str
      Copernicus password
  shape : str
      Path to shapefile containing regions to mask
  shape_var : str
      Attribute (column) of the shapefile used to filter (e.g. MRGID)
  shape_filter : str or numeric
      Code defining the filter on the selected attribute
  time_range : list, str
      A list with two values defining the time range. Should be on the format
      YYYY-MM-DD. E.g. ['2010-01-15', '2010-01-17']
  depth_range : list, str
      A list with two values defining the depth range. Should be on the format
      min-max. E.g. [0, 10]. Note that in some cases the depth may be stored 
      as negative, in which case it should be [-10, 0]. Inspect your dataset
      to know which one to use.
  ret_type : str
      Which value to return. There are 2 types available:
        'raster' - a netcdf raster is saved in save_path containing the whole
        sliced (and optionally masked) dataset.
        'df' - a data frame like object is returned. Note that depending on the
        size of your dataset this may be problematic! Use with care.
  agg_depth : bool
      If True, then first the dataset is aggregated over depth (i.e. the mean
      over the depths for each time step). Default is False.
  metric : str or list, optional
      Metric to summarise the results. The metric is usually applied over the layers,
      but can also produce a single value by changing the group_by_lonlat argument. Note
      that if you set group_by_lonlat = True, then you should also define ret_type to
      'df', otherwise the function will stop.
      Metric can be mean, min, max, median, sum, and std. There is one special available
      which calculate both the mean and standard deviation: 'meanstd' If you supply
      a list of metrics, then all will be calculated and merged in a single file or
      data frame. The saved file (if one is saved) will be named "{name}_summarised.nc"
      instead of "{name}_{metric}.nc".
  group_by_lonlat : bool
      If set to True and a metric is chosen, then the summary is applied over
      the lon-lat grouped by time (i.e. a single value per time step). The time
      step may be grouped differently by changing the 'grouper' argument. In any case,
      if the dataset was not aggregated by depth, metrics are calculated for each depth.
  res_by_month : bool
      If set to True, daily information is resampled to monthly values (mean of 
      each month). This option is only relevant if you have daily data. Note that
      after resampling, metrics are calculated with the resampled version.
  grouper : str
      The grouper before calculating the metrics. It is usually 'time', but can 
      also be 'time.month' to aggregate by month or other valid Xarray value. 
      If group is done by lon-lat (group_by_lonlat=True) and the grouper is 
      different than 'time', then the data is first aggregated by the grouper 
      (e.g. get the month means over all years) before calculating the lon-lat 
      metric per time step.
  save_path : str, optional
      Path to save the netcdf if ret_type = 'raster'. Should be on the format 
      "path/folder/", i.e. with the final "/"
  mask_by_shape : bool
      If True, then the values of the Xarray are masked by the shape, otherwise
      it is used only to crop the Xarray to the bounds of the shape.
  k_to_celsius : bool
      If True, then the values are converted from Kelvin to Celsius (C = k - 273.15).
      Default is False, so the users need to really check if the dataset is in Kelvin.
  plot : bool
      Should the function plot the maps during the process? Default is False.
      
  Returns
  -------
  saved netcdf
      files are saved on path, or:
  dataframe
      containing the values
      
  Typical usage example:
  -------
    user = 'myuser'
    pass = 'mypassword'
    shp_path = 'path/to/my_shape.shp'
    out_path = 'data/saved/'
    retrieve_cop('dataset-armor-3d-rep-monthly', 'to',
                  user, pass, shp_path, 'MRGID', '26836', ['2010-01-15', '2010-01-17'],
                  [0, 10])
                  
  Depends on the following packages:
  -------                
  matplotlib, xarray, numpy, pandas, geopandas and regionmask
  For saving netcdf, also netCDF4
  You may need to install lxml
  
  Possible combinations of parameters:
  -------
  Depending on the way you combine the group_by_lonlat, res_by_month and grouper
  parameters you can generate a large array of results. 
  
  If you have a daily dataset, spanning several years, you can get monthly maps
  of the chosen metric by first using res_by_month = True and then grouper = 'time.month'.
  If instead you want the maximum over all months use grouper = 'time'.
  
  For getting a metric over the area instead, you set group_by_lonlat = True. That
  way, the metric is calculated over the pixels (lon-lat coordinates). So, the maximum
  will return the pixel with maximum value. Because the way this calculation is done,
  grouping by different time steps necessarily involve first aggregating the results by
  the refered time step. So using grouper = 'time.year' in this case will make the
  function to first get the year means, and then the metric over each year (e.g.
  pixel with maximum temperature of that year).
  
  If you need more flexibility on the calculations you have to do, then the best
  way is to first open the dataset using inspect_cop or inspect_ds and then
  perform the filterings/combinations/calculations using the returned object.
  For that, you must use the xarray functions, which are well described in the
  docs: https://docs.xarray.dev/en/latest/index.html
  """
  
  import matplotlib.pyplot as plt
  import xarray as xr
  import numpy as np
  import pandas as pd
  import geopandas as gpd
  import regionmask
  
  if metric != None and group_by_lonlat == True:
    if ret_type != 'df':
      raise ValueError('When grouping by lonlat, ret_type should be "df". Check arguments.')
    
  if ret_type != 'raster' and ret_type != 'df':
    raise ValueError('ret_type should be one of "raster" or "df"')
  
  #! /usr/bin/env python3
  # -*- coding: utf-8 -*-
  __author__ = "Copernicus Marine User Support Team"
  __copyright__ = "(C) 2022 E.U. Copernicus Marine Service Information"
  __credits__ = ["E.U. Copernicus Marine Service Information"]
  __license__ = "MIT License - You must cite this source"
  __version__ = "202104"
  __maintainer__ = "D. Bazin, E. DiMedio, C. Giordan"
  __email__ = "servicedesk dot cmems at mercator hyphen ocean dot eu"
  
  def copernicusmarine_datastore(ds, user, pwd):
      from pydap.client import open_url
      from pydap.cas.get_cookies import setup_session
      cas_url = 'https://cmems-cas.cls.fr/cas/login'
      session = setup_session(cas_url, user, pwd)
      session.cookies.set("CASTGC", session.cookies.get_dict()['CASTGC'])
      database = ['my', 'nrt']
      url = f'https://{database[0]}.cmems-du.eu/thredds/dodsC/{ds}'
      try:
          data_store = xr.backends.PydapDataStore(open_url(url, session=session)) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits 
      except:
          url = f'https://{database[1]}.cmems-du.eu/thredds/dodsC/{dataset}'
          data_store = xr.backends.PydapDataStore(open_url(url, session=session)) # needs PyDAP >= v3.3.0 see https://github.com/pydap/pydap/pull/223/commits
      return data_store
  
  if isinstance(dataset, list):
    print('Dataset provided', flush = True) 
    DS = dataset[1]
    dataset = dataset[0]
  else:
    print('Loading dataset...', flush = True) 
    data_store = copernicusmarine_datastore(dataset, username, password)
    DS = xr.open_dataset(data_store)
  
  print('Loading shapefile', flush = True)
  shape_area = gpd.read_file(shape)
  
  print('Slicing time and area for region '+shape_var+' '+shape_filter, flush = True)
  sel_area = shape_area[getattr(shape_area, shape_var) == str(shape_filter)]
  
  if plot:
    f, ax = plt.subplots()
    sel_area.plot(ax=ax)
    ax.set(title="MRGID "+str(shape_filter))
    plt.show()
    
  sel_area_lat = [float(sel_area.total_bounds[1]), float(sel_area.total_bounds[3])]
  sel_area_lon = [float(sel_area.total_bounds[0]), float(sel_area.total_bounds[2])]
  
  if DS.attrs["easternmost_longitude"] > 190:
    sel_area_lon[0] = sel_area_lon[0] + 360
    sel_area_lon[1] = sel_area_lon[1] + 360
  
  start_date = time_range[0]
  end_date = time_range[1]
  
  start_depth = depth_range[0]
  end_depth = depth_range[1]
  
  if 'longitude' in DS.dims:
    DS = DS.rename({'longitude' : 'lon', 'latitude' : 'lat'})
  if 'Longitude' in DS.dims:
    DS = DS.rename({'Longitude' : 'lon', 'Latitude' : 'lat'})
  
  if max(DS.lon) > 190:
    print('Converting from 0/360 to -180/180. Check plot to see if it is ok.')
    lon_name = 'lon'  
    DS['_lon_cor'] = xr.where(
        DS[lon_name] > 180,
        DS[lon_name] - 360,
        DS[lon_name])
    DS = (
        DS
        .swap_dims({lon_name: '_lon_cor'})
        .sel(**{'_lon_cor': sorted(DS._lon_cor)})
        .drop(lon_name))
    
    DS = DS.rename({'_lon_cor': lon_name})
    
    if plot:
      to_plot = DS[variable].isel(time=1)
      f, ax = plt.subplots()
      to_plot.plot(ax=ax)
      plt.show()
  
  sliced_data = DS[variable].sel(
    depth=slice(start_depth, end_depth),
    time=slice(start_date, end_date),
    lon=slice(sel_area_lon[0], sel_area_lon[1]),
    lat=slice(sel_area_lat[0], sel_area_lat[1]))
    
  if k_to_celsius:
    sliced_data = sliced_data - 273.15
    
  if mask_by_shape:
    area_mask = regionmask.mask_3D_geopandas(sel_area, sliced_data.lon, sliced_data.lat)
    #area_mask = regionmask.mask_geopandas(sel_area, sliced_data)
    
    sliced_data = sliced_data.where(area_mask)
    #sliced_data = sliced_data.where(area_mask == 0)
    
    if plot:
      sliced_data.isel(time=1, depth=1).plot(figsize=(10, 10))
      plt.show()
  
  else:
    if plot:
      sliced_data.isel(time=1, depth=1).plot(figsize=(10, 10))
      plt.show()
  
  if res_by_month:
    print("Resampling (mean by month).", flush = True)
    sliced_data = sliced_data.resample(time='1MS').mean('time',keep_attrs=True,skipna=False)
    
  if agg_depth:
    sliced_data = sliced_data.mean('depth')
      
  if metric != None:
    print("Calculating metric.", flush = True)
    
    if group_by_lonlat:
      print("Grouping by "+grouper+" and calculating metric over lon-lat.", flush = True)
      sliced_data = sliced_data.groupby(grouper)
      if grouper != 'time':
        sliced_data = sliced_data.mean()
      sli_by = ['lat', 'lon']
    else:
      sli_by = 'time'
      if grouper != 'time':
        sliced_data = sliced_data.groupby(grouper)
    
    if not isinstance(metric, list):
      metric = [metric]
  
    for i in range(len(metric)):
      smetric = metric[i]
      print('Getting metric '+smetric)
      
      if smetric == 'mean':
        int_data = sliced_data.mean(sli_by).rename('mean')
      elif smetric == 'min':
        int_data = sliced_data.min(sli_by).rename('min')
      elif smetric == 'max':
        int_data = sliced_data.max(sli_by).rename('max')
      elif smetric == 'median':
        int_data = sliced_data.median(sli_by).rename('median')
      elif smetric == 'sum':
        int_data = sliced_data.sum(sli_by).rename('sum')
      elif smetric == 'std':
        int_data = sliced_data.std(sli_by).rename('sd')
      elif smetric == 'meanstd':
        int_data_a = sliced_data.mean(sli_by).rename('mean')
        int_data_b = sliced_data.std(sli_by).rename('sd')
        int_data = xr.merge([int_data_a, int_data_b])
      else:
        print("Metric not found, returning mean and standard deviation.", flush = True)
        int_data_a = sliced_data.mean(sli_by).rename('mean')
        int_data_b = sliced_data.std(sli_by).rename('sd')
        int_data = xr.merge([int_data_a, int_data_b])
        smetric = 'meanstd'
      
      if i == 0:
        final_data = int_data
      else:
        final_data = xr.merge([final_data, int_data])
    
    if len(metric) > 1:
      send = "_summarised.nc"
    else:
      send = "_"+smetric+".nc"
      
    if ret_type == 'df':
      final_data = final_data.to_dataframe()
  
  if metric == None:
    final_data = sliced_data
    send = ".nc"
      
  if ret_type == 'raster':
    final_data.to_netcdf(save_path+dataset+"_mrgid"+str(shape_filter)+send)
    print("File saved.", flush = True)
    return save_path+dataset+"_mrgid"+str(shape_filter)+send
  else:
    print("Data retrieving done.", flush = True)
    return final_data

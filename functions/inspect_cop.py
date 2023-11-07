############# eDNA Expeditions - UNESCO World Heritage Marine Sites ############
################# Data analysis of the project - Thermal niche #################
# August of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############### Inspect satellite dataset from Marine Copernicus ###############
def inspect_cop(dataset, username, password, plot = False, variable = None):
  
  """Inspect dataset from Copernicus Marine data products
  
  This function enables to inspect Copernicus satellite data prior to retrieving
  the full dataset or part of that.
  
  Parameters
  ----------
  dataset : str
      Name of the dataset (e.g. METOFFICE-GLO-SST-L4-REP-OBS-SST)
  username : str
      Copernicus username
  password : str
      Copernicus password
  plot : bool
      Should the function plot the maps during the process? Default is False.
  variable : str, optional
      Name of the variable (e.g. analyzed_sst). Should be supplied if plot = True.
      
  Returns
  -------
  xarray.Dataset
    This object can be supplied to the function retrieve_cop
      
  Typical usage example:
  -------
    user = 'myuser'
    pass = 'mypassword'
    shp_path = 'path/to/my_shape.shp'
    out_path = 'data/saved/'
    inspect_cop('METOFFICE-GLO-SST-L4-REP-OBS-SST', 'analyzed_sst',
                  user, pass)
                  
  Depends on the following packages:
  -------                
  matplotlib, xarray
  You may need to install lxml
  """
  
  import matplotlib.pyplot as plt
  import xarray as xr
  
  if plot and variable == None:
    raise ValueError('When plot=True, variable should be supplied. Check arguments.')
  
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
  
  print('Loading dataset...', flush = True) 
  data_store = copernicusmarine_datastore(dataset, username, password)
  DS = xr.open_dataset(data_store)
  
  print(DS, flush = True)
  
  if plot:
    sliced = DS[variable].isel(time=1)
    if 'depth' in DS.dims:
      sliced = sliced.isel(depth=1)
    plt.figure()
    sliced.plot()
    plt.show()

  return DS

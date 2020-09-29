import os
import sys
import xarray as xr
import geopandas as gpd
import zipfile

sys.path.append("../../../NHFLO/NHFLOPY")
from modules import mgrid, util

datadir = '../../data'
modelname = "heel_gebied"
regis_nc = f"regis_ugw.nc"  # regis NetCDF filename

# %% copy more specific files from google drive folder of Artesia
# the waterinfo file
fname = os.path.join(datadir, '20200717_026.zip')
util.download_file_from_google_drive(
    '1ui1_wZJoN2oh2xwbdLHStVVO7qv0rQmj', fname)

# modflow_heel_gebied.zip
file_id = '1AHp8LVo0JwSze58s3vCiEOO8au2b_Oys'
fname = os.path.join(datadir, 'modflow_sfw.zip')
util.download_file_from_google_drive(file_id, fname)
# extract to modflow_sfw
with zipfile.ZipFile(fname, 'r') as zip_ref:
    zip_ref.extractall(os.path.join(datadir, f'modflow_sfw'))
os.remove(fname)

# bathymetry
fname = os.path.join(datadir, 'Bathymetry.zip')
util.download_file_from_google_drive(
    '1k0O6FiOpsjy8-wA0cVgLzuNmAYBwz4IX', fname)
with zipfile.ZipFile(fname, 'r') as zip_ref:
    zip_ref.extractall(datadir)
os.remove(fname)

# %% extent
shp = os.path.join(datadir, f"modflow_sfw/waterareas.shp")
gdf = gpd.read_file(shp)
bounds = gdf.geometry.total_bounds
extent = [bounds[0], bounds[2], bounds[1], bounds[3]]

# %% regis
cs_regis = 100.  # cell size regis
url = r'http://www.dinodata.nl:80/opendap/REGIS/REGIS.nc'

# redefine extent, nrow & ncol (fit to regis)
extent, nrow, ncol = mgrid.fit_extent_to_regis(extent, cs_regis, cs_regis)

# get number of layers
regis_ds_raw = xr.open_dataset(url).sel(x=slice(extent[0], extent[1]),
                                        y=slice(extent[2], extent[3]))
regis_ds_raw = regis_ds_raw[['top', 'bottom', 'kD', 'c', 'kh', 'kv']]

nlay, lay_sel = mgrid.get_lay_from_ml_layers(regis_ds_raw)
regis_ds_raw = regis_ds_raw.sel(layer=lay_sel)

# %% write netcdf
regis_ds_raw.to_netcdf(os.path.join(datadir, regis_nc))

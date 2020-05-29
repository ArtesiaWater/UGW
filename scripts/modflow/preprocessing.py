import os
import sys
import xarray as xr
import geopandas as gpd

sys.path.insert(1, "../../mfutil")
import mgrid

# %% extent
extent = (112000.0, 115200.0, 444800.0, 447000.0)

shp = "../../data/modflow_sfw_schoonhoven/waterareas.shp"
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

nlay, lay_sel = mgrid.get_number_of_layers_from_regis(regis_ds_raw)
regis_ds_raw = regis_ds_raw.sel(layer=lay_sel)

# %% write netcdf
datadir = '../../data'
regis_ds_raw.to_netcdf(os.path.join(datadir, 'regis_ugw_test2.nc'))

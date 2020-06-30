import os
import sys
import xarray as xr
import geopandas as gpd
import zipfile

sys.path.append("../../../../NHFLO/NHFLOPY")
from modules import mgrid, util

datadir = '../../../data'
modelname = "heel_gebied"
regis_nc = f"regis_ugw_{modelname}.nc"  # regis NetCDF filename

# %% copy more specific files from google drive folder of Artesia
if False:
    # the waterinfo file
    fname = os.path.join(datadir, '20200603_044.zip')
    util.download_file_from_google_drive(
        '1jG7KvCSQH1PLGI00HYo4EpjYxTM1UCAr', fname)

    # modflow.zip
    if modelname == 'schoonhoven':
        file_id = '1k0QiOlMLAg5KS-B9Tra4ncpRnBvOylX0'
    elif modelname == 'heel_gebied':
        file_id = '1-OcX9R-j_wBQQWYHbXdQNZaw7-YHdcoM'
    else:
        raise(Exception(f'Modelname {modelname} not supported'))
    fname = os.path.join(datadir, 'modflow.zip')
    util.download_file_from_google_drive(file_id, fname)
    # extract to modflow_sfw_schoonhoven
    with zipfile.ZipFile(fname, 'r') as zip_ref:
        zip_ref.extractall(os.path.join(datadir, f'modflow_sfw_{modelname}'))
    os.remove(fname)

    # bathiemetry
    fname = os.path.join(datadir, 'Bathymetry.zip')
    util.download_file_from_google_drive(
        '1k0O6FiOpsjy8-wA0cVgLzuNmAYBwz4IX', fname)
    with zipfile.ZipFile(fname, 'r') as zip_ref:
        zip_ref.extractall(datadir)

# %% extent
shp = os.path.join(datadir, f"modflow_sfw_{modelname}/waterareas.shp")
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
regis_ds_raw.to_netcdf(os.path.join(datadir, regis_nc))

# %% download DINO
if False:
    import hydropandas as hpd

    oc_dino = hpd.ObsCollection.from_dino(extent=extent, verbose=True)
    oc_dino.to_pickle(os.path.join(datadir, f'oc_dino_{modelname}.pklz'))

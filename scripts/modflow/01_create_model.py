import os
import platform
import sys
from shutil import copyfile
from timeit import default_timer

import flopy as fp
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import shapely
import xarray as xr
from matplotlib.patches import Patch
from tqdm.auto import tqdm

sys.path.insert(2, "../../../NHFLO/NHFLOPY")
from modules import ahn, layer_models, mgrid, mtime, rws, surface_water, util

# mpl.interactive(True)

start = default_timer()

# %% Model settings

# modelnaam en map
model_name = 'ugw_agg_dl'
model_ws = f'../../models/{model_name}'

# method for parsing surface water
agg_method = "de_lange"  # 'individual', 'de_lange', 'area_weighted' or 'max_area'

# extent en gridgrootte
# extent = [115900, 121000, 436600, 442000]  # Schoonhoven
# extent = [151990, 172800, 467700, 489300.]  # Nijkerk
extent = None  # extent is determined from surface-water shape

delr = 500.            # zelfde als dx
delc = 500.            # zelfde als dy

# tijdsdiscretisatie
steady_state = True  # steady state flag, if True, parameters below are not used
steady_start = True  # if True start transient model with steady timestep
nstp = 1  # no. of timesteps per stress period
perlen = 10  # length of timestep in time_units (default is days)
start_time = '2019-01-01'  # start time (after the steady state period)
transient_timesteps = int(365 / 10)  # no. of transient steps

# geef hier paths op en of de cache gebruikt moet worden
datadir = '../../data'
figdir = os.path.join(model_ws, 'figure')
use_cache = True
cachedir = os.path.join(model_ws, 'cache')

# files
regis_nc = 'regis_ugw.nc'
water_shp = os.path.join(datadir, "modflow_sfw", "waterareas.shp")
lines_shp = os.path.join(datadir, "modflow_sfw", "waterlines.shp")
rivobs_fname = os.path.join(datadir, '20200717_026.zip')

# verander dit niet
if not os.path.exists(model_ws):
    os.makedirs(model_ws)

if not os.path.exists(figdir):
    os.mkdir(figdir)

if not os.path.exists(cachedir):
    os.mkdir(cachedir)

surfwat_pkgs = []  # collect surface water packages in this list

try:
    os.remove(os.path.join(model_ws, os.path.basename(__file__)))
except:
    pass
try:
    copyfile(__file__, os.path.join(model_ws, os.path.basename(__file__)))
except NameError:
    print("Could not copy script to model directory.")

# %% Shapefile (for RIV and extent)

# read water areas shapefile
print("Read waterareas shapefile...")
if extent is None:
    sfw = gpd.read_file(water_shp)
else:
    bbox = (extent[0], extent[2], extent[1], extent[3])
    sfw = gpd.read_file(water_shp, bbox=bbox)

# drop duplicates
mask = sfw["geometry"].apply(lambda geom: geom.wkb)
dupes = mask.drop_duplicates()
sfw = sfw.loc[dupes.index]
print(f"- Dropped {sfw.index.size - dupes.index.size} duplicates!")

# read water lines shapefile
if extent is None:
    slope_lines = gpd.read_file(lines_shp)
else:
    bbox = (extent[0], extent[2], extent[1], extent[3])
    slope_lines = gpd.read_file(lines_shp, bbox=bbox)

# %% time discretization
time_units = 'DAYS'
tsmult = 1.0

model_ds = mtime.get_model_ds_time(model_name,
                                   start_time=start_time,
                                   steady_state=steady_state,
                                   steady_start=steady_start,
                                   time_units=time_units,
                                   transient_timesteps=transient_timesteps,
                                   perlen=perlen,
                                   nstp=nstp,
                                   tsmult=tsmult)

tdis_perioddata = [(model_ds.perlen,
                    model_ds.nstp,
                    model_ds.tsmult)] * model_ds.nper

# %% SIM
# Create the Flopy simulation object
exe_name = '../../tools/mf6'
if platform.system() in 'Windows':
    exe_name += '.exe'

sim = fp.mf6.MFSimulation(sim_name=model_name,
                          exe_name=exe_name,
                          version='mf6',
                          sim_ws=model_ws)

# %% TDIS
tdis = fp.mf6.ModflowTdis(sim,
                          pname='tdis',
                          time_units=model_ds.time_units,
                          nper=model_ds.nper,
                          start_date_time=model_ds.start_time,
                          perioddata=tdis_perioddata)

# %% GWF
model_nam_file = '{}.nam'.format(model_name)
gwf = fp.mf6.ModflowGwf(sim,
                        modelname=model_name,
                        model_nam_file=model_nam_file)

# %% IMS
ims = fp.mf6.ModflowIms(sim,
                        pname='ims',
                        complexity='SIMPLE')

# %% Define modflow grid
if extent is None:
    bounds = sfw.geometry.total_bounds
    extent = (bounds[0], bounds[2], bounds[1], bounds[3])
else:
    mlpol = shapely.geometry.box(extent[0], extent[2], extent[1], extent[3])
    sfw = sfw[sfw.intersects(mlpol)].reset_index(drop=True)
    sfw.geometry = sfw.intersection(mlpol)

length_units = 'METERS'

# redefine extent, nrow & ncol (fit to regis)
extent, nrow, ncol = mgrid.fit_extent_to_regis(list(extent), delr, delc)

# get regis dataset
regis_path = os.path.join(datadir, regis_nc)
regis_ds_raw = xr.open_dataset(regis_path).sel(x=slice(extent[0], extent[1]),
                                               y=slice(extent[2], extent[3]))

nlay, lay_sel = mgrid.get_lay_from_ml_layers(regis_ds_raw)
regis_ds_raw = regis_ds_raw.sel(layer=lay_sel)

# convert regis dataset to grid
regis_ds = util.get_ml_layer_dataset_struc(gridtype='structured',
                                           raw_ds=regis_ds_raw,
                                           extent=extent,
                                           delr=delr,
                                           delc=delc,
                                           interp_method="nearest",
                                           cachedir=cachedir,
                                           fname_netcdf=regis_nc,
                                           use_cache=use_cache)

# %% get model_ds, add idomain, top & bot
model_ds = mgrid.update_model_ds_from_ml_layer_ds(model_ds, regis_ds,
                                                  keep_vars=['x', 'y'],
                                                  verbose=True)
model_ds = mgrid.add_idomain_from_bottom_to_dataset(regis_ds['bottom'],
                                                    model_ds)
model_ds = layer_models.add_top_bot_to_model_ds(regis_ds,
                                                model_ds,
                                                gridtype='structured')

# %% flow parameters
confined = True
fill_value_kh = 1.
fill_value_kv = 0.1
anisotropy = 10

if confined:
    icelltype = 0
else:
    icelltype = 1

model_ds = layer_models.add_kh_kv_from_ml_layer_to_dataset(regis_ds,
                                                           model_ds,
                                                           anisotropy,
                                                           fill_value_kh,
                                                           fill_value_kv)


# %% DIS

# update idomain on adjusted tops and bots
model_ds['thickness'] = mgrid.get_thickness_from_topbot(
    model_ds['top'], model_ds['bot'])
model_ds['idomain'] = mgrid.update_idomain_from_thickness(
    model_ds['idomain'], model_ds['thickness'], 1)
model_ds['first_active_layer'] = mgrid.get_first_active_layer_from_idomain(
    model_ds['idomain'])

# Create DIS package
dis = fp.mf6.ModflowGwfdis(gwf,
                           pname='dis',
                           length_units=length_units,
                           xorigin=model_ds.extent[0],
                           yorigin=model_ds.extent[2],
                           angrot=0.0,
                           nlay=model_ds.dims['layer'],
                           nrow=model_ds.dims['y'],
                           ncol=model_ds.dims['x'],
                           delr=model_ds.delr,
                           delc=model_ds.delc,
                           top=model_ds['top'].data,
                           botm=model_ds['bot'].data,
                           idomain=model_ds['idomain'].data,
                           filename=f'{model_name}.dis')

# %% Add information about the surface level (also bathymetry)
print("\nLoad AHN and bathymetry data...")

# add the surface level of each grid cell
model_ds['area'] = (('y', 'x'), mgrid.get_surface_area(gwf))

# get the minimum ahn level in each cell
ahn_fname = ahn.get_ahn_within_extent(model_ds.attrs['extent'],
                                      return_fname=True, cache_dir=cachedir)
resampling = rasterio.enums.Resampling.min
model_ds['ahn_min'] = mgrid.raster_to_quadtree_grid(ahn_fname, model_ds,
                                                    resampling=resampling)
resampling = rasterio.enums.Resampling.average
model_ds['ahn_average'] = mgrid.raster_to_quadtree_grid(ahn_fname, model_ds,
                                                        resampling=resampling)
resampling = rasterio.enums.Resampling.max
model_ds['ahn_max'] = mgrid.raster_to_quadtree_grid(ahn_fname, model_ds,
                                                    resampling=resampling)

# read Bathymetry of river data
fname = os.path.join(datadir, 'Bathymetry', 'bathymetry_masks.shp')
bathshp = gpd.read_file(fname)
bathshp["FILE"] = bathshp["FILE"].apply(lambda fp: fp.replace(
    "\\", "/") if isinstance(fp, str) else None)
extent_polygon = shapely.geometry.box(model_ds.attrs['extent'][0],
                                      model_ds.attrs['extent'][2],
                                      model_ds.attrs['extent'][1],
                                      model_ds.attrs['extent'][3])
mask = bathshp.intersects(extent_polygon)
bathshp = bathshp[mask]
bath = xr.full_like(model_ds['top'], np.NaN)
for file in bathshp['FILE'].dropna():
    fname = os.path.join(datadir, file.replace('../data/sources/', ''))
    # get the minimum bathemetry-level in each cell
    resampling = rasterio.enums.Resampling.min
    zt = mgrid.raster_to_quadtree_grid(fname, model_ds, resampling=resampling)
    # update bath when zt is lower
    bath = bath.where(np.isnan(zt) | (bath < zt), zt)
# apparently bathemetry is in mm (need to check if this is always the case)
model_ds['bathymetry'] = bath = bath / 1000.

# TODO: ensure ahn option always works
# fill bathemetry by ahn, so there is allways data
model_ds['bathymetry'] = model_ds['bathymetry'].fillna(model_ds['ahn_min'])
# and otherwise fill by the top of the model
model_ds['bathymetry'] = model_ds['bathymetry'].fillna(model_ds['top'])

# %% NPF
npf = fp.mf6.ModflowGwfnpf(gwf,
                           pname='npf',
                           icelltype=icelltype,
                           k=model_ds['kh'].data,
                           k33=model_ds['kv'].data,
                           save_flows=True,
                           save_specific_discharge=True)

# %% STO
sy = 0.2
ss = 1e-5

if not model_ds.steady_state:

    if model_ds.steady_start:
        sts_spd = {0: True}
        trn_spd = {1: True}
    else:
        sts_spd = None
        trn_spd = {0: True}

    sto = fp.mf6.ModflowGwfsto(gwf,
                               pname='sto',
                               save_flows=True,
                               iconvert=1,
                               ss=ss,
                               sy=sy,
                               steady_state=sts_spd,
                               transient=trn_spd)

# %% IC
starting_head = model_ds["top"].data

# Create the initial conditions array
layer_store_type = [
    fp.mf6.data.mfdatastorage.DataStorageType.internal_array
]

starting_head = fp.mf6.ModflowGwfic.strt.empty(gwf,
                                               layered=False,
                                               data_storage_type_list=layer_store_type,
                                               default_value=starting_head)

# Create IC package
ic = fp.mf6.ModflowGwfic(gwf,
                         pname='ic',
                         strt=1.0)

# %% RCH
rech = 1e-3  # m/day

rch = fp.mf6.ModflowGwfrcha(gwf,
                            pname="rch",
                            recharge=rech)

# %% Surface water
print("\nIntersect surface water data with modelgrid...")

# get boundary type
sfw["bc"] = surface_water.get_bc_type(sfw)

# intersect with grid (or load from cache if cached)
sfw_grid = util.get_cache_gdf(use_cache, cachedir, "sfw_grid.pkl",
                              surface_water.gdf2grid, model_ds=model_ds,
                              check_grid=True, check_time=True,
                              get_args=(sfw, gwf), method="vertex",
                              verbose=True)

sfw_grid.reset_index(inplace=True)

print(f"\nSFW: Total features source data: {sfw.index.size}")
print(f"\nSFW: Total features after grid intersect: {sfw_grid.index.size}")

# %% Get slope information

# get water areas per gridcell with slope
mask_slope = sfw_grid["bc"] == "riv_slp"
sfw_slope = sfw_grid.loc[mask_slope].copy()

# update source data boundary condition based on missing data
m1a = (sfw["bc"] == "riv_slp") & (sfw["src_id_wla"].isna())
m2a = (sfw["bc"] == "riv_slp") & (sfw["BL"].isna())
ma = m1a | m2a
sfw.loc[ma, "bc"] = "none"

print(f"\nRIV_slp: Parsing {sfw_slope.index.size} features")

if not sfw_slope.empty:

    sfw_slope["name_largest"] = sfw_slope.src_id_wla.apply(
        lambda s: s.replace(".", "_")
        .replace(" ", "_")
        .replace("/", "_") if isinstance(s, str) else s)

    # loop over
    for idx, linedf in tqdm(slope_lines.iterrows(),
                            total=slope_lines.index.size,
                            desc="RIV w slope"):

        line = linedf.geometry  # line for slope
        src_id_wl = linedf.src_id_wl  # name of water areas to consider

        # get associated water areas
        idf = sfw_slope.loc[sfw_slope.src_id_wl == src_id_wl]

        # if any water areas are associated with line
        if not idf.empty:

            # loop per cellid to determine rbot, stage
            gr = idf.groupby(by="cellid")
            for cid, group in gr:

                # interpolate rbot at centroid
                weights = surface_water.calc_interpolation_weights(
                    group.centroid, line)
                z_arr = np.c_[linedf.BLU, linedf.BLD].T
                rbot = weights.dot(z_arr)

                sfw_slope.loc[group.index, "BL"] = rbot.squeeze()

                # stage, if lower than rbot, use rbot
                for col in ["ZP", "WP"]:
                    m = ((sfw_slope[col] < sfw_slope["BL"]) & mask)
                    idx = group.loc[group[col] < rbot.squeeze()].index
                    sfw_slope.loc[idx, col] = sfw_slope.loc[idx, "BL"]

# update dataset
maskdrn = sfw_slope.SUB == 0
sfw_slope[maskdrn,"bc"] = "riv_drn"
maskinf = sfw_slope.SUB == 1
sfw_slope[maskinf,"bc"] = "riv_inf"
sfw_grid.update(sfw_slope)
print("Features added to RIV_inf dataset.")

# %% RIV for surface water with water supply

# get the bottom height from the bathymetry-data
mask_bath = sfw_grid.has_bath == 1
row, col = zip(*sfw_grid.loc[mask_bath, 'cellid'])
sfw_grid.loc[mask_bath, 'BL'] = model_ds['bathymetry'].values[row, col]

mask_riv = sfw_grid["bc"] == "riv_inf"
sfw_riv_inf = sfw_grid.loc[mask_riv]

# check which features have all requisite information
no_wla = sfw_riv_inf.src_id_wla.isna()
no_bl = np.isnan(sfw_riv_inf.BL)
print(f"\nRIV_inf: Parsing {sfw_riv_inf.index.size} features")
print(f"- {no_wla.sum()} features without water level area info!")
print(f"- {no_bl.sum()} features without bottom level!")
print(f"- Skipping {(no_wla | no_bl).sum()} features!")

# drop features without requisite info
sfw_riv_inf = sfw_riv_inf.loc[~(no_wla | no_bl)]

# update main dataset that these entries are not parsed
m1b = (sfw["bc"] == "riv_inf") & (sfw.src_id_wla.isna())
m2b = (sfw["bc"] == "riv_inf") & (sfw.BL.isna())
mb = m1b | m2b
sfw.loc[mb, "bc"] = "none"

# riv
if not sfw_riv_inf.empty:
    riv_data = surface_water.aggregate_surface_water(
        sfw_riv_inf, model_ds=model_ds, method=agg_method)
    riv_spd = surface_water.build_spd(riv_data, "RIV", model_ds)
else:
    riv_spd = []

if len(riv_spd) > 0:
    riv = fp.mf6.ModflowGwfriv(gwf,
                               stress_period_data=riv_spd,
                               save_flows=True,
                               maxbound=len(riv_spd),
                               pname="RIV_inf",
                               filename=f"{model_name}.inf.riv")
    surfwat_pkgs.append("RIV_inf")
else:
    print("No RIV_inf data!")

# model period
mstart = pd.Timestamp(start_time)
mend = model_ds.time.isel(time=-1).to_pandas()

if not steady_state:
    if len(riv_spd) > 0:
        # RIV timeseries
        tseries_list = surface_water.create_timeseries(riv_data, mstart, mend)

        ts0 = tseries_list[0]
        riv.ts.initialize(filename=f'RIV_{ts0.name}.ts',
                          timeseries=list(
                              zip(ts0.index.to_list(), ts0.to_list())),
                          time_series_namerecord=ts0.name,
                          interpolation_methodrecord='stepwise')
        for its in tqdm(tseries_list[1:], desc="Adding RIV timeseries"):
            riv.ts.append_package(filename=f'RIV_{its.name}.ts',
                                  timeseries=list(
                                      zip(its.index.to_list(), its.to_list())),
                                  time_series_namerecord=its.name,
                                  interpolation_methodrecord='stepwise')

# %% DRN for surface water without water supply

mask_riv_drn = sfw_grid["bc"] == "riv_drn"
sfw_riv_drn = sfw_grid.loc[mask_riv_drn]

# check which features have all requisite information
no_wla = sfw_riv_drn.src_id_wla.isna()
no_bl = np.isnan(sfw_riv_drn.BL)
print(f"\nRIV_drn: Parsing {sfw_riv_drn.index.size} features")
print(f"- {no_wla.sum()} features without water level area info!")
print(f"- {no_bl.sum()} features without bottom level!")
print(f"- Skipping {(no_wla | no_bl).sum()} features!")

# select features with requisite info
sfw_riv_drn = sfw_riv_drn.loc[~(no_wla | no_bl)]

# update main dataset that these entries are not parsed
m1c = (sfw["bc"] == "riv_drn") & (sfw.src_id_wla.isna())
m2c = (sfw["bc"] == "riv_drn") & (sfw.BL.isna())
mc = m1c | m2c
sfw.loc[mc, "bc"] = "none"

# drn
if not sfw_riv_drn.empty:
    drn_data = surface_water.aggregate_surface_water(
        sfw_riv_drn, model_ds=model_ds, method=agg_method)
    drn_spd = surface_water.build_spd(drn_data.round(5), "DRN", model_ds)
else:
    drn_spd = []

if len(drn_spd) > 0:
    drn = fp.mf6.ModflowGwfdrn(gwf,
                               stress_period_data=drn_spd,
                               save_flows=True,
                               maxbound=len(drn_spd),
                               pname="RIV_drn",
                               filename=f"{model_name}.drn")
    surfwat_pkgs.append("RIV_drn")
else:
    print("No RIV_drn data!")

# model period
mstart = pd.Timestamp(start_time)
mend = model_ds.time.isel(time=-1).to_pandas()

if not steady_state:
    if len(riv_spd) > 0:
        # RIV timeseries
        tseries_list = surface_water.create_timeseries(drn_data, mstart, mend)

        ts0 = tseries_list[0]
        drn.ts.initialize(filename=f'DRN_{ts0.name}.ts',
                          timeseries=list(
                              zip(ts0.index.to_list(), ts0.to_list())),
                          time_series_namerecord=ts0.name,
                          interpolation_methodrecord='stepwise')
        for its in tqdm(tseries_list[1:], desc="Adding RIV timeseries"):
            drn.ts.append_package(filename=f'DRN_{its.name}.ts',
                                  timeseries=list(
                                      zip(its.index.to_list(), its.to_list())),
                                  time_series_namerecord=its.name,
                                  interpolation_methodrecord='stepwise')
# %% GHB (surface water with bathymetry)

mask_rws = sfw_grid.src_wla == "rws_krw"
gdf_rws = sfw_grid[mask_rws].copy()

if False:
    # make a plot of the surface_water
    ax = gdf_rws.plot('src_id_wla', categorical=True, legend=True, colormap='tab20',
                      legend_kwds={'prop': {"size": 7.0}, 'loc': (0.6, 0.34)},
                      figsize=(10, 10))
    df, meta = rws.read_waterinfo_zip(rivobs_fname, freq='d', metadata=True)
    ax.plot(meta.x, meta.y, linestyle='none', marker='o')
    for name, row in meta.iterrows():
        ax.text(row.x, row.y, name)
    ax.figure.tight_layout(pad=0.0)

gdf_rws = gdf_rws.set_index('src_id_wla')

riv2stn = {}
# Amsterdam-Rijnkanaal Betuwepand
riv2stn['NL.1.NL86_5'] = ['Tiel Kanaal']
# Amsterdam-Rijnkanaal Noordpand
riv2stn['NL.1.NL86_6'] = ['Amsterdam Surinamekade', 'Weesp West', 'Maarssen',
                          'Nieuwegein', 'Wijk bij Duurstede kanaal']
# Noordzeekanaal
riv2stn['NL.1.NL87_1'] = ['IJmuiden binnen', 'Buitenhuizen (kilometer 10)',
                          'Amsterdam Surinamekade', 'Weesp West']
# Bedijkte Maas
riv2stn['NL.1.NL91BM'] = ['Lith boven', 'Megen dorp', 'Grave beneden']
# Markermeer
riv2stn['NL.1.NL92_MARKERMEER'] = ['Schellingwouderbrug', 'Hollandse brug']
# 'Randmeren-oost'
riv2stn['NL.1.NL92_RANDMEREN_OOST'] = ['Nijkerk Nuldernauw', 'Elburg']
# 'Randmeren-zuid'
riv2stn['NL.1.NL92_RANDMEREN_ZUID'] = ['Hollandse brug', 'Nijkerk west']
# Nederrijn, Lek
riv2stn['NL.1.NL93_7'] = ['Hagestein boven', 'Culemborg brug',
                          'Amerongen beneden', 'Amerongen boven', 'Grebbe',
                          'Driel beneden']
# Bovenrijn, Waal
riv2stn['NL.1.NL93_8'] = ['Vuren', 'Zaltbommel', 'Tiel Waal', 'Dodewaard',
                          'Nijmegen haven']
# stukje Boven- en Beneden Merwede
riv2stn['NL.1.NL94_2'] = ['Dordrecht', 'Werkendam buiten', 'Vuren']
# Boven- en Beneden Merwede
riv2stn['NL.1.NL94_3'] = ['Dordrecht', 'Werkendam buiten', 'Vuren']
# Oude Maas
riv2stn['NL.1.NL94_4'] = ['Krimpen a/d IJssel', 'Krimpen a/d Lek',
                          'Schoonhoven', 'Hagestein beneden']
# Beneden Maas
riv2stn['NL.1.NL94_5'] = ['Heesbeen', 'Empel beneden', 'Lith dorp']
# Hollandsche IJssel
riv2stn['NL.1.NL94_7'] = ['Krimpen a/d IJssel', 'Gouda brug']
# stukje Amsterdam-Rijnkanaal Noordpand
riv2stn['NL.14.NL86_6'] = ['Maarssen', 'Nieuwegein']
# stukje Boven- en Beneden Merwede
riv2stn['NL.9.NL94_3'] = ['Dordrecht', 'Werkendam buiten', 'Vuren']

# gdfv = rws.get_river_polygons(extent)
# gdfv['geometry'] = gdfv.intersection(box(extent[0], extent[1], extent[2], extent[3]))
# gdfv[~gdfv.is_empty].reset_index().plot('OWMNAAM', legend=True, cmap='tab20')

if not os.path.isfile(rivobs_fname):
    msg = ('Connot find file {0}. Please run below code (after changing '
           'your e-mail adress) to request data and place requested '
           'file in {0}.')
    raise(Exception(msg.format(rivobs_fname)))
    mailadres = "d.brakenhoff@artesia-water.nl"
    surface_water.request_waterinfo_waterlevels(riv2stn, mailadres,
                                                tmin=model_ds.time.values[0],
                                                tmax=(model_ds.time.values[-1]
                                                      + np.timedelta64(perlen, "D")))

surface_water.waterinfo_to_ghb(rivobs_fname, riv2stn, gdf_rws, gwf, model_ds,
                               intersect_method="vertex")
surfwat_pkgs.append("GHB")

# %% OC

# Create the output control package
headfile = f'{model_name}.hds'
head_filerecord = [headfile]
budgetfile = f'{model_name}.cbb'
budget_filerecord = [budgetfile]
saverecord = [('HEAD', 'ALL'),
              ('BUDGET', 'ALL')]

oc = fp.mf6.ModflowGwfoc(gwf,
                         pname='oc',
                         saverecord=saverecord,
                         head_filerecord=head_filerecord,
                         budget_filerecord=budget_filerecord)

# %% timing
postmodel = default_timer()
print(f"\nElapsed time up to model solve: {postmodel-start:.1f} s")
1/0
# %% Write simulation files
sim.write_simulation()

# %% Run model
success, buff = sim.run_simulation()
print("\nModel ran successfully:", success)

if not success:
    msg = "Model did not run succesfully!"
    raise Exception(msg)

# %% End script
end = default_timer()
print(f"\nElapsed time: {end-start:.1f} s")

# %% plot input
plot_input = False

if plot_input:

    # %% Plot model BC
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_aspect("equal", adjustable="box")

    for ilay in np.arange(gwf.modelgrid.nlay):
        mv = fp.plot.PlotMapView(gwf, layer=ilay, ax=ax)
        mv.plot_bc("RIV_DRN", color=(1, 1, 0, 0.5))
        mv.plot_bc("RIV_INF", color=(0, 0, 128 / 255., 0.5))
        mv.plot_bc("GHB", color=(173 / 255., 216 / 255., 230 / 255., 0.6))

    sfw.plot(color="red", alpha=0.5, ax=ax, zorder=3)

    patches = []
    colors = ["red", "lightblue", "yellow", "navy"]
    labels = ["source data", "large rivers", "drainage", "river"]

    for lbl, c in zip(labels, colors):
        patches.append(Patch(facecolor=c, label=lbl, edgecolor=c))

    ax.legend(patches, labels, loc="upper left", framealpha=1.0)

    plt.yticks(rotation=90, va="center")
    ax.set_xlabel("X (m RD)")
    ax.set_ylabel("Y (m RD)")

    ax.set_title("Oppervlaktewater randvoorwaarden in model")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "BC_types_in_gwmodel.png"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # %% Plot BC types in source data
    fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=150)
    ax.set_aspect("equal", adjustable="box")
    maskrws = sfw.bc == "ghb"
    sfw.loc[maskrws].plot(color="lightblue", ax=ax)  # GHB, large rivs
    maskrivdrn = sfw.bc == "riv_drn"
    sfw.loc[maskrivdrn].plot(color="orange", ax=ax)  # DRN, drainage
    maskrivinf = sfw.bc == "riv_inf"
    sfw.loc[maskrivinf].plot(color="navy", ax=ax)  # RIV, can supply water
    maskslope = sfw.bc == "riv_slp"
    sfw.loc[maskslope].plot(color="limegreen", ax=ax)  # RIV, has slope
    masknone = sfw.bc == "none"
    sfw.loc[masknone].plot(color="red", ax=ax)  # NOT PARSED!

    patches = []
    colors = ["lightblue", "orange", "navy", "limegreen", "red"]
    labels = ["large rivers", "drainage",
              "river", "sloping river", "not parsed"]
    for lbl, c in zip(labels, colors):
        patches.append(Patch(facecolor=c, label=lbl, edgecolor=c))

    xmin, ymin, xmax, ymax = sfw.total_bounds

    ax.legend(patches, labels, loc="best")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.yticks(rotation=90, va="center")
    ax.set_xlabel("X (m RD)")
    ax.set_ylabel("Y (m RD)")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "BC_types_in_source_data.png"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # %% No bottom level
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_aspect("equal", adjustable="box")

    sfw.plot(color="gray", alpha=0.8, ax=ax)
    
    m = (masknone) & (sfw.BL.isna())
    sfw.loc[m].plot(color="red", ax=ax)
    
    xmin, ymin, xmax, ymax = sfw.total_bounds
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.yticks(rotation=90, va="center")
    ax.set_xlabel("X (m RD)")
    ax.set_ylabel("Y (m RD)")
    ax.set_title("No 'bottom_level' defined")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "BC_no_bottom_levels.png"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # %% No src_id_wla (no information about summer/winter levels)
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_aspect("equal", adjustable="box")
    sfw.plot(color="gray", alpha=0.8, ax=ax)
    
    m = (masknone) & (sfw.src_id_wla.isna())
    sfw.loc[m].plot(color="red", ax=ax)
    
    xmin, ymin, xmax, ymax = sfw.total_bounds
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.yticks(rotation=90, va="center")
    ax.set_xlabel("X (m RD)")
    ax.set_ylabel("Y (m RD)")
    ax.set_title("No 'src_id_wla' defined (= no water levels)")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "BC_no_waterlevel-areas.png"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)

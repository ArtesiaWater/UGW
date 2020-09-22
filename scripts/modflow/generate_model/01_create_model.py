import os
from shutil import copyfile
from timeit import default_timer
import platform

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from mpl_toolkits.axes_grid1 import make_axes_locatable
import shapely
from tqdm import tqdm

import flopy as fp
# import modules from NHFLO repo (for now)
import sys
sys.path.insert(2, "../../../../NHFLO/NHFLOPY")
from modules import ahn, mgrid, mtime, rws, subsurface, surface_water, util


start = default_timer()
mpl.interactive(True)

# %% Model settings

use_cache = True
model_ws = '../../model/schoonhoven_ind'
model_name = 'schnhvn_ind'

# method
riv_method = "individual"  # individual or aggregated
agg_method = "delange"  # 'delange', 'area' or 'max'
add_riv_slope = True  # add sloping surface water from line shapefile
surfwat_pkgs = []  # collect surface water packages in this list

# grid
extent = None  # extent is determined from surface-water shape
# extent = [105900., 172800., 425400., 489300.] # entire model area
extent = [116000, 120000, 438000, 442000]  # Schoonhoven
delr = 50.            # zelfde als dx
delc = 50.            # zelfde als dy

# geef hier paths op
subname = "schoonhoven"
datadir = '../../../data'
figdir = os.path.join(model_ws, 'figure')
cachedir = os.path.join(model_ws, 'cache')

# files
regis_nc = f'regis_ugw_{subname}.nc'
dinofile = os.path.join(datadir, f'oc_dino_{subname}.pklz')
water_shp = os.path.join(datadir, f"modflow_sfw", "waterareas.shp")
lines_shp = os.path.join(datadir, f"modflow_sfw", "waterlines.shp")
rivobs_fname = os.path.join(datadir, '20200717_026.zip')

# verander dit niet
if not os.path.exists(model_ws):
    os.makedirs(model_ws)

if not os.path.exists(figdir):
    os.mkdir(figdir)

if not os.path.exists(cachedir):
    os.mkdir(cachedir)

f = open(os.path.join(model_ws, f"{model_name}_log.txt"), "w")

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
if extent is None:
    sfw = gpd.read_file(water_shp)
else:
    bbox = (extent[0], extent[2], extent[1], extent[3])
    sfw = gpd.read_file(water_shp, bbox=bbox)
# drop duplicates
mask = sfw["geometry"].apply(lambda geom: geom.wkb)
sfw = sfw.loc[mask.drop_duplicates().index]

# read water lines shapefile
try:
    slope_lines = gpd.read_file(lines_shp)
except:
    add_riv_slope = False

if add_riv_slope:
    # TODO: this should not be necessary, but since has_slope isn't correct in
    # source we need to set it manually.
    mask_slopes = sfw.src_id_wl.isin(slope_lines.src_id_wl)
    sfw.loc[mask_slopes, "has_slope"] = 1

# check for odd rbots
mask = (sfw["BL"] > sfw[["ZP", "WP"]].min(axis=1))
if mask.sum() > 0:
    print(f"Warning! RBOT above waterlevel in {mask.sum()} cases!")

# %% Time discretization
# general
time_units = 'DAYS'
nstp = 1
tsmult = 1.0

# steady-state/transient
steady_state = True  # steady state flag
start_time = '2019-01-01'  # start time (after the steady state period)

# no. of transient time steps (only if steady_state is False)
transient_timesteps = int(365 / 10)
steady_start = True  # if True start transient model with steady timestep
perlen = 10  # length of timestep in time_units (see below)

# %% time discretization
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
exe_name = r'..\..\..\tools\mf6'
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
regis_path = os.path.join(datadir, f'regis_ugw_heel_gebied.nc')
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
model_ds = subsurface.add_top_bot_to_model_ds(regis_ds,
                                              model_ds,
                                              gridtype='structured')

# %% flow parameters
anisotropy = 10
confined = True
fill_value_kh = 1.
fill_value_kv = 0.1

# berekenen
if confined:
    icelltype = 0
else:
    icelltype = 1

model_ds = subsurface.add_kh_kv_from_ml_layer_to_dataset(regis_ds,
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
                           filename='{}.dis'.format(model_name))

# %% Add information about the surface level (also bathymetry)

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
                         strt=starting_head)

# %% RCH
rech = 1e-3  # m/day

rch = fp.mf6.ModflowGwfrcha(gwf,
                            pname="rch",
                            recharge=rech)

# %% Surface water


def check_surface_water_shape(sfw):
    print("Checking surface water shape")
    print("--------------------------------")
    mask_riv = ((sfw.SUB == 1) &
                (sfw.has_bath == 0) &
                (sfw.has_slope == 0) &
                (sfw.src_wla != "rws_legger") &
                (~sfw.BL.isna()) &
                (~sfw.src_id_wla.isna()))
    mask_drn = ((sfw.SUB == 0) &
                (sfw.has_bath == 0) &
                (sfw.has_slope == 0) &
                (sfw.src_wla != "rws_legger") &
                (~sfw.BL.isna()) &
                (~sfw.src_id_wla.isna()))
    mask_slope = ((sfw.has_slope == 1) &
                  (~sfw.loc[:, ["ZP", "WP"]].isna().any(axis=1)))
    mask_bath = sfw.has_bath == 1
    mask_all = mask_riv | mask_drn | mask_slope | mask_bath

    maskrivdrn = (sfw.has_bath == 0) & (sfw.has_slope == 0)
    mask_bl_nan = sfw.loc[maskrivdrn, "BL"].isna()
    mask_src_nan = sfw.loc[maskrivdrn, "src_id_wla"].isna()
    mask_zpwp_nan = sfw.loc[sfw.has_slope ==
                            1, ["ZP", "WP"]].isna().any(axis=1)

    print(f"No. of features        : {sfw.index.size}")
    print(f"Not parsed             : {sfw.index.size - mask_all.sum()}")
    print("Problems for RIV/DRN   :")
    print(f"- 'BL' is NaN          : {mask_bl_nan.sum()}")
    print(f"- 'src_id_wla' is NaN  : {mask_src_nan.sum()}")
    print("Problems for RIV slope :")
    print(f"- 'ZP' or 'WP' are NaN : {mask_zpwp_nan.sum()}")


check_surface_water_shape(sfw)

# intersect with grid (or load from cache if cached)
sfw_grid = util.get_cache_gdf(use_cache, cachedir, "sfw_grid.pkl",
                              surface_water.gdf2grid, model_ds=model_ds,
                              check_grid=True, check_time=True,
                              get_args=(sfw, gwf), method="vertex",
                              verbose=True)

# %% RIV/DRN for surface water

# get the bottom height from the bathemetry-data
mask = sfw_grid.has_bath == 1
row, col = zip(*sfw_grid.loc[mask, 'cellid'])
sfw_grid.loc[mask, 'BL'] = model_ds['bathymetry'].values[row, col]

# do not parse if src_id_wla is NaN.
# do not parse if bottom level is NaN
mask_riv = ((sfw_grid.SUB == 1) &
            (sfw_grid.has_slope == 0) &
            (sfw_grid.src_wla != "rws_legger") &
            (~sfw_grid.BL.isna()) &
            (~sfw_grid.src_id_wla.isna()))
mask_drn = ((sfw_grid.SUB == 0) &
            (sfw_grid.has_slope == 0) &
            (sfw_grid.src_wla != "rws_legger") &
            (~sfw_grid.BL.isna()) &
            (~sfw_grid.src_id_wla.isna()))
sfw_riv = sfw_grid.loc[mask_riv]
sfw_drn = sfw_grid.loc[mask_drn]

if riv_method == "aggregated":
    boundnames = False

    # riv
    if not sfw_riv.empty:
        riv_data = surface_water.aggregate_surface_water(
            sfw_riv, model_ds=model_ds, method=agg_method)
        riv_data2 = riv_data.loc[~riv_data.ZP.isna()]
        riv_spd = surface_water.build_spd(riv_data2, "RIV", model_ds, f)
    else:
        riv_spd = []

    # drn
    if not sfw_drn.empty:
        drn_data = surface_water.aggregate_surface_water(
            sfw_drn, model_ds=model_ds, method=agg_method)
        drn_data2 = drn_data.loc[~drn_data.ZP.isna()]
        drn_spd = surface_water.build_spd(drn_data2, "DRN", model_ds, f)
    else:
        drn_spd = []

elif riv_method == "individual":

    cbot = 1.0  # days
    boundnames = False  # option for later

    riv_data = sfw_riv.copy()
    drn_data = sfw_drn.copy()

    # make name empty string if name does not exist
    riv_data.loc[sfw_riv['name'].isna(), 'name'] = ''
    drn_data.loc[drn_data['name'].isna(), 'name'] = ''

    # cond
    riv_data["cond"] = riv_data.area / cbot
    drn_data["cond"] = drn_data.area / cbot

    # stage
    riv_data["name_largest"] = riv_data.src_id_wla
    drn_data["name_largest"] = drn_data.src_id_wla

    # build spd
    riv_spd = surface_water.build_spd(
        riv_data.set_index("cellid"), "RIV", model_ds, f)
    drn_spd = surface_water.build_spd(
        drn_data.set_index("cellid"), "DRN", model_ds, f)

else:
    raise ValueError(
        f"Invalid method for creating RIV package: '{riv_method}'")

if len(riv_spd) > 0:
    riv = fp.mf6.ModflowGwfriv(gwf,
                               stress_period_data=riv_spd,
                               boundnames=boundnames,
                               save_flows=True,
                               maxbound=len(riv_spd))
    surfwat_pkgs.append("RIV")

if len(drn_spd) > 0:
    drn = fp.mf6.ModflowGwfdrn(gwf,
                               stress_period_data=drn_spd,
                               boundnames=boundnames,
                               save_flows=True,
                               maxbound=len(drn_spd))
    surfwat_pkgs.append("DRN")

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

    if len(drn_spd) > 0:
        # DRN timeseries
        tseries_list = surface_water.create_timeseries(drn_data, mstart, mend)

        ts0 = tseries_list[0]
        drn.ts.initialize(filename=f'DRN_{ts0.name}.ts',
                          timeseries=list(
                              zip(ts0.index.to_list(), ts0.to_list())),
                          time_series_namerecord=ts0.name,
                          interpolation_methodrecord='stepwise')
        for its in tqdm(tseries_list[1:], desc="Adding DRN timeseries"):
            drn.ts.append_package(filename=f'DRN_{its.name}.ts',
                                  timeseries=list(
                                      zip(its.index.to_list(), its.to_list())),
                                  time_series_namerecord=its.name,
                                  interpolation_methodrecord='stepwise')

# %% RIV slope

# get water areas per gridcell with slope
mask_slope = ((sfw_grid.has_slope == 1) &
              (~sfw_grid.loc[:, ["ZP", "WP"]].isna().any(axis=1)))
if (mask_slope.sum() > 0) and (add_riv_slope):

    sfw_grid["name_largest"] = sfw_grid.src_id_wla.apply(
        lambda s: s.replace(".", "_")
        .replace(" ", "_")
        .replace("/", "_") if isinstance(s, str) else s)
    sfw_slope_grid = sfw_grid.loc[mask_slope].copy()

    # stress period data
    spd = []

    cbot = 1.0  # days

    # store stage data here
    riv_slp_stage_data = pd.DataFrame(
        index=sfw_slope_grid.name_largest.unique())

    # loop over
    for idx, linedf in tqdm(slope_lines.iterrows(),
                            total=slope_lines.index.size,
                            desc="RIV w slope"):

        line = linedf.geometry  # line for slope
        src_id_wl = linedf.src_id_wl  # name of water areas to consider

        # get associated water areas
        idf = sfw_slope_grid.loc[sfw_slope_grid.src_id_wl == src_id_wl]

        if idf.empty:
            continue

        # loop per cellid to determine rbot, stage, cond
        gr = idf.groupby(by="cellid")
        for cid, group in gr:

            geom = shapely.ops.unary_union(group.geometry.values)

            # rbot
            weights = surface_water.calc_interpolation_weights(
                [geom.centroid], line)
            z_arr = np.c_[linedf.BLU, linedf.BLD].T
            rbot = weights.dot(z_arr)[0, 0]

            # stage
            # set based on mean ZP/WP or rbot if rbot is higher.
            idxmax = group.area.idxmax()
            stage = group.loc[idxmax, ["ZP", "WP"]].mean()
            if np.isnan(stage):
                continue
            if stage < rbot:
                stage = rbot

            # cond
            cond = geom.area / cbot

            # distribute according to rbot and layer info
            lays, conds = surface_water.distribute_cond_over_lays(cond,
                                                                  cid,
                                                                  rbot,
                                                                  model_ds.top,
                                                                  model_ds.bot,
                                                                  model_ds.idomain,
                                                                  model_ds.kh)

            # write SPD
            for lay, cond in zip(lays, conds):
                cellid = (lay,) + cid
                spd.append([cellid, stage, cond, rbot])

    # create riv package
    riv_slp = fp.mf6.ModflowGwfriv(gwf,
                                   filename=f"{model_name}_riv_slp",
                                   stress_period_data=spd,
                                   save_flows=True,
                                   maxbound=len(spd))

# %% GHB (surface water from rws)
mask = sfw_grid.src_wla == "rws_krw"
gdf_rws = sfw_grid[mask].copy()
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

#gdfv = rws.get_river_polygons(extent)
#gdfv['geometry'] = gdfv.intersection(box(extent[0], extent[1], extent[2], extent[3]))
#gdfv[~gdfv.is_empty].reset_index().plot('OWMNAAM', legend=True, cmap='tab20')

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

# close log file
f.close()

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
print(f"Elapsed time up to model solve: {postmodel-start:.1f} s")

# %% Write simulation files
sim.write_simulation()

# %% Run model
success, buff = sim.run_simulation()
print("Model ran successfully:", success)

if not success:
    msg = "Model did not run succesfully!"
    raise Exception(msg)

# %% Post-processing

# open head file
fname = os.path.join(model_ws, headfile)
hds = fp.utils.HeadFile(fname)
h = hds.get_data(kstpkper=(0, 0))
h[h > 1e5] = np.nan

# open budgetfile
fname = os.path.join(model_ws, budgetfile)
cbc = fp.utils.CellBudgetFile(fname)

t = 0
spdis = cbc.get_data(text="SPDIS")

q3D = 0.0

for txt in surfwat_pkgs:
    q = cbc.get_data(kstpkper=(0, 0), text=txt)[t]
    q3Di = cbc.create3D(q, gwf.modelgrid.nlay, gwf.modelgrid.nrow,
                        gwf.modelgrid.ncol)
    q3D += q3Di.data

qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(gwf,
                                                            fname,
                                                            precision="double")

# %% plot heads

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(h[0], cmap="RdBu")
# qv = mapview.plot_vector(qx, qy, istep=1, jstep=1, normalize=False,
#                          scale=2, alpha=0.75, width=.00175,
#                          headwidth=3, headlength=3, headaxislength=2,
#                          pivot="mid")

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("head (m NAP)")
ax.set_title("Grondwaterstand en stroming")
fig.savefig(os.path.join(figdir, "head_layer0.png"),
            bbox_inches="tight", dpi=150)

# %% plot river leakage

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qzmm = q3D / (delr * delc) * 1e3
qzmax = np.max(np.abs(qzmm[0]))
qm = mapview.plot_array(qzmm, cmap="RdBu_r", vmin=-qzmax, vmax=qzmax)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("Kwel/Infiltratie (mm/d)")
fig.savefig(os.path.join(figdir, "river_leakage.png"),
            bbox_inches="tight", dpi=150)

# %% compare to obs

if os.path.isfile(dinofile):

    oc_dino = pd.read_pickle(dinofile)
    oc_dino[['start_date', 'end_date']
            ] = oc_dino.stats.get_first_last_obs_date()
    oc_dino = oc_dino[oc_dino['end_date'] > model_ds.time[0].data]
    oc_dino.loc[:, 'modellayer'] = oc_dino.gwobs.get_modellayers(
        gwf, verbose=False)

    counter = 0
    for i, row in oc_dino.iterrows():
        x = row['x']
        y = row['y']

        lay = row['modellayer']
        if np.isnan(lay):
            continue
        else:
            lay = int(lay)

        # get observations in model period
        obs_plot = row['obs'].loc[model_ds.time.data[0]:model_ds.time.data[-1],
                                  'stand_m_tov_nap']
        obs_plot.columns = [i]

        if obs_plot.dropna().empty:
            continue

        # plot if there are any observations
        if not obs_plot.empty:
            if ((obs_plot.index[-1] > model_ds.time.data[0]) and
                    (obs_plot.index[0] < model_ds.time.data[-1])):

                # create figure
                fig, ax = plt.subplots(figsize=(12, 6))

                # get model heads
                idx = gwf.modelgrid.intersect(x, y)
                hds_idx = hds.get_ts((lay,) + idx)[:, 1]
                hds_df = pd.DataFrame(index=pd.to_datetime(model_ds.time.data),
                                      data={f'model head lay {lay}': hds_idx})

                # plot heads
                hds_df.plot(ax=ax, marker='o', lw=0.2, x_compat=True)

                # plot obs
                obs_day = obs_plot.groupby(obs_plot.index).mean()
                obs_day.plot(ax=ax, marker='.', lw=1.0, markersize=5,
                             label=i, x_compat=True)

                # set axes
                ax.set_xlim(model_ds.time.data[0], model_ds.time.data[-1])
                ax.legend(loc=2)
                ax.grid(b=True)
                ax.set_ylabel("Stijghoogte (m NAP)")

                fig.savefig(os.path.join(figdir, f"model_vs_{i}.png"),
                            bbox_inches="tight", dpi=150)
                plt.close(fig)

# %% End script
end = default_timer()
print(f"Elapsed time: {end-start:.1f} s")

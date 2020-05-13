import os
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from timeit import default_timer

# import custom flopy
import sys
sys.path.insert(1, "/home/david/Github/flopy_db")
import flopy as fp

# import custom module
sys.path.insert(1, "../../mfutil")
import mtime
import mgrid
import subsurface
import util


start = default_timer()

# %% Model settings

use_cache = True
model_ws = r'../../model/test001'
model_name = 'test001'

# geef hier paths op
datadir = r'../../data'
figdir = os.path.join(model_ws, 'figure')
cachedir = os.path.join(model_ws, 'cache')

# verander dit niet
if not os.path.exists(model_ws):
    os.makedirs(model_ws)

if not os.path.exists(figdir):
    os.mkdir(figdir)

if not os.path.exists(cachedir):
    os.mkdir(cachedir)

# %% Time discretization
# general
time_units = 'DAYS'
nstp = 1
tsmult = 1.0

# steady-state/transient
steady_state = False  # steady state flag
start_time = '2020-01-01'  # start time

# no. of transient time steps (only if steady is False)
transient_timesteps = 24
steady_start = False  # if True start transient model with steady timestep
perlen = 30  # length of timestep in time_units (see below)

# %% time discretization
model_ds = mtime.get_model_ts(start_time=start_time,
                              steady_state=steady_state,
                              steady_start=steady_start,
                              time_units=time_units,
                              transient_timesteps=transient_timesteps,
                              perlen=perlen,
                              nstp=nstp,
                              tsmult=tsmult)

tdis_perioddata = [(model_ds.perlen, model_ds.nstp,
                    model_ds.tsmult)] * model_ds.nper

# %% SIM
# Create the Flopy simulation object
sim = fp.mf6.MFSimulation(sim_name=model_name,
                          exe_name='mf6',
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

# extent = (111900.0, 116450.0, 442700.0, 447450.0)  # full riv shape
extent = (112000.0, 115200.0, 444800.0, 447000.0)

delr = 10.            # zelfde als dx
delc = 10.            # zelfde als dy
nlay = 36             # alle actieve regis lagen
angrot = 0            # nog niet geimplementeerd
length_units = 'METERS'

# redefine extent, nrow & ncol (fit to regis)
extent, nrow, ncol = mgrid.fit_extent_to_regis(list(extent), delr, delc)

# get regis dataset
regis_path = os.path.join(datadir, 'regis_ugw_test.nc')
regis_ds_raw = xr.open_dataset(regis_path).sel(x=slice(extent[0], extent[1]),
                                               y=slice(extent[2], extent[3]))

# gebruik dit alleen als je het aantal modellagen wil
# aanpassen n.a.v. het aantal actieve regis lagen
nlay, lay_sel = mgrid.get_number_of_layers_from_regis(regis_ds_raw)
regis_ds_raw = regis_ds_raw.sel(layer=lay_sel)

# convert regis dataset to grid
regis_ds = util.get_regis_dataset(gridtype='structured',
                                  regis_ds_raw=regis_ds_raw,
                                  extent=extent,
                                  delr=delr,
                                  delc=delc,
                                  interp_method="nearest",
                                  cachedir=cachedir,
                                  fname_netcdf='regis_ugw_test.nc',
                                  use_cache=use_cache)

# %% get model_ds, add idomain, top & bot
model_ds = mgrid.update_model_ds_from_regis_ds(model_ds, regis_ds,
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

model_ds = subsurface.add_kh_kv_from_regis_to_dataset(regis_ds,
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
                           angrot=angrot,
                           nlay=model_ds.dims['layer'],
                           nrow=model_ds.dims['y'],
                           ncol=model_ds.dims['x'],
                           delr=model_ds.delr,
                           delc=model_ds.delc,
                           top=model_ds['top'].data,
                           botm=model_ds['bot'].data,
                           idomain=model_ds['idomain'].data,
                           filename='{}.dis'.format(model_name))

# %% NPF

npf = fp.mf6.ModflowGwfnpf(gwf,
                           pname='npf',
                           icelltype=icelltype,
                           k=model_ds['kh'].data,
                           k33=model_ds['kv'].data,
                           save_flows=True,
                           save_specific_discharge=True)

# %% storage package

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

starting_head = 1.0

# Create the initial conditions array
layer_store_type = [
    fp.mf6.data.mfdatastorage.DataStorageType.internal_constant
]

starting_head = fp.mf6.ModflowGwfic.strt.empty(gwf,
                                               layered=False,
                                               data_storage_type_list=layer_store_type,
                                               default_value=1.0)

# Create IC package
ic = fp.mf6.ModflowGwfic(gwf,
                         pname='ic',
                         strt=starting_head)

# %% RCH

rech = 1e-3  # m/day

rch = fp.mf6.ModflowGwfrcha(gwf,
                            pname="rch",
                            recharge=rech)


# %% RIV

# read shapefile
water_shp = os.path.join(datadir, "modflow_sfw", "waterareas.shp")
sfw = gpd.read_file(water_shp)

# check implausible rbots
mask = (sfw["BL"] > sfw[["ZP", "WP"]].min(axis=1))
if mask.sum() > 0:
    print(f"Warning! RBOT above waterlevel in {mask.sum()} cases!")
    print("... setting RBOT 1 m below lowest water level")
    sfw.loc[mask, "BL"] = sfw.loc[mask, ["ZP", "WP"]].min() - 1.0

# intersection
ix = fp.utils.GridIntersect(gwf.modelgrid, method="vertex", rtree="strtree")

# intersect with modelgrid and store attributes
keep_cols = ["CAT", "Z", "ZP", "WP", "BL", "BB", "src_id_wla"]
collect_ix = []
for irow, ishp in tqdm(iterable=sfw.iterrows(), total=sfw.index.size):
    r = ix.intersect_polygon(ishp.geometry)
    idf = gpd.GeoDataFrame(r, geometry="ixshapes")
    # add attributes
    for icol in keep_cols:
        idf[icol] = ishp[icol]
        idf["ishp"] = irow  # id to original shape
    if ishp["name"] is not None:
        idf["name"] = ishp["name"]
    else:
        idf["name"] = ""
    collect_ix.append(idf)

sfw_grid = pd.concat(collect_ix, axis=0)
sfw_grid = sfw_grid.reset_index(drop=True).astype({"BL": np.float,
                                                   "BB": np.float})

# individual method
spd = []
cbot = 1.0
for i, row in sfw_grid.iterrows():

    cid = (0,) + row["cellids"]
    if steady_state:
        stage = row[["ZP", "WP"]].mean()  # mean level summer/winter
    else:
        stage = row["src_id_wla"]
    cond = row["areas"] / cbot
    rbot = row["BL"] if row["BL"] < row["WP"] else row["WP"] - 1.0
    name = row["name"].replace(" ", "_")

    spd.append([cid, stage, cond, rbot, name])

riv = fp.mf6.ModflowGwfriv(gwf,
                           stress_period_data=spd,
                           boundnames=True,
                           save_flows=True,
                           maxbound=len(spd))

# build obs data
riv_obs = {'riv_flows.csv': [
    ('H_IJssel', 'RIV', 'Gekanaliseerde_Hollandsche_IJssel'),
    ('Gr_Keulevaart', 'RIV', 'Tiendwegwetering_Groot_Keulevaart'),
    ('Vlist', 'RIV', 'Vlist')]}

# initialize obs package
riv.obs.initialize(filename=f'{model_name}.riv.obs',
                   digits=9, print_input=True,
                   continuous=riv_obs)

# model period
mstart = pd.Timestamp(start_time)
mend = model_ds.time.isel(time=-1).to_pandas()

if not steady_state:
    tseries_list = []

    for peilvak in sfw.src_id_wla.unique():

        peilen = sfw.loc[sfw.src_id_wla == peilvak, ["ZP", "WP"]].iloc[0]

        dt = pd.date_range(mstart, mend, freq="MS")
        dt_apr_oct = [i for i in dt if i.month in [4, 10]]
        doffset = pd.tseries.offsets.DateOffset(months=6)
        dt_apr_oct.insert(0, dt_apr_oct[0] - doffset)
        dt_apr_oct.append(dt_apr_oct[-1] + doffset)
        dt = pd.DatetimeIndex(dt_apr_oct)
        dtnum = ((dt - mstart).days).to_numpy()
        dtnum[dtnum < 0] = 0.0
        ts = pd.Series(index=dtnum, dtype=float)

        ts.where(dt.month == 4, peilen.ZP, inplace=True)
        ts.where(dt.month == 10, peilen.WP, inplace=True)
        ts.name = peilvak
        tseries_list.append(ts)

    ts0 = tseries_list[0]
    riv.ts.initialize(filename=f'{ts0.name}.ts',
                      timeseries=list(zip(ts0.index.to_list(), ts0.to_list())),
                      time_series_namerecord=ts0.name,
                      interpolation_methodrecord='stepwise')
    for its in tseries_list[1:]:
        riv.ts.append_package(filename=f'{its.name}.ts',
                              timeseries=list(
                                  zip(its.index.to_list(), its.to_list())),
                              time_series_namerecord=its.name,
                              interpolation_methodrecord='stepwise')

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

# %% Write simulation files
sim.write_simulation()

# %% Run model
success, buff = sim.run_simulation()
print("Model ran successfully:", success)

# %% plot riv

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.interactive(True)

# open head file
fname = os.path.join(model_ws, headfile)
hds = fp.utils.HeadFile(fname)
h = hds.get_data(kstpkper=(0, 0))

# open budgetfile
fname = os.path.join(model_ws, budgetfile)
cbc = fp.utils.CellBudgetFile(fname)

spdis = cbc.get_data(text="SPDIS")
qriv = cbc.get_data(kstpkper=(0, 0), text="RIV")[0]
qriv3D = cbc.create3D(qriv, 1, gwf.modelgrid.nrow, gwf.modelgrid.ncol)

qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(gwf,
                                                            fname,
                                                            precision="double")

fname = os.path.join(model_ws, f"riv_flows.csv")
obs = fp.utils.Mf6Obs(fname, isBinary=False, verbose=True)
# obs.get_data()

# %% plot heads

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(h[0], cmap="RdBu")
qv = mapview.plot_vector(qx, qy, istep=2, jstep=2, normalize=False,
                         scale=2, alpha=0.75, width=.00175,
                         headwidth=3, headlength=3, headaxislength=2,
                         pivot="mid")

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("head (m NAP)")
ax.set_title("Grondwaterstand en stroming")
fig.savefig(os.path.join(figdir, "head_layer0_quiver.png"),
            bbox_inches="tight", dpi=150)

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(h[0], cmap="RdBu")
qv = mapview.plot_vector(qx, qy, istep=2, jstep=2, normalize=True,
                         width=.001, headwidth=4, headlength=5,
                         headaxislength=4,
                         pivot="mid")

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("head (m NAP)")
ax.set_title("Grondwaterstand en stroming (genormaliseerd)")
fig.savefig(os.path.join(figdir, "head_layer0_norm-quiver.png"),
            bbox_inches="tight", dpi=150)

# %% plot river leakage

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(qriv3D / (delr * delc) * 1e3, cmap="RdBu_r")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("Kwel/Infiltratie (mm/d)")
fig.savefig(os.path.join(figdir, "river_leakage.png"),
            bbox_inches="tight", dpi=150)

end = default_timer()
print(f"Elapsed time: {end-start:.1f} s")

import sys
import os
import warnings
from timeit import default_timer

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import rasterio
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
import flopy as fp  # latest develop branch

# import modules from NHFLO repo (for now)
sys.path.insert(2, "../../../../NHFLO/NHFLOPY")
from modules import mgrid, mtime, subsurface, util, rws, surface_water, ahn
from utils import de_lange

start = default_timer()

# %% Model settings

use_cache = True
model_ws = r'../../model/test003'
model_name = 'test003'

# method
riv_method = "aggregated"
delange = True

# geef hier paths op
datadir = r'../../../data'
figdir = os.path.join(model_ws, 'figure')
cachedir = os.path.join(model_ws, 'cache')

# verander dit niet
if not os.path.exists(model_ws):
    os.makedirs(model_ws)

if not os.path.exists(figdir):
    os.mkdir(figdir)

if not os.path.exists(cachedir):
    os.mkdir(cachedir)

# %% Shapefile (for RIV and extent)
# read shapefile
water_shp = os.path.join(datadir, "modflow_sfw_schoonhoven", "waterareas.shp")
sfw = gpd.read_file(water_shp)

# %% Time discretization
# general
time_units = 'DAYS'
nstp = 1
tsmult = 1.0

# steady-state/transient
steady_state = False  # steady state flag
start_time = '2019-01-01'  # start time (after the steady state period)

# no. of transient time steps (only if steady is False)
transient_timesteps = int(365 / 10)
steady_start = True  # if True start transient model with steady timestep
perlen = 10  # length of timestep in time_units (see below)

# %% time discretization
model_ds = mtime.get_model_ds_time(start_time=start_time,
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

# extent = (111900.0, 116450.0, 442700.0, 447450.0)
# extent = (112000.0, 115200.0, 444800.0, 447000.0)
bounds = sfw.geometry.total_bounds
extent = (bounds[0], bounds[2], bounds[1], bounds[3])

# geef hier waarden op
delr = 50.            # zelfde als dx
delc = 50.            # zelfde als dy
angrot = 0            # nog niet geimplementeerd
length_units = 'METERS'

# redefine extent, nrow & ncol (fit to regis)
extent, nrow, ncol = mgrid.fit_extent_to_regis(list(extent), delr, delc)

# get regis dataset
regis_path = os.path.join(datadir, 'regis_ugw_test2.nc')
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
                                  fname_netcdf='regis_ugw_test2.nc',
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

# %% Add information about the surface level (also bathymetry)
# add the surface level of each grid cell
model_ds['area'] = (('y', 'x'), mgrid.get_surface_area(gwf))

# get the minimum ahn level in each cell
ahn_fname = ahn.get_ahn_within_extent(model_ds.attrs['extent'],
                                      return_fname=True)
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
extent_polygon = surface_water.extent2polygon(model_ds.attrs['extent'])
mask = bathshp.intersects(extent_polygon)
bathshp = bathshp[mask]
bath = xr.full_like(model_ds['top'], np.NaN)
for file in bathshp['FILE']:
    fname = os.path.join(datadir, file.replace('..\\data\\sources\\', ''))
    # get the minimum bathemetry-level in each cell
    resampling = rasterio.enums.Resampling.min
    zt = mgrid.raster_to_quadtree_grid(fname, model_ds, resampling=resampling)
    # update bath when zt is lower
    bath = bath.where(np.isnan(zt) | (bath < zt), zt)
# apparently bathemetry is in mm (need to check if this is always the case)
model_ds['bathymetry'] = bath = bath / 1000.

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

mask_bathymetry = ((sfw.admin != "RWS") &
                   (sfw.src_id_wla != "NL.9.Lek") &
                   (sfw.src_id_wla != "NL.39.Lek") &
                   (~sfw.src_id_wla.isna()))
sfw = sfw.loc[mask_bathymetry]

# check implausible rbots
mask = (sfw["BL"] > sfw[["ZP", "WP"]].min(axis=1))
if mask.sum() > 0:
    print(f"Warning! RBOT above waterlevel in {mask.sum()} cases!")
    print("... setting RBOT 1 m below lowest water level")
    sfw.loc[mask, "BL"] = (sfw.loc[mask, ["ZP", "WP"]].min() - 1.0).values

# aggregated
if riv_method == "aggregated":

    boundnames = False
    sfw_grid = surface_water.gdf2grid(sfw, gwf, method="vertex")

    # Post process intersection result
    gr = sfw_grid.groupby(by="cellid")
    calc_cols = ["ZP", "WP"]
    mdata = pd.DataFrame(index=gr.groups.keys())
    for igr, group in tqdm(gr):
        idf = group

        for icol in calc_cols:
            # area-weighted
            mdata.loc[igr, icol] = \
                (idf.area * idf[icol]).sum(skipna=False) / idf.area.sum()

        # lowest rbot
        lowest_rbot = idf["BL"].min()
        lowest_lvl = mdata.loc[igr, calc_cols].min()
        if np.isnan(lowest_lvl) or np.isnan(lowest_rbot):
            lowest_rbot = np.nan
            raise ValueError("RBOT is NaN!")
        else:
            if lowest_rbot >= lowest_lvl:
                lowest_rbot = lowest_lvl - 0.5
        mdata.loc[igr, "BL"] = lowest_rbot

        # estimate length from polygon (for shapefactor > 4)
        shape_factor = idf.length / np.sqrt(idf.area)
        len_est1 = (
            idf.length - np.sqrt(idf.length**2 - 16 * idf.area)) / 4
        len_est2 = (
            idf.length + np.sqrt(idf.length**2 - 16 * idf.area)) / 4
        len_est = pd.concat([len_est1, len_est2], axis=1).max(axis=1)

        # estimate length from minimum rotated rectangle (for shapefactor < 4)
        min_rect = idf.geometry.apply(lambda g: g.minimum_rotated_rectangle)
        xy = min_rect.apply(lambda g: np.sqrt(
            (np.array(g.exterior.xy[0]) - np.array(g.exterior.xy[0][0]))**2 +
            (np.array(g.exterior.xy[1]) - np.array(g.exterior.xy[1][0]))**2))
        len_est3 = xy.apply(lambda a: np.partition(a.flatten(), -2)[-2])
        
        # update length estimate where shape factor is lower than 4
        len_est.loc[shape_factor < 4] = len_est3.loc[shape_factor < 4]

        mdata.loc[igr, "len_estimate"] = len_est.sum()

        # area
        mdata.loc[igr, "area"] = idf.area.sum()

        # De Lange Params
        laytop = model_ds.top.isel(x=igr[1], y=igr[0])
        laybot = model_ds.bot.isel(x=igr[1], y=igr[0])

        A = model_ds.delr * model_ds.delc  # cell area
        H0 = laytop.data - laybot.data[0]  # layer thickness
        kv = model_ds['kv'].data[0, igr[0], igr[1]]
        kh = model_ds['kh'].data[0, igr[0], igr[1]]
        c1 = 0.0  # resistance between phreatic layer and regional aquifer

        li = mdata.loc[igr, "len_estimate"]  # length
        B = mdata.loc[igr, "area"] / li  # width
        c0 = 1.0  # river bottom resistance
        p = idf.loc[idf.area.idxmax(), ["ZP", "WP"]].mean()  # waterlevel
        N = 1e-3  # recharge

        pstar, cstar, cond = de_lange(A, H0, kv, kh, c1, li, B, c0, p, N)
        if cond < 0 :
            warnings.warn("Calculated conductance (De Lange) is < 0!")
        mdata.loc[igr, "pstar"] = pstar
        mdata.loc[igr, "cstar"] = cstar
        mdata.loc[igr, "cond"] = cond

        # wla of largest surface water feature
        name_largest = idf.loc[idf.area.idxmax(), "src_id_wla"]
        mdata.loc[igr, "name_largest"] = name_largest

    spd = []
    cbot = 1.0  # bottom resistance, days

    for cellid, row in mdata.iterrows():

        rbot = row["BL"]
        laytop = model_ds.top.isel(x=cellid[1], y=cellid[0])
        laybot = model_ds.bot.isel(x=cellid[1], y=cellid[0])

        layers = []
        if laytop.data < rbot:
            layers = [0]  # weirdness, rbot above top of model
        else:
            layers = [0]
        if (laybot.data > rbot).sum() > 0:
            layers += list(range(1, (laybot.data > rbot).sum() + 1))

        if steady_state:
            stage = row[["ZP", "WP"]].mean()  # mean level summer/winter
        else:
            stage = row["name_largest"]
            if stage is None:
                # if delange:
                #     stage = row["pstar"]
                # else:
                stage = row[["ZP", "WP"]].mean()
            elif isinstance(stage, str):
                stage = stage.replace(".", "_")
            elif stage.isna():
                continue
        if delange:
            cond = row["cond"]
        else:
            cond = row["area"] / cbot
        # rbot = row["BL"] if row["BL"] < stage else stage - 1.0

        for nlay in layers:
            cid = (nlay,) + cellid
            spd.append([cid, stage, cond, rbot])

elif riv_method == "individual":
    boundnames = True

    # intersection
    ix = fp.utils.GridIntersect(gwf.modelgrid, method="vertex")

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

else:
    raise ValueError(
        f"Invalid method for creating RIV package: '{riv_method}'")


riv = fp.mf6.ModflowGwfriv(gwf,
                           stress_period_data=spd,
                           boundnames=boundnames,
                           save_flows=True,
                           maxbound=len(spd))

if boundnames:
    # build obs data
    riv_obs = {'riv_flows.csv': [
        ('H_IJssel', 'RIV', 'Gekanaliseerde_Hollandsche_IJssel'),
        ('Gr_Keulevaart', 'RIV', 'Tiendwegwetering_Groot_Keulevaart'),
        ('Vlist', 'RIV', 'Vlist')]}

    # initialize obs package
    riv.obs.initialize(filename=f'{model_name}.riv.obs',
                       digits=9,
                       print_input=True,
                       continuous=riv_obs)

# model period
mstart = pd.Timestamp(start_time)
mend = model_ds.time.isel(time=-1).to_pandas()

if not steady_state:
    tseries_list = []

    for peilvak in sfw.src_id_wla.unique():
        if peilvak is None:
            continue
        peilen = sfw.loc[sfw.src_id_wla == peilvak, ["ZP", "WP"]].iloc[0]

        if peilen.isna().any():
            continue

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
        ts.name = peilvak.replace(".", "_")
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

# %% RIV2
fname = os.path.join(datadir, '20200603_044.zip')

riv2stn = {}
riv2stn['Amsterdam-Rijnkanaal Betuwepand'] = ['Tiel Kanaal']
riv2stn['Amsterdam-Rijnkanaal Noordpand'] = ['Amsterdam Surinamekade',
                                             'Weesp West', 'Maarssen',
                                             'Nieuwegein',
                                             'Wijk bij Duurstede kanaal']
riv2stn['Bedijkte Maas'] = ['Lith boven', 'Megen dorp', 'Grave beneden']
riv2stn['Beneden Maas'] = ['Heesbeen', 'Empel beneden', 'Lith dorp']
riv2stn['Bovenrijn, Waal'] = ['Vuren', 'Zaltbommel', 'Tiel Waal', 'Dodewaard',
                              'Nijmegen haven']
riv2stn['Boven- en Beneden Merwede'] = ['Dordrecht',
                                        'Werkendam buiten',
                                        'Vuren']
riv2stn['Dordtse Biesbosch'] = ['Moerdijk', 'Werkendam buiten']
riv2stn['Hollandsche IJssel'] = ['Krimpen a/d IJssel', 'Gouda brug']
riv2stn['Markermeer'] = ['Schellingwouderbrug', 'Hollandse brug']
riv2stn['Nederrijn, Lek'] = ['Hagestein boven', 'Culemborg brug',
                             'Amerongen beneden', 'Amerongen boven',
                             'Grebbe', 'Driel beneden']
riv2stn['Noordzeekanaal'] = ['IJmuiden binnen', 'Buitenhuizen (kilometer 10)',
                             'Amsterdam Surinamekade', 'Weesp West']
riv2stn['Oude Maas'] = ['Krimpen a/d IJssel', 'Krimpen a/d Lek',
                        'Schoonhoven', 'Hagestein beneden']
riv2stn['Randmeren-zuid'] = ['Hollandse brug', 'Nijkerk west']
riv2stn['Randmeren-oost'] = ['Nijkerk Nuldernauw', 'Elburg']

gdfv = rws.get_river_polygons(extent)
gdfl = rws.get_river_lines(extent)
if not os.path.isfile(fname):
    msg = ('Connot find file {0}. Please run below code (after changing '
           'your e-mail adress) to request data and place requested '
           'file in {0}.')
    raise(Exception(msg.format(fname)))
    surface_water.request_waterinfo_waterlevels(riv2stn, '<mail_adress>',
                                                tmin=model_ds.time.values[0],
                                                tmax=model_ds.time.values[-1])
surface_water.waterinfo_to_ghb(fname, riv2stn, gdfv, gwf, gdfl=gdfl,
                               steady_start=steady_start,
                               intersect_method="vertex")

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
qghb = cbc.get_data(kstpkper=(0, 0), text="GHB")[0]
qghb3D = cbc.create3D(qghb, 1, gwf.modelgrid.nrow, gwf.modelgrid.ncol)

q3D = qriv3D.data + qghb3D.data

qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(gwf,
                                                            fname,
                                                            precision="double")

# %% plot heads and flow quiver

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(h[0], cmap="RdBu")
qv = mapview.plot_vector(qx, qy, istep=1, jstep=1, normalize=False,
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
qv = mapview.plot_vector(qx, qy, istep=1, jstep=1, normalize=True,
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
qm = mapview.plot_array(q3D / (delr * delc) * 1e3, cmap="RdBu_r")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("Kwel/Infiltratie (mm/d)")
fig.savefig(os.path.join(figdir, "river_leakage.png"),
            bbox_inches="tight", dpi=150)

# %%
end = default_timer()
print(f"Elapsed time: {end-start:.1f} s")

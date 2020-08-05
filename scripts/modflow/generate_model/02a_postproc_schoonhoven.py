import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import flopy as fp


model_ws = '../../model/schoonhoven_agg'
model_name = 'schoonhoven_agg'
sim_name = 'mfsim'
figdir = os.path.join(model_ws, 'figure')
headfile = f'{model_name}.hds'
budgetfile = f'{model_name}.cbb'


sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=model_ws)
gwf = sim.get_model(model_name)

delr = gwf.dis.delr.data[0]
delc = gwf.dis.delc.data[0]

# open head file
fname = os.path.join(model_ws, headfile)
hds = fp.utils.HeadFile(fname)
h = hds.get_data(kstpkper=(0, 0))
h[h > 1e5] = np.nan

# open budgetfile
fname = os.path.join(model_ws, budgetfile)
cbc = fp.utils.CellBudgetFile(fname)

spdis = cbc.get_data(text="SPDIS")

q3D = 0.0
t = 0

for txt in ["RIV", "DRN", "GHB"]:
    try:
        q = cbc.get_data(kstpkper=(0, 0), text=txt)[t]
    except Exception:
        print(f"{txt} not in budget file.")
        continue
    q3Di = cbc.create3D(q, gwf.modelgrid.nlay, gwf.modelgrid.nrow,
                        gwf.modelgrid.ncol)
    q3D += q3Di.data

qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(gwf,
                                                            fname,
                                                            precision="double")

# %% plot head and quiver
fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(h[0], cmap="RdBu", vmin=-2, vmax=2)
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

# fig, ax = plt.subplots(1, 1, figsize=(14, 10))
# mapview = fp.plot.PlotMapView(model=gwf)
# qm = mapview.plot_array(h[0], cmap="RdBu", vmin=-2, vmax=2)
# qv = mapview.plot_vector(qx, qy, istep=1, jstep=1, normalize=True,
#                          width=.001, headwidth=4, headlength=5,
#                          headaxislength=4,
#                          pivot="mid")

# divider = make_axes_locatable(ax)
# cax = divider.append_axes('right', size='3%', pad=0.05)
# cbar = plt.colorbar(qm, cax=cax)
# cbar.set_label("head (m NAP)")
# ax.set_title("Grondwaterstand en stroming (genormaliseerd)")
# fig.savefig(os.path.join(figdir, "head_layer0_norm-quiver.png"),
#             bbox_inches="tight", dpi=150)

# %% plot river leakage

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(q3D / (delr * delc) * 1e3, cmap="RdBu_r",
                        vmin=-40, vmax=40)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("Kwel/Infiltratie (mm/d)")
fig.savefig(os.path.join(figdir, "river_leakage.png"),
            bbox_inches="tight", dpi=150)


#%% show riv and shape
# import geopandas as gpd
# datadir = '../../../data'
# water_shp = os.path.join(datadir, "modflow_sfw_schoonhoven", "waterareas.shp")
# sfw = gpd.read_file(water_shp)
# mask = sfw["geometry"].apply(lambda geom: geom.wkb)
# sfw = sfw.loc[mask.drop_duplicates().index]

# fig, ax = plt.subplots(1, 1, figsize=(12, 10))
# mapview = fp.plot.PlotMapView(model=gwf, layer=0)
# mapview.plot_bc("RIV", color="DarkGrey")
# sfw.plot(color="CornFlowerBlue", zorder=10, ax=ax)

# fig.savefig("RIV_vs_shape_l1.png", bbox_inches="tight", dpi=150)
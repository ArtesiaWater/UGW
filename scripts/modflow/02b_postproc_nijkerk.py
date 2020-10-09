import os

import flopy as fp
import hydropandas as hpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Polygon

# mpl.interactive(True)

model_name = 'nkrk_ind'
figtitle_prefix = "Individueel: "
model_ws = f'../../models/{model_name}'
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
h = hds.get_alldata()
h[h > 1e5] = np.nan

h_first_active_layer = h[0, 0].copy()
first_active_layer = np.zeros((gwf.dis.nrow.data, gwf.dis.ncol.data))

counter = 0
while np.sum(np.isnan(h_first_active_layer)) > 0:

    irow, icol = np.where(np.isnan(h_first_active_layer))
    h_first_active_layer[irow, icol] = h[0, counter, irow, icol]
    first_active_layer[irow, icol] = counter
    counter += 1

# add nan data from other layers to first layer in original array
h[0] = h_first_active_layer

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
qm = mapview.plot_array(h[0], cmap="RdBu", vmin=-6, vmax=6)
# qv = mapview.plot_vector(qx, qy, istep=1, jstep=1, normalize=False,
#                          scale=2, alpha=0.75, width=.00175,
#                          headwidth=3, headlength=3, headaxislength=2,
#                          pivot="mid")

plt.yticks(rotation=90, va="center")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("head (m NAP)")
ax.set_title(figtitle_prefix + "grondwaterstand")
ax.set_xlabel("X (m RD)")
ax.set_ylabel("Y (m RD)")

fig.tight_layout()
fig.savefig(os.path.join(figdir, "head_layer0.png"),
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
                        vmin=-30, vmax=30)
plt.yticks(rotation=90, va="center")

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("Kwel/Infiltratie (mm/d)")

ax.set_title(figtitle_prefix + "flux van/naar het oppervlaktewater")
ax.set_xlabel("X (m RD)")
ax.set_ylabel("Y (m RD)")

fig.tight_layout()
fig.savefig(os.path.join(figdir, "river_leakage.png"),
            bbox_inches="tight", dpi=150)

# %% load piezometers
xmin, xmax, ymin, ymax = gwf.modelgrid.extent
model_extent = Polygon(
    [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)])

dinodir = "../../data/piezometers"

files = [f for f in os.listdir(dinodir) if f.startswith("dino_download__")]

collections = []
parsed_files = []

for f in files:
    extent_str = f.split("__")[-1].split(".")[0]
    xmin, ymin, xmax, ymax = np.array(extent_str.split("_"), dtype=np.float)
    p = Polygon([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)])

    if p.intersects(model_extent):
        parsed_files.append(f)
        collections.append(pd.read_pickle(
            os.path.join(dinodir, f), compression="zip"))

pbc = pd.concat(collections, axis=0)

# at least 3 consecutive years of 12 observations
mask = np.any(pbc.stats.consecutive_obs_years() > 3., axis=0)
pbc = pbc.loc[mask]

pbc["modellayer"] = pbc.gwobs.get_modellayers(gwf, verbose=True)

if not pbc["modellayer"].isna().all():
    hobs_model = hpd.ObsCollection.from_modflow(pbc, gwf, h,
                                                [gwf.modeltime.start_datetime])
    hobs_model["hmod"] = hobs_model.stats.mean_in_period(col=0)

pbc["hobs"] = pbc.stats.mean_in_period()

colors = {}
cmap = mpl.cm.get_cmap("tab20")
norm = mpl.colors.Normalize(vmin=-0.5, vmax=19.5)
for i, mlay in enumerate(pbc["modellayer"].dropna().unique()):
    colors[mlay] = cmap(i)


def apply_colors(s):
    if np.isnan(s):
        return (0, 0, 0, 0)
    else:
        return colors[s]


c = pbc["modellayer"].apply(apply_colors)

fig, ax = plt.subplots(1, 1, figsize=(10, 8), dpi=50)
ax.scatter(pbc["hobs"], hobs_model["hmod"], c=np.array(c.values),
           marker="o", edgecolor="k",
           facecolor="C0", s=100)
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = plt.colorbar(sm)
cbar.set_ticks(range(20))
ticks = list(colors.keys())
ticks.sort()
cbar.set_ticklabels(np.array(ticks, dtype=int))
cbar.set_label("Model layer number")
ax.set_aspect("equal", adjustable="box")
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
ax.autoscale(False)
ax.plot([-100, 100], [-100, 100], ls="dashed", color="k", lw=2)
ax.set_xlim(np.min([xmin, ymin]), np.max([xmax, ymax]))
ax.set_ylim(np.min([xmin, ymin]), np.max([xmax, ymax]))
ax.grid(b=True)
ax.set_xlabel("Mean observed head (2015-2019) (m NAP)")
ax.set_ylabel("Modelled head (m NAP)")
fig.tight_layout()

fig.savefig(os.path.join(figdir, "modelled_vs_observed.png"),
            bbox_inches="tight", dpi=150)

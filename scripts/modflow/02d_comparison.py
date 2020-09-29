import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Polygon

import flopy as fp
import hydropandas as hpd

sim_name = 'mfsim'
model_names = ['nkrk_ind', 'nkrk_agg_dl']
model_wss = [f'../../model/{i}' for i in model_names]
figdir = "../../model/comparisons"

heads = []
qlist = []

for i in range(len(model_names)):
    headfile = f'{model_names[i]}.hds'
    budgetfile = f'{model_names[i]}.cbb'

    sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=model_wss[i])
    gwf = sim.get_model(model_names[i])

    delr = gwf.dis.delr.data[0]
    delc = gwf.dis.delc.data[0]

    # open head file
    fname = os.path.join(model_wss[i], headfile)
    hds = fp.utils.HeadFile(fname)
    h = hds.get_alldata()
    # h = hds.get_data(kstpkper=(0, 0))
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

    heads.append(h)

    # open budgetfile
    fname = os.path.join(model_wss[i], budgetfile)
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
        q3Di = cbc.create3D(q,
                            gwf.modelgrid.nlay,
                            gwf.modelgrid.nrow,
                            gwf.modelgrid.ncol)
        q3D += q3Di.data

    qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(
        gwf, fname, precision="double")

    # qlist.append((qx, qy, qz))
    qlist.append(q3D)

# %% Head difference

h1, h2 = heads
dh = h1 - h2

lay = 0
lim = np.nanmax(np.abs(dh[:, lay]))

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(dh[0], cmap="PuOr", vmin=-lim, vmax=lim)
plt.yticks(rotation=90, va="center")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("head difference (m)")
ax.set_title(f"Verschil in stijghoogtes: "
             f"{model_names[0]} - {model_names[1]}")
ax.set_xlabel("X (m RD)")
ax.set_ylabel("Y (m RD)")
fig.tight_layout()
fig.savefig(os.path.join(figdir, f"dhead_{model_names[0]}-{model_names[1]}.png"),
            bbox_inches="tight", dpi=150)

# %% plot river leakage
q3D1, q3D2 = qlist
dq = (q3D1 - q3D2) / (delr * delc) * 1e3

lay = 0
lim = np.nanmax(np.abs(dq[lay]))

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf, layer=lay)
qm = mapview.plot_array(dq, cmap="PuOr", vmin=-lim, vmax=lim)
plt.yticks(rotation=90, va="center")
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("Verschil flux (mm/d)")
ax.set_xlabel("X (m RD)")
ax.set_ylabel("Y (m RD)")
ax.set_title("Verschil in opp. water flux: "
             f"{model_names[0]} - {model_names[1]}")
fig.tight_layout()
fig.savefig(os.path.join(figdir, f"dq_{model_names[0]}-{model_names[1]}.png"),
            bbox_inches="tight", dpi=150)

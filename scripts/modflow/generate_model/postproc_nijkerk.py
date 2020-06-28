import os
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import flopy as fp


model_ws = '../../model/ugw'
model_name = 'ugw'
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

h_first_active_layer = h[0]

counter = 0
while np.sum(np.isnan(h_first_active_layer)) > 0:

    irow, icol = np.where(np.isnan(h_first_active_layer))
    h_first_active_layer[irow, icol] = h[counter, irow, icol]
    counter += 1

# add nan data from other layers to first layer in original array
h[0] = h_first_active_layer

# open budgetfile
fname = os.path.join(model_ws, budgetfile)
cbc = fp.utils.CellBudgetFile(fname)

spdis = cbc.get_data(text="SPDIS")
qriv = cbc.get_data(kstpkper=(0, 0), text="RIV")[0]
qriv3D = cbc.create3D(qriv, gwf.modelgrid.nlay, gwf.modelgrid.nrow,
                      gwf.modelgrid.ncol)
qghb = cbc.get_data(kstpkper=(0, 0), text="GHB")[0]
qghb3D = cbc.create3D(qghb, gwf.modelgrid.nlay, gwf.modelgrid.nrow,
                      gwf.modelgrid.ncol)

q3D = qriv3D.data + qghb3D.data

qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(gwf,
                                                            fname,
                                                            precision="double")

# %% plot head and quiver
fig, ax = plt.subplots(1, 1, figsize=(14, 10))
mapview = fp.plot.PlotMapView(model=gwf)
qm = mapview.plot_array(h[0], cmap="RdBu", vmin=-6, vmax=15)
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
                        vmin=-50, vmax=50)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05)
cbar = plt.colorbar(qm, cax=cax)
cbar.set_label("Kwel/Infiltratie (mm/d)")
fig.savefig(os.path.join(figdir, "river_leakage.png"),
            bbox_inches="tight", dpi=150)

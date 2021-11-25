# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 11:12:00 2020

@author: oebbe
"""

import os

import art_tools
import flopy
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from . import surface_water


def plot_recharge(model_ds, gwf):
    raise NotImplementedError('this does not work yet')
    rch_array = model_ds['recharge'][:, 0]
    fig, ax = plot_array(gwf, rch_array)

    # rch_xarray = xr.zeros_like(model_ds['rch_name'], dtype=float)
    # for i in range(62):
    #     rch_xarray = xr.where(model_ds['rch_name']==f'recharge_{i}', model_ds[f'recharge_{i}'].data[0], rch_xarray)

    # fig, ax = plt.subplots()
    # rch_xarray.plot(ax=ax)

    return fig, ax


def plot_timeseries_result_vertex_grid(hdobj, gwf, model_ds, x, y, lay=0):

    idx = gwf.modelgrid.intersect(x, y)

    # get model heads
    hds_idx = [hdobj.get_data(totim=tim, mflay=lay)[0][idx]
               for tim in hdobj.times]
    hds_df = pd.DataFrame(index=pd.to_datetime(model_ds.time.data),
                          data={f'model head lay {lay}': hds_idx})

    # get maaiveld
    ahn_lvl = model_ds['ahn'].data[idx]

    # get recharge
    rch_idx = model_ds['recharge'].data[idx]
    rch_df = pd.DataFrame(index=pd.to_datetime(model_ds.time.data),
                          data={'recharge': rch_idx})

    # create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # plot maaiveld
    ax.plot_date([hds_df.index[0], hds_df.index[-1]], [ahn_lvl, ahn_lvl], ls=':', marker=None,
                 label='maaiveld')

    # plot model heads
    hds_df.plot(ax=ax, marker='o', lw=0.2)
    ylim1 = ax.get_ylim()

    # plot recharge
    ax2 = ax.twinx()
    xlim = ax2.get_xlim()
    ax2.hlines(0, xlim[0], xlim[1], color='silver')
    ax2.set_xlim(xlim)
    ax2.bar(rch_df.index, rch_df.values[:, 0], label=rch_df.columns[0])

    # add layout stuff
    ax.legend(loc=2)
    ax.grid()
    ax2.legend()
    ylim2 = ax2.get_ylim()
    ax.set_ylim(ylim1[0] - 2 * (ylim1[1] - ylim1[0]), ylim1[1])
    ax.set_ylabel('m NAP')
    ax2.set_ylim(ylim2[0], ylim2[1] + (ylim2[1] - ylim2[0]))
    ax2.set_ylabel('recharge m/dag')

    return ax, ax2


def plot_surface_water(extent):
    # read shapes of water boards (only used in figure for now)
    # ws = geo.get_waterschappen(extent=extent)

    top10_gdf_dic = surface_water.get_top10(extent)

    # plot the model area
    f, ax = plt.subplots(figsize=(8, 7))
    ax.axis('scaled')
    ax.axis(extent)
    plt.yticks(rotation=90, va='center')
    f.tight_layout(pad=0.0)
    art_tools.OpenTopo(ax=ax).plot(alpha=0.5)

    top10_gdf_dic['waterdeel_vlak'].plot(ax=ax)
    top10_gdf_dic['waterdeel_lijn'].plot(ax=ax, linewidth=1)
    #ws.plot(ax=ax, edgecolor='k', facecolor='none')

    return f, ax


def plot_raster(data_array, **kwargs):

    fig, ax = plt.subplots()
    data_array.plot(ax=ax, **kwargs)
    ax.axis('scaled')

    return fig, ax


def plot_array(gwf, array, figsize=(8, 8), colorbar=True, **kwargs):
    f, ax = plt.subplots(figsize=figsize)
    yticklabels = ax.yaxis.get_ticklabels()
    plt.setp(yticklabels, rotation=90, verticalalignment='center')
    ax.axis('scaled')
    pmv = flopy.plot.PlotMapView(modelgrid=gwf.modelgrid, ax=ax)
    pcm = pmv.plot_array(array, **kwargs)
    if colorbar:
        plt.colorbar(pcm)
    f.tight_layout(pad=0.0)
    return f, ax


def plot_map_cross_section(x=None, y=None, shp_fname=None, **kwargs):

    if not (x is None):
        fig, ax = plot_map_cross_section_at_x(x, **kwargs)

    if not (y is None):
        fig, ax = plot_map_cross_section_at_y(y, **kwargs)

    if not (shp_fname is None):
        fig, ax = plot_map_cross_section_at_line_shp(shp_fname, **kwargs)

    return fig, ax


def plot_map_cross_section_at_line_shp(shp_fname, model_ds, gwf, gdf_opp_water, panden_shp,
                                       figdir, kh, ylim=(-100, 30),
                                       cross_section_name='D',
                                       plot_mapview=False):

    shp_doorsnede_a = gpd.read_file(shp_fname)
    line = flopy.plot.plotutil.shapefile_get_vertices(shp_fname)

    # Let's plot the model grid in map view to look at it
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_title(f"Model Grid with cross sectional line {cross_section_name}")

    # use PlotMapView to plot a DISV (vertex) model
    if plot_mapview:
        mapview = flopy.plot.PlotMapView(gwf, layer=0)
        linecollection = mapview.plot_grid(alpha=0.5)
    panden_shp.plot(ax=ax, color='green')

    # plot the line over the model grid
    shp_doorsnede_a.plot(ax=ax, color='red')
    gdf_opp_water.plot(ax=ax, color='blue', alpha=0.8)

    fig.savefig(os.path.join(
        figdir, f'map_view_cross_section_{cross_section_name}.png'))

    # plot cross section
    fig = plt.figure(figsize=(20, 5))
    ax = fig.add_subplot(1, 1, 1)
    xsect = flopy.plot.PlotCrossSection(gwf, line={'line': line[0][:-1]})
    csa = xsect.plot_array(kh, vmin=0.0, vmax=20)
    #linecollection = xsect.plot_grid()
    ax.set_ylim(*ylim)
    t = ax.set_title(
        f"Cross-Section {cross_section_name} with Horizontal hydraulic conductivity")

    cb = plt.colorbar(csa, shrink=0.75)

    fig.savefig(os.path.join(
        figdir, f'cross_section_{cross_section_name}.png'))

    return fig, ax


def plot_map_cross_section_at_y(y, model_ds, gwf, gdf_opp_water,
                                panden_shp,
                                figdir, kh, ylim=(-100, 30),
                                cross_section_name='D',
                                plot_mapview=False,
                                gdf=None):

    row = gwf.modelgrid.intersect(model_ds.extent[1].values - 1, y)[0]

    line_vert = np.array([(model_ds.extent[0], y),
                          (model_ds.extent[1], y)])

    # Let's plot the model grid in map view to look at it
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_title("Model Grid with cross sectional line")

    # use PlotMapView to plot a DISV (vertex) model
    if plot_mapview:
        mapview = flopy.plot.PlotMapView(gwf, layer=0)
        linecollection = mapview.plot_grid()
    panden_shp.plot(ax=ax, color='green', alpha=0.8)

    # plot the line over the model grid
    lc = plt.plot(line_vert.T[0], line_vert.T[1], 'r--', lw=2)
    gdf_opp_water.plot(ax=ax, color='blue', alpha=0.5)
    if gdf is not None:
        gdf.plot(ax=ax, facecolor="none", edgecolor='black', lw=3)

    fig.savefig(os.path.join(
        figdir, f'map_view_cross_section_{cross_section_name}.png'))

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(1, 1, 1)

    xsect = flopy.plot.PlotCrossSection(model=gwf, line={'Row': row})
    csa = xsect.plot_array(kh, vmin=0.0, vmax=20)
    linecollection = xsect.plot_grid()
    ax.set_ylim(*ylim)
    t = ax.set_title(
        f"y: {y} Cross-Section {cross_section_name} with Horizontal hydraulic conductivity")
    cb = plt.colorbar(csa, shrink=0.75)

    fig.savefig(os.path.join(
        figdir, f'cross_section_{cross_section_name}.png'))

    return fig, ax


def plot_map_cross_section_at_x(x, model_ds, gwf, gdf_opp_water, panden_shp,
                                figdir, kh, ylim=(-100, 30), cross_section_name='A',
                                plot_mapview=False, gdf=None):

    print(f'creating figure for cross section -> {cross_section_name}')

    column = gwf.modelgrid.intersect(x + 1, model_ds.extent[2].values)[1]

    line_vert = np.array([(x, model_ds.extent[2]),
                          (x, model_ds.extent[3])])

    # Let's plot the model grid in map view to look at it
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_title("Model Grid with cross sectional line")

    # use PlotMapView to plot a DISV (vertex) model
    if plot_mapview:
        mapview = flopy.plot.PlotMapView(gwf, layer=0)
        linecollection = mapview.plot_grid()

    # plot the line over the model grid
    lc = plt.plot(line_vert.T[0], line_vert.T[1], 'r--', lw=2)
    gdf_opp_water.plot(ax=ax, color='blue', alpha=0.5)
    panden_shp.plot(ax=ax, color='green', alpha=0.5)
    if gdf is not None:
        gdf.plot(ax=ax, facecolor="none", edgecolor='black', lw=3)

    fig.savefig(os.path.join(figdir,
                             f'map_view_cross_section_{cross_section_name}.png'))

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(1, 1, 1)
    # plot the horizontal hydraulic conductivities

    xsect = flopy.plot.PlotCrossSection(model=gwf, line={'Column': column})
    csa = xsect.plot_array(kh, vmin=0.0, vmax=20)
    linecollection = xsect.plot_grid()
    ax.set_ylim(*ylim)
    t = ax.set_title(
        f"x: {x} Cross-Section {cross_section_name} with Horizontal hydraulic conductivity")
    cb = plt.colorbar(csa, shrink=0.75)

    # flip x axis
    #cur_xlim = ax.get_xlim()
    #ax.set_xlim(cur_xlim[1], 10000)

    fig.savefig(os.path.join(
        figdir, f'cross_section_{cross_section_name}.png'))

    return fig, ax


def plot_layer_stitches_y(y, model_ds, gwf, gdf_opp_water,
                          panden_shp, figdir, lay_no,
                          cross_section_name='',
                          ylim=(-400, 30),
                          plot_mapview=True,
                          gdf_modgeb_pwn=None):

    row = gwf.modelgrid.intersect(model_ds.extent[1].values - 1, y)[0]
    line_vert = np.array([(model_ds.extent[0], y),
                          (model_ds.extent[1], y)])

    # Let's plot the model grid in map view to look at it
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_title("Model Grid with cross sectional line")

    # use PlotMapView to plot a DISV (vertex) model
    mapview = flopy.plot.PlotMapView(gwf, layer=0)
    linecollection = mapview.plot_grid()

    # plot the line over the model grid
    lc = plt.plot(line_vert.T[0], line_vert.T[1], 'r--', lw=2)
    gdf_opp_water.plot(ax=ax, color='blue', alpha=0.5)
    panden_shp.plot(ax=ax, color='green', alpha=0.5)

    gdf_modgeb_pwn.plot(ax=ax, facecolor="none", edgecolor='black', lw=3)

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(1, 1, 1)

    xsect = flopy.plot.PlotCrossSection(model=gwf, line={'Row': row})
    csa = xsect.plot_array(lay_no, cmap='tab20', alpha=0.5)
    #linecollection = xsect.plot_grid()

    # plot laagscheidingen
    cmap = mpl.cm.get_cmap('tab20')

    for i, lay in enumerate(model_ds.layer):
        lay_bot = model_ds['bot'].sel(layer=lay)[row, :]
        lay_bot_dis = lay_bot.x.values - lay_bot.x.values.min()
        color_num = lay_no[i][0][0]
        ax.plot(lay_bot_dis, lay_bot.values, color=cmap(color_num / 10),
                lw=2, ls=':')

    # voeg verticale lijnen toe op grens van de modellen
    x_bound_pwn = gdf_modgeb_pwn.bounds.loc[0, 'maxx']
    x_bound_res = x_bound_pwn - model_ds.extent[0]
    ax.vlines(x_bound_res, ylim[0], ylim[1], ls=':', color='red')
    ax.set_ylim(*ylim)

    t = ax.set_title(
        f"y: {y} Cross-Section {cross_section_name} with layer numbers")
    cb = plt.colorbar(csa, shrink=0.75)

    fig.savefig(os.path.join(figdir, f'stitches_{cross_section_name}.png'),
                bbox_inches='tight', dpi=300)

    return fig, ax


def plot_layer_stitches_x(x, model_ds, gwf, gdf_opp_water,
                          panden_shp, figdir, lay_no,
                          cross_section_name='',
                          ylim=(-400, 30),
                          plot_mapview=True,
                          gdf_modgeb_pwn=None):

    column = gwf.modelgrid.intersect(x + 1, model_ds.extent[2].values)[1]
    line_vert = np.array([(x, model_ds.extent[2]),
                          (x, model_ds.extent[3])])

    # Let's plot the model grid in map view to look at it
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_title("Model Grid with cross sectional line")

    # use PlotMapView to plot a DISV (vertex) model
    mapview = flopy.plot.PlotMapView(gwf, layer=0)
    linecollection = mapview.plot_grid()

    # plot the line over the model grid
    lc = plt.plot(line_vert.T[0], line_vert.T[1], 'r--', lw=2)
    gdf_opp_water.plot(ax=ax, color='blue', alpha=0.5)
    panden_shp.plot(ax=ax, color='green', alpha=0.5)

    gdf_modgeb_pwn.plot(ax=ax, facecolor="none", edgecolor='black', lw=3)

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(1, 1, 1)

    xsect = flopy.plot.PlotCrossSection(model=gwf, line={'Column': column})
    csa = xsect.plot_array(lay_no, cmap='tab20', alpha=0.5)
    #linecollection = xsect.plot_grid()

    # plot laagscheidingen
    cmap = mpl.cm.get_cmap('tab20')

    for i, lay in enumerate(model_ds.layer):
        lay_bot = model_ds['bot'].sel(layer=lay)[:, column]
        lay_bot_dis = (lay_bot.y.values - lay_bot.y.values.min())[::-1]
        color_num = lay_no[i][0][0]
        ax.plot(lay_bot_dis, lay_bot.values, color=cmap(color_num / 10),
                lw=2, ls=':')

    # voeg verticale lijnen toe op grens van de modellen
    if x < gdf_modgeb_pwn.bounds.loc[0, 'maxx']:
        y_min_bound_pwn = gdf_modgeb_pwn.bounds.loc[0, 'miny']
        y_min_bound_res = model_ds.extent[3] - y_min_bound_pwn

        y_max_bound_pwn = gdf_modgeb_pwn.bounds.loc[0, 'maxy']
        y_max_bound_res = model_ds.extent[3] - y_max_bound_pwn

        ax.vlines(y_min_bound_res, ylim[0], ylim[1], ls=':', color='red')
        ax.vlines(y_max_bound_res, ylim[0], ylim[1], ls=':', color='red')

    ax.set_ylim(*ylim)

    t = ax.set_title(
        f"x: {x} Cross-Section {cross_section_name} with layer numbers")
    cb = plt.colorbar(csa, shrink=0.75)

    fig.savefig(os.path.join(figdir, f'stitches_{cross_section_name}.png'),
                bbox_inches='tight', dpi=300)

    return fig, ax


def compare_layer_models_cross_section(y, lay_mod1, lay_mod2, lay_mod1_bot_name='bottom',
                                       layer_mod1=0, layer_mod2=0,
                                       lay_mod2_bot_name='bottom', ylim=(-200, 20)):
    #y = 502000.
    row_lay_mod1 = np.where(lay_mod1.y == (
        lay_mod1.y.sel(y=y, method='nearest').data))[0]
    row_lay_mod2 = np.where(lay_mod2.y == (
        lay_mod2.y.sel(y=y, method='nearest').data))[0]

    fig = plt.figure(figsize=(10, 10))

    z1 = lay_mod1[lay_mod1_bot_name][:, row_lay_mod1]
    plt.plot(lay_mod1['top'][0, row_lay_mod1].x.data,
             lay_mod1['top'][0, row_lay_mod1][0], c="k")
    for zi in z1:
        #plt.plot(range(z1.shape[1]), zi, c="k")
        plt.plot(zi.x.data, zi[0], c="k")
    for i in range(len(z1) - 1):
        if i == layer_mod1:
            c = "lightblue"
            plt.fill_between(
                z1.x.data, z1[i][0], lay_mod1['top'][i, row_lay_mod1][0], alpha=0.5, color=c)
        # else:

        # c = "lightblue"

        # plt.fill_between(z1.x.data, z1[i][0], z1[i+1][0], alpha=0.5, color=c)

    z2 = lay_mod2[lay_mod2_bot_name][:, row_lay_mod2]
    plt.plot(lay_mod2['top'][0, row_lay_mod2].x.data,
             lay_mod2['top'][0, row_lay_mod2][0], c="r")
    for zi in z2:
        plt.plot(zi.x.data, zi[0], c="r", ls="dashed")
    for i in range(len(z2) - 1):
        if i == layer_mod2:
            c = "lightcoral"
            plt.fill_between(
                z2.x.data, z2[i][0], lay_mod2['top'][i, row_lay_mod2][0], alpha=0.5, color=c)

    plt.ylim(*ylim)
    plt.grid()
    plt.ylabel("elevation")

    return fig


def plot_comparison_modellayers_topview(lay_mod1, lay_mod2,
                                        layer_mod1=0, layer_mod2=0,
                                        lay_mod1_bot_name='bottom',
                                        lay_mod2_bot_name='bottom',
                                        xsel=None, ysel=None):

    compare = compare_layer_models_top_view(lay_mod1, lay_mod2,
                                            layer_mod1, layer_mod2,
                                            lay_mod1_bot_name,
                                            lay_mod2_bot_name,
                                            xsel, ysel)

    cmap = mpl.cm.Paired
    norm = mpl.colors.Normalize(vmin=0.5, vmax=12.5)

    extent = [compare.x[0], compare.x[-1], compare.y[-1], compare.y[0]]
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    im = ax.matshow(compare, cmap=cmap, norm=norm,
                    extent=extent, aspect='auto')
    cbar = plt.colorbar(im, shrink=0.75)
    names = [
        "equal",
        "top inside",
        "bot inside",
        "inside",
        "outside",
        "top outside",
        "bot outside",
        "below",
        "shifted down",
        "shifted up",
        "above",
        "nan",
    ]
    cbar.set_ticks(range(1, 13))
    cbar.ax.set_yticklabels(names)

    return fig, ax, compare


def compare_layer_models_top_view(lay_mod1, lay_mod2,
                                  layer_mod1=0, layer_mod2=0,
                                  lay_mod1_bot_name='bottom',
                                  lay_mod2_bot_name='bottom',
                                  xsel=None, ysel=None):

    if ysel is None:
        ysel = lay_mod1.y.data
    if xsel is None:
        xsel = lay_mod1.x.data

    lay_mod1 = lay_mod1.sel(x=xsel, y=ysel)
    lay_mod2 = lay_mod2.sel(x=xsel, y=ysel)

    bot1 = lay_mod1[lay_mod1_bot_name][layer_mod1]
    if len(lay_mod1['top'].dims) == 2:
        if layer_mod1 == 0:
            top1 = lay_mod1['top']
        else:
            top1 = lay_mod1[lay_mod1_bot_name][layer_mod1 - 1]
    elif len(lay_mod1['top'].dims) == 3:
        top1 = lay_mod1['top'][layer_mod1]

    bot2 = lay_mod2[lay_mod2_bot_name][layer_mod2]
    if len(lay_mod2['top'].dims) == 2:
        if layer_mod2 == 0:
            top2 = lay_mod2['top']
        else:
            top2 = lay_mod2[lay_mod2_bot_name][layer_mod2 - 1]
        top2 = lay_mod2['top']
    elif len(lay_mod2['top'].dims) == 3:
        top2 = lay_mod2['top'][layer_mod2]

    compare = compare_top_bots(bot1, top1, bot2, top2)

    return compare


def compare_layer_models_border(lay_mod1, lay_mod2, border_mask,
                                layer_mod1=0, layer_mod2=0,
                                lay_mod1_bot_name='bottom',
                                lay_mod2_bot_name='bottom',
                                xsel=None, ysel=None):

    if ysel is None:
        ysel = lay_mod1.y.data
    if xsel is None:
        xsel = lay_mod1.x.data

    lay_mod1 = lay_mod1.sel(x=xsel, y=ysel)
    lay_mod2 = lay_mod2.sel(x=xsel, y=ysel)

    bot1 = lay_mod1[lay_mod1_bot_name][layer_mod1]
    if len(lay_mod1['top'].dims) == 2:
        if layer_mod1 == 0:
            top1 = lay_mod1['top']
        else:
            top1 = lay_mod1[lay_mod1_bot_name][layer_mod1 - 1]
    elif len(lay_mod1['top'].dims) == 3:
        top1 = lay_mod1['top'][layer_mod1]

    bot2 = lay_mod2[lay_mod2_bot_name][layer_mod2]
    if len(lay_mod2['top'].dims) == 2:
        if layer_mod2 == 0:
            top2 = lay_mod2['top']
        else:
            top2 = lay_mod2[lay_mod2_bot_name][layer_mod2 - 1]
        top2 = lay_mod2['top']
    elif len(lay_mod2['top'].dims) == 3:
        top2 = lay_mod2['top'][layer_mod2]

    bot1 = xr.where(border_mask, bot1, np.nan)
    bot2 = xr.where(border_mask, bot2, np.nan)
    top1 = xr.where(border_mask, top1, np.nan)
    top2 = xr.where(border_mask, top2, np.nan)

    compare = compare_top_bots(bot1, top1, bot2, top2)

    return compare


def compare_top_bots(bot1, top1, bot2, top2):

    compare = xr.zeros_like(top1)

    # 1: equal
    mask_eq = (top1 == top2) & (bot1 == bot2)

    # 2: top within: top2 lower & bot2 equal
    mask_top_within = (top1 > top2) & (bot1 == bot2)

    # 3: bottom within: bot2 higher & top2 equal
    mask_bot_within = (top1 == top2) & (bot1 < bot2)

    # 4: within: top2 lower & bot2 higher
    mask_within = (top1 > top2) & (bot1 < bot2)

    # 5: outside: top2 higher & bot2 lower
    mask_outside = (top1 < top2) & (bot1 > bot2)

    # 6: top outside: top2 higher & bot2 equal
    mask_top_oustide = (top1 < top2) & (bot1 == bot2)

    # 7: bot outside: bot2 lower & top2 equal
    mask_bot_outside = (top1 == top2) & (bot1 > bot2)

    # 8: under: bot1 >= top2
    mask_under = (bot1 >= top2)

    # 9: shifted down: (top1 > top2 > bot1) & (bot1 > bot2)
    mask_shift_down = ((top1 > top2) & (top2 > bot1)) & (bot1 > bot2)

    # 10: shifted up: (top1 < top2) & (bot1 < bot2 < top1)
    mask_shift_up = (top1 < top2) & ((bot1 < bot2) & (bot2 < top1))

    # 11: above: top1 <= bot2
    mask_above = (top1 <= bot2)

    # 12: bot1 is nan
    mask_botnan = np.logical_or(np.isnan(bot1), np.isnan(bot2))
    mask_topnan = np.logical_or(np.isnan(top1), np.isnan(top2))
    mask_nan = np.logical_or(mask_botnan, mask_topnan)

    compare = xr.where(mask_eq, 1, compare)
    compare = xr.where(mask_top_within, 2, compare)
    compare = xr.where(mask_bot_within, 3, compare)
    compare = xr.where(mask_within, 4, compare)
    compare = xr.where(mask_outside, 5, compare)
    compare = xr.where(mask_top_oustide, 6, compare)
    compare = xr.where(mask_bot_outside, 7, compare)
    compare = xr.where(mask_under, 8, compare)
    compare = xr.where(mask_shift_down, 9, compare)
    compare = xr.where(mask_shift_up, 10, compare)
    compare = xr.where(mask_above, 11, compare)
    compare = xr.where(mask_nan, 12, compare)

    return compare

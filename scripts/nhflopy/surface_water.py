# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 11:18:22 2020

@author: oebbe
"""
import json
import os
import time
import warnings

import fiona
import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.features
import requests
import xarray as xr
from osgeo import gdal
from shapely.geometry import LineString, Point, Polygon, mapping, shape
from shapely.strtree import STRtree
from tqdm import tqdm

from . import mgrid, rws, util


def get_general_head_boundary(model_ds, gdf,
                              modelgrid, name,
                              gridtype='structured',
                              cachedir=None,
                              use_cache=False,
                              verbose=False):
    """ Get general head boundary from surface water geodataframe

    Parameters
    ----------
    model_ds : xr.DataSet
        dataset containing relevant model grid information
    gdf : geopandas.GeoDataFrame
        polygon shapes with surface water.
    modelgrid : flopy grid
        model grid.
    name : str
        name of the polygon shapes, name is used to store data arrays in 
        model_ds
    gridtype : str, optional
        type of grid, options are 'structured' and 'unstructured'. The default is 'structured'.
    cachedir : str, optional
        directory to store cached values, if None a temporary directory is
        used. default is None
    use_cache : bool, optional
        if True the cached ghb data is used. The default is False.
    verbose : bool, optional
        print additional information to the screen. The default is False.

    Returns
    -------
    model_ds : xr.DataSet
        dataset with spatial model data including the ghb rasters

    """
    model_ds = util.get_cache_netcdf(use_cache, cachedir, 'ghb_model_ds.nc',
                                     gdf_to_model_dataset,
                                     model_ds, verbose=verbose, gdf=gdf,
                                     modelgrid=modelgrid, name=name,
                                     gridtype=gridtype)

    return model_ds


def gdf_to_model_dataset(model_ds, gdf, modelgrid, name, gridtype='structured'):
    """ create 3 data-arrays from a geodataframe with oppervlaktewater:
    - area: with the area of the geodataframe in the cell
    - cond: with the conductance based on the area and bweerstand column in gdf
    - peil: with the surface water lvl based on the peil column in the gdf


    Parameters
    ----------
    model_ds : xr.DataSet
        xarray with model data
    gdf : geopandas.GeoDataFrame
        polygon shapes with surface water.
    modelgrid : flopy grid
        model grid.
    name : str
        name of the polygon shapes, name is used to store data arrays in 
        model_ds

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with modelgrid data. Has 

    """
    area = xr.zeros_like(model_ds['top'])
    cond = xr.zeros_like(model_ds['top'])
    peil = xr.zeros_like(model_ds['top'])
    for i, row in gdf.iterrows():
        area_pol = mgrid.polygon_to_area(modelgrid, row['geometry'],
                                         xr.ones_like(model_ds['top']),
                                         gridtype)
        cond = xr.where(area_pol > area, area_pol / row['bweerstand'], cond)
        peil = xr.where(area_pol > area, row['peil'], peil)
        area = xr.where(area_pol > area, area_pol, area)

    model_ds_out = util.get_model_ds_empty(model_ds)
    model_ds_out[f'{name}_area'] = area
    model_ds_out[f'{name}_cond'] = cond
    model_ds_out[f'{name}_peil'] = peil

    return model_ds_out


def get_modelgrid_sea(gdf_sea, mfgrid, model_ds,
                      gridtype='structured',
                      cachedir=None, use_cache=False,
                      verbose=False):
    """ Get DataArray which is 1 at sea and 0 overywhere else.
    Sea is defined by the geometries in gdf_sea
    grid is defined by mfgrid and model_ds

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        polygon shapes with surface water.
    modelgrid : flopy grid
        model grid.
    model_ds : xr.DataSet
        xarray with model data
    cachedir : str, optional
        directory to store cached values, if None a temporary directory is
        used. default is None
    use_cache : bool, optional
        if True the cached sea data is used. The default is False.
    verbose : bool, optional
        print additional information to the screen. The default is False.

    Returns
    -------
    model_ds : xr.DataSet
        dataset with 'sea' DataVariable.

    """
    model_ds = util.get_cache_netcdf(use_cache, cachedir, 'sea_model_ds.nc',
                                     mgrid.gdf_to_bool_dataset, model_ds,
                                     mfgrid=mfgrid,
                                     gdf=gdf_sea, da_name='sea',
                                     check_time=False,
                                     verbose=verbose)

    return model_ds


def cut_gdf_by_polygon_gdf(gdf_ow, gdf_pg, geom_type='LineString'):
    add_columns = gdf_pg.columns[~gdf_pg.columns.isin(gdf_ow.columns)]
    # set the index as as attribute to the(multi)-polygons
    geometries = gdf_pg['geometry']
    for ind in geometries.index:
        geometries[ind].ind = ind
    s = STRtree(geometries)
    if geom_type == 'LineString':
        ignore = ['Point', 'MultiPoint']
    else:
        ignore = ['Point', 'MultiPoint', 'LineString', 'MultiLineString']

    def add_row_to_list(i, row, props, rows, geom_type, ignore):
        # this code is in a method, so it can call itself for multi-geometries
        if i.geom_type == geom_type:
            rown = row.copy()
            # set the geoemtry to the intersection
            rown['geometry'] = i
            # add the columns from gdf_pg
            rown = rown.append(props)
            rown.name = row.name
            rows.append(rown)
        elif i.geom_type in ['GeometryCollection', 'Multi{}'.format(geom_type)]:
            # call this method on each of the children
            for ii in i:
                add_row_to_list(ii, row, props, rows, geom_type, ignore)
        elif i.geom_type in ignore:
            # the intersection consists of only one or multiple points
            pass
        else:
            raise(Exception('{} not yet supported'.format(i.geom_type)))
    rows = []
    rows_out = []
    for index, row in tqdm(gdf_ow.iterrows(), total=gdf_ow.shape[0],
                           desc="cut gdf by polygons"):
        row_out = row.copy()
        for r in s.query(row.geometry):
            if not r.is_valid:
                r = r.buffer(0)
            i = row.geometry.intersection(r)

            if i.is_empty:
                continue
            props = gdf_pg.loc[r.ind, add_columns]
            add_row_to_list(i, row, props, rows, geom_type, ignore)

            # store all geometries that are outside any of the polygons
            if not row_out.geometry.is_empty:
                row_out.geometry = row_out.geometry.difference(r)
        if not row_out.geometry.is_empty:
            rows_out.append(row_out)

    # make a GeoDataFrame from all the geometries inside gdf_pg
    gdf_hhnk = gpd.GeoDataFrame(rows)

    # add all watercourses that are outside the gdf_pg
    gdf_ow_out = gpd.GeoDataFrame(rows_out)
    gdf_hhnk = gdf_hhnk.append(gdf_ow_out)

    return gdf_hhnk


def remove_intersecting_features(bgt, other):
    """Remove features in bgt that intersect with another geodataframe"""
    keep = pd.Series(True, bgt.index)
    for polygon in other.geometry:
        mask = bgt.intersects(polygon)
        keep[mask] = False
    return bgt[keep]


def calculate_min_ahn_in_polygons(bgt, fname, verbose=True):
    """Get the minimum ahn-depth in the polygons of the BGT"""
    # da = xr.open_rasterio(fname)
    # da = da.isel(band=0).drop_vars('band')

    try:
        from rasterstats import zonal_stats
        stats = zonal_stats(bgt['geometry'], fname, stats='min',
                            all_touched=True)
        return [x['min'] for x in stats]
    except:
        print('Rasterstats failed. Using slower method.')
    ds = rasterio.open(fname)
    shapes = list(zip(bgt['geometry'], range(bgt.shape[0])))
    index = rasterio.features.rasterize(shapes, ds.shape, fill=np.NaN,
                                        transform=ds.transform)
    values = ds.read()[0]
    values[values == ds.nodata] = np.NaN
    ahn_min = []
    for i in tqdm(range(bgt.shape[0]), desc="Min AHN"):
        mask = index == i
        if mask.any():
            ahn_min.append(np.nanmin(values[mask]))
        else:
            ahn_min.append(np.NaN)
    return ahn_min


def gdf2grid(gdf, ml, method="vertex", **kwargs):
    """
    Cut a geodataframe gdf by the grid of a flopy modflow model ml. This method
    is just a wrapper around the GridIntersect method from flopy

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        A GeoDataFrame that needs to be cut by the grid. The GeoDataFrame can
        consist of multiple types (Point, LineString, Polygon and the Multi-
        variants).
    ml : flopy.modflow.Modflow or flopy.mf6.ModflowGwf
        The flopy model that defines the grid.
    method : string, optional
        Method passed to the GridIntersect-class. The default is None, which
        makes GridIntersect choose the best method.
    **kwargs : keyword arguments
        keyword arguments are passed to the intersect_*-methods.

    Returns
    -------
    geopandas.GeoDataFrame
        The GeoDataFrame with the geometries per grid-cell.

    """
    ix = flopy.utils.GridIntersect(ml.modelgrid, method=method)
    shps = []
    for _, shp in tqdm(gdf.iterrows(), total=gdf.shape[0],
                       desc="Intersecting with grid"):

        r = ix.intersect(shp.geometry, **kwargs)

        for i in range(r.shape[0]):
            shpn = shp.copy()
            shpn['cellid'] = r['cellids'][i]
            shpn.geometry = r['ixshapes'][i]
            shps.append(shpn)
    return gpd.GeoDataFrame(shps)


def extent2polygon(extent):
    """Make a Polygon of the extent of a matplotlib axes"""
    nw = (extent[0], extent[2])
    no = (extent[1], extent[2])
    zo = (extent[1], extent[3])
    zw = (extent[0], extent[3])
    polygon = Polygon([nw, no, zo, zw])
    return polygon


def get_bgt(extent, layer="waterdeel", cut_by_extent=True):
    """
    Get geometries within an extent or polygon from the Basis Registratie Grootschalige Topografie (BGT)

    Parameters
    ----------
    extent : list or tuple of length 4 or shapely Polygon
        The extent (xmin, xmax, ymin, ymax) or polygon for which shapes are
        requested.
    layer : string, optional
        The layer for which shapes are requested. The default is "waterdeel".
    cut_by_extent : bool, optional
        Only return the intersection with the extent if True. The default is True

    Returns
    -------
    gdf : GeoPandas GeoDataFrame
        A GeoDataFrame containing all geometries and properties.

    """

    api_url = 'https://api.pdok.nl'
    url = '{}/lv/bgt/download/v1_0/full/custom'.format(api_url)
    body = {"format": "citygml",
            "featuretypes": [layer]}

    if isinstance(extent, Polygon):
        polygon = extent
    else:
        polygon = extent2polygon(extent)

    body['geofilter'] = polygon.to_wkt()

    headers = {'content-type': 'application/json'}

    response = requests.post(url, headers=headers, data=json.dumps(body))

    # check api-status, if completed, download
    if response.status_code in range(200, 300):
        running = True
        href = response.json()["_links"]["status"]["href"]
        url = '{}{}'.format(api_url, href)

        while running:
            response = requests.get(url)
            if response.status_code in range(200, 300):
                status = response.json()['status']
                if status == "COMPLETED":
                    running = False
                else:
                    time.sleep(2)
            else:
                running = False
    else:
        msg = 'Download of bgt-data failed: {}'.format(response.text)
        raise(Exception(msg))

    href = response.json()["_links"]["download"]["href"]
    response = requests.get('{}{}'.format(api_url, href))

    vsif = fiona.ogrext.buffer_to_virtual_file(bytes(response.content))
    vsiz = vsif + '.zip'
    gdal.Rename(vsif, vsiz)

    fc = fiona.Collection(vsiz, vsi='zip')
    gdf = gpd.GeoDataFrame.from_features(
        [feature for feature in fc], crs='epsg:28992')

    # remove double features by removing features with an eindRegistratie
    gdf = gdf[gdf['eindRegistratie'].isna()]

    # re-order columns
    columns = [col for col in gdf.columns if not col ==
               'geometry'] + ['geometry']
    gdf = gdf[columns]

    if cut_by_extent:
        gdf.geometry = gdf.intersection(polygon)
        gdf = gdf[~gdf.is_empty]

    return gdf


def get_brt(extent, layer='waterdelen', crs="epsg:28992", pagesize=1000,
            key='ca5f0a1b-64ae-49ed-83fd-81d4ee128aad',
            return_raw_response=False, only_keep_first_name=True):
    """
    Get geometries within an extent or polygon from the Basis Registratie Topografie (BRT)
    https://brt.basisregistraties.overheid.nl/
    https://www.pdok.nl/introductie/-/article/basisregistratie-topografie-brt-topnl

    Parameters
    ----------
    extent : list or tuple of length 4 or shapely Polygon
        The extent (xmin, xmax, ymin, ymax) or polygon for which shapes are
        requested.
    layer : string, optional
        The layer for which shapes are requested. The default is 'waterdelen'.
    crs : string, optional
        The coordinate reference system. The default is "epsg:28992".
    pagesize : int, optional
        The number of features that is requsted per page. The default is 1000.
    key : string, optional
        The BGT api-key. The default is 'ca5f0a1b-64ae-49ed-83fd-81d4ee128aad'.
    return_raw_response : bool, optional
        Return the raw request-response of the first page when True. The
        default is False.
    only_keep_first_name : bool, optional
        The attribute 'naamNL' (or 'naamFries') can contain multiple names,
        which makes processing more difficult. Only keep the first name of a
        geometry when only_keep_first_name is True. The default is True.

    Returns
    -------
    gdf : GeoPandas GeoDataFrame or requests.models.Response
        A GeoDataFrame containing all geometries and properties.

    """
    #url = 'http://brt.basisregistraties.overheid.nl/def/top10nl'
    url = 'https://brt.basisregistraties.overheid.nl/api/v2/'
    url = '{}{}'.format(url, layer)

    if isinstance(extent, Polygon):
        polygon = extent
    else:
        polygon = extent2polygon(extent)

    headers = {'Accept-Crs': crs,
               'Content-Crs': crs,
               'content-type': 'application/json',
               'X-Api-Key': key,
               'pageSize': str(pagesize)}
    body = {'geometrie': {'intersects': mapping(polygon)}}
    response = requests.post(url, headers=headers, data=json.dumps(body))

    if response.status_code in range(200, 300):
        pass
    else:
        msg = 'Download of brt-data failed: {}'.format(response.text)
        raise(Exception(msg))

    if return_raw_response:
        return response

    df = pd.DataFrame(response.json()['_embedded'][layer])
    page = 1
    while 'next' in response.json()['_links']:
        page = page + 1
        headers['page'] = str(page)

        response = requests.post(url, headers=headers, data=json.dumps(body))

        if response.status_code in range(200, 300):
            pass
        else:
            msg = 'Download of brt-data failed: {}'.format(response.text)
            raise(Exception(msg))
        df = df.append(pd.DataFrame(response.json()['_embedded'][layer]))
    df['link'] = [links['self']['href'] for links in df['_links']]
    geometry = [shape(geojson['geometrie']) for geojson in df['_embedded']]
    df = df.drop(columns=['_embedded', '_links'])
    gdf = gpd.GeoDataFrame(df, geometry=geometry)
    if only_keep_first_name and 'naamNL' in gdf.columns:
        gdf['naamNL'] = [x[0] if len(x) > 0 else '' for x in gdf['naamNL']]
    if only_keep_first_name and 'naamFries' in gdf.columns:
        gdf['naamFries'] = [x[0] if len(
            x) > 0 else '' for x in gdf['naamFries']]
    return gdf


def get_brt_in_parts(extent, dx=1000, dy=None, column='link', **kwargs):
    if dy is None:
        dy = dx
    x = np.arange(extent[0], extent[1], dx)
    y = np.arange(extent[2], extent[3], dy)
    with tqdm(total=len(x) * len(y)) as pbar:
        for x0 in x:
            for y0 in y:
                pbar.update()
                ext = [x0, min(x0 + dx, extent[1]), y0,
                       min(y0 + dx, extent[3])]
                # print(ext)
                gdft = get_brt(ext, **kwargs)
                if x0 == extent[0] and y0 == extent[2]:
                    gdf = gdft
                else:
                    gdf = gdf.append(gdft)
    _, iloc = np.unique(gdf[column], return_index=True)
    gdf = gdf.iloc[iloc]
    return gdf


def read_top10nl_from_zipfile(zipfile, path, files, layer='waterdeel',
                              extent=None):
    """
    Read top10nl (brt) shapes from a zipfile that can be downloaden from
    https://www.pdok.nl/downloads/-/article/basisregistratie-topografie-brt-topnl

    Parameters
    ----------
    zipfile : string
        THe zipfile containing the top10nl files.
    path : string
        THe path within the zipfile that contains the gml-tiles.
    files : list
        THe files within path within the zipfile.
    layer : string, optional
        THe layer for which to read the shapes. The default is 'waterdeel'.

    Returns
    -------
    gdft : TYPE
        DESCRIPTION.

    Example
    -------
    from pathlib import Path
    zipfile = os.path.join(Path.home(),'Downloads','top10nl.zip')
    path = r'TOP10NL_GML_50d_Blokken_april_2020/TOP10NL_GML_50d_Blokken'
    files = ['TOP10NL_25W.gml','TOP10NL_25O.gml']
    gdf = read_top10nl_from_zipfile(zipfile, path, files)

    """

    if extent is not None:
        bbox = (extent[0], extent[2], extent[1], extent[3])
    else:
        bbox = None
    if not isinstance(files, list):
        files = [files]
    for i, file in tqdm(enumerate(files), total=len(files),
                        desc='Reading gml-files'):
        filename = "zip:{}!{}".format(zipfile, os.path.join(path, file))
        filename = filename.replace('\\', '/')
        gdft = gpd.read_file(filename, layer=layer, bbox=bbox)
        if i == 0:
            gdf = gdft
        else:
            gdf = gdf.append(gdft)
    gdf = gdf.set_index('lokaalID')
    if len(files) > 1:
        # drop duplicate shapes, caused by instances in multiple tiles
        gdf[~gdf.index.duplicated()]
    if extent is not None:
        polygon = extent2polygon(extent)
        gdf = gdf[gdf.intersects(polygon)]
        gdf['geometry'] = gdf.intersection(polygon)
    return gdf


def request_waterinfo_waterlevels(riv2stn, email_adress, tmin=None, tmax=None,
                                  check_available_measurements=True):
    # request measurements by mail
    names = np.unique(np.hstack([riv2stn[key] for key in riv2stn]))

    locations = rws.get_locations()
    locations = locations[locations['Naam'].isin(names)]
    if check_available_measurements:
        keep = pd.Series(False, index=locations.index)
        for ind in tqdm(locations.index, desc='Checking if there are measurements'):
            if rws.check_waarnemingen_aanwezig(locations.loc[ind], tmin=tmin,
                                               tmax=tmax):
                keep.at[ind] = True
        locations = locations[keep]

        # check if all names are in locations:
        mask = ~pd.Index(names).isin(locations['Naam'])
        if np.any(mask):
            msg = 'The following names are unknown: {}'.format(names[mask])
            raise(Exception(msg))

    rws.aanvragen_bulk_waarnemingen(locations, email_adress, tmin=tmin,
                                    tmax=tmax)


def waterinfo_to_ghb(fname, riv2stn, gdfv, gwf, model_ds, gdfl=None,
                     intersect_method="vertex"):
    """
    Use waterlevels from waterinfo for a boundary condition, using the ghb
    package. 

    Parameters
    ----------
    fname : string
        The filename of the waterinfo zip-file.
    riv2stn : dict
        A dictionary explaining which stations belong to the indexes in gdfv.
    gdfv : geopandas.GeoDataFrame
        DESCRIPTION.
    gwf : flopy.mf6.ModflowGwf
        The Modflow 6 flopy model to which the ghb-package should be added.
    model_ds : xarray.DataSet
        model dataset, containing information about bathemetry, the layer
        elevation, idomain, and the horizontal conductivity.
    gdfl : geopandas.GeoDataFrame, optional
        Add centerlines of the rivers along which to interpolate the stage.
        Interpolate between the stations when None. The default is None.

    Returns
    -------
    ghb : flopy.mf6.ModflowGwfghb
        The general Head Boundary package that is created for the rivers.

    """

    if isinstance(fname, pd.DataFrame):
        meta = fname
    else:
        if not os.path.isfile(fname):
            msg = 'Connot find file {}. Please assign the right file-name or run {} to request data.'
            raise(Exception(msg.format(fname, 'request_waterinfo_waterlevels')))
        df, meta = rws.read_waterinfo_zip(fname, freq='d', metadata=True)
        # fill nan's
        df = df.interpolate(method='linear', axis=0)
        if np.any(df.fillna('linear').isna()):
            raise(Exception('There are still NaNs in the DataFrame'))
        meta['series'] = None
        for stn in meta.index:
            meta.at[stn, 'series'] = df[stn]

    spd = {}
    if model_ds.steady_start:
        spd[0] = []
        ftp = 1  # first transient period
    else:
        ftp = 0  # first transient period
    spd[ftp] = []

    used_stns = set()
    if 'cellid' not in gdfv.columns:
        gdfv = gdf2grid(gdfv, gwf, method=intersect_method)

    for key in gdfv.index.unique():
        # where along the line are the measuring stations?
        # check if all stations are in the file
        mask = np.array([stn in meta.index for stn in riv2stn[key]])
        if not mask.all():
            msg = 'The following stations are missing: {}'
            # raise(Exception(msg.format(np.array(riv2stn[key])[~mask])))
            print(msg.format(np.array(riv2stn[key])[~mask]))

        if len(riv2stn[key]) == 0:
            raise(Exception(f'No stations supplied for {key}'))
        if len(riv2stn[key]) == 1:
            s = 0.0
        else:
            if (gdfl is not None and key in gdfl.index and
                    isinstance(gdfl.at[key, 'geometry'], LineString)):
                line = gdfl.at[key, 'geometry']
            else:
                # this takes too long and does not always work
                #line = get_centerline_simple(gdfv.at[key,'geometry'], simplify=10)
                coordinates = [(x, y) for x, y in zip(
                    meta.loc[riv2stn[key], 'x'], meta.loc[riv2stn[key], 'y'])]
                line = LineString(coordinates)

            s = [line.project(x)
                 for x in meta.loc[riv2stn[key], 'geometry'].values]
        # make a Series and sort the stations along the line
        ss = pd.Series(s, index=riv2stn[key]).sort_values()

        # where along the line are the river cells?
        for row in gdfv.loc[[key]].itertuples():
            # determine conductance
            cond = row.geometry.area
            if 'bathymetry' in model_ds:
                lays, conds = distribute_cond_over_lays(cond,
                                                        row.cellid,
                                                        model_ds['bathymetry'],
                                                        model_ds['top'],
                                                        model_ds['bot'],
                                                        model_ds['idomain'],
                                                        model_ds['kh'])
            else:
                warnings.warn('No bathymetry found. Adding river to top layer')
                lays = [model_ds.first_active_layer[row.cellid]]
                conds = [cond]

            if np.isnan(conds).any():
                raise ValueError("Conductance is NaN!")

            if len(riv2stn[key]) == 1:
                sm = 0.0
            else:
                sm = line.project(row.geometry.centroid)
            if sm <= ss[0]:
                # the level is equal to that of the first station
                stn = ss.index[0]
                used_stns.add(stn)
                if model_ds.steady_start:
                    for lay, cond in zip(lays, conds):
                        elev = meta.at[stn, 'series'].mean()
                        spd[0].append([(lay,) + row.cellid, elev, cond])
                ts = stn.replace('/', '').replace(' ', '_')
                for lay, cond in zip(lays, conds):
                    spd[ftp].append([(lay,) + row.cellid, ts, cond])
            elif sm >= ss[-1]:
                # the level is equal to that of the last station
                stn = ss.index[-1]
                used_stns.add(stn)
                if model_ds.steady_start:
                    for lay, cond in zip(lays, conds):
                        elev = meta.at[stn, 'series'].mean()
                        spd[0].append([(lay,) + row.cellid, elev, cond])
                ts = stn.replace('/', '').replace(' ', '_')
                for lay, cond in zip(lays, conds):
                    spd[ftp].append([(lay,) + row.cellid, ts, cond])
            else:
                # the level is a combination of two stations
                cross = np.where(ss > sm)[0][0]
                stns = ss.index[cross - 1:cross + 1]
                # determine the multiplication factors of both stations
                fcts = 1 - (np.abs(sm - ss[stns]) /
                            (ss[stns[1]] - ss[stns[0]]))
                for i in range(2):
                    used_stns.add(stns[i])
                    if model_ds.steady_start:
                        for lay, cond in zip(lays, conds):
                            elev = meta.at[stns[i], 'series'].mean()
                            condi = cond * fcts[i]
                            spd[0].append([(lay,) + row.cellid, elev, condi])
                    ts = stns[i].replace('/', '').replace(' ', '_')
                    for lay, cond in zip(lays, conds):
                        condi = cond * fcts[i]
                        spd[ftp].append([(lay,) + row.cellid, ts, condi])

    ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=spd, save_flows=True)
    for i, stn in enumerate(used_stns):
        # make time series
        s = meta.at[stn, 'series']
        index = (s.index - pd.Timestamp(gwf.modeltime.start_datetime)
                 ) / pd.Timedelta(1, 'd')
        timeseries = list(zip(index, s.values))
        stn = stn.replace('/', '').replace(' ', '_')
        ts = {'filename': '{}.ts'.format(stn), 'timeseries': timeseries,
              'time_series_namerecord': stn,
              'interpolation_methodrecord': 'linear'}
        if i == 0:
            ghb.ts.initialize(**ts)
        else:
            ghb.ts.append_package(**ts)
    return ghb


def distribute_cond_over_lays(cond, cellid, rivbot, laytop, laybot,
                              idomain=None, kh=None, stage=None):
    if isinstance(rivbot, np.ndarray) or isinstance(rivbot, xr.DataArray):
        rivbot = float(rivbot[cellid])
    if len(laybot.shape) == 3:
        # the grid is structured grid
        laytop = laytop[cellid[0], cellid[1]]
        laybot = laybot[:, cellid[0], cellid[1]]
        if idomain is not None:
            idomain = idomain[:, cellid[0], cellid[1]]
        if kh is not None:
            kh = kh[:, cellid[0], cellid[1]]
    elif len(laybot.shape) == 2:
        # the grid is a vertex grid
        laytop = laytop[cellid]
        laybot = laybot[:, cellid]
        if idomain is not None:
            idomain = idomain[:, cellid]
        if kh is not None:
            kh = kh[:, cellid]

    if stage is None or isinstance(stage, str):
        lays = np.arange(int(np.sum(rivbot < laybot)) + 1)
    elif np.isfinite(stage):
        lays = np.arange(int(np.sum(stage < laybot)),
                         int(np.sum(rivbot < laybot)) + 1)
    else:
        lays = np.arange(int(np.sum(rivbot < laybot)) + 1)
    if idomain is not None:
        # only distribute conductance over active layers
        lays = lays[idomain.values[lays] > 0]
    topbot = np.hstack((laytop, laybot))
    topbot[topbot < rivbot] = rivbot
    d = -1 * np.diff(topbot)
    if kh is not None:
        kd = kh * d
    else:
        kd = d
    if np.all(kd <= 0):
        # when for some reason the kd is 0 in all layers (for example when the
        # river bottom is above all the layers), add to the first active layer
        if idomain is not None:
            first_active = np.where(idomain == 1)[0][0]
        else:
            first_active = 0
        lays = [first_active]
        kd[first_active] = 1.
    conds = cond * kd[lays] / np.sum(kd[lays])
    return np.array(lays), np.array(conds)


def calc_interpolation_weights(xy, linestring):
    """Calculate interpolation weights for points relative to 
    endpoints linestring.

    Parameters
    ----------
    xy : list
        list of shapely.Point
    linestring : shapely.geometry.LineString
        line to interpolate points

    Returns
    -------
    weights : np.array
        array containing weights per xy point, shape is (Nxy, 2)
    """

    if hasattr(linestring, "geoms"):
        x0, y0 = linestring[0].xy[0][0], linestring[0].xy[1][0]
        xf, yf = linestring[-1].xy[0][-1], linestring[-1].xy[1][-1]
    else:
        x0, y0 = linestring.xy[0][0], linestring.xy[1][0]
        xf, yf = linestring.xy[0][-1], linestring.xy[1][-1]

    zpts = [Point(x0, y0),
            Point(xf, yf)]

    dist_along_ls_z = [linestring.project(z) for z in zpts]

    # projected distance of cell centers along linestring
    dist_along_ls_cells = [linestring.project(p) for p in xy]

    # calculate linear weights between z-pts
    linear_weights = np.interp(
        dist_along_ls_cells, dist_along_ls_z, range(len(dist_along_ls_z)))

    # fill weights matrix
    weights = np.zeros((len(xy), len(zpts)))
    for irow in range(len(xy)):
        w = linear_weights[irow]
        w1 = 1 - (w - int(w)) if w > 0 else 1
        w2 = 1 - w1 if w < len(zpts) else 1
        weights[irow, int(w)] = w1

        if int(w) + 1 < len(zpts):
            weights[irow, int(w) + 1] = w2

    return weights


def aggregate_surface_water(sfw_grid, method, model_ds=None):
    """Aggregate surface water features.

    Parameters
    ----------
    sfw_grid : geopandas.GeoDataFrame
        GeoDataFrame containing surfacewater polygons per grid cell
    method : str, optional
        "de_lange" for De Lange formula for conductance (default),
        "area_weighted" for area-weighted params, 
        "max_area" for max area params
        "individual" for no aggregation
    model_ds : xarray.DataSet, optional
        DataSet containing layer information (only required for 
        method='delange')

    Returns
    -------
    mdata : pd.DataFrame
        DataFrame with aggregated surface water parameters per grid cell
    """
    # params
    c0 = 1.0  # river bottom resistance

    # Post process intersection result
    if method == "individual":
        mdata = pd.DataFrame(index=sfw_grid.cellid)
        stage, cond, rbot = get_surfacewater_params(
            sfw_grid, method, c0=c0)

        mdata.loc[:, "stage"] = stage.values
        mdata.loc[:, "cond"] = cond.values
        mdata.loc[:, "rbot"] = rbot.values

        mdata.loc[:, "area"] = sfw_grid.area.values
        mdata.loc[:, "name_largest"] = sfw_grid.src_id_wla.values

    else:
        gr = sfw_grid.groupby(by="cellid")
        mdata = pd.DataFrame(index=gr.groups.keys())

        for cid, group in tqdm(gr, desc="Aggregate surface water data"):

            stage, cond, rbot = get_surfacewater_params(
                group, method, c0=c0, cid=cid, model_ds=model_ds)

            mdata.loc[cid, "stage"] = stage
            mdata.loc[cid, "cond"] = cond
            mdata.loc[cid, "rbot"] = rbot

            mdata.loc[cid, "area"] = group.area.sum()

            # largest surface water feature
            name_largest = group.loc[group.area.idxmax(), "src_id_wla"]
            mdata.loc[cid, "name_largest"] = name_largest

    return mdata


def build_spd(celldata, pkg, model_ds):
    """Build RIV stress period data.

    Parameters
    ----------
    celldata : geopandas.GeoDataFrame
        GeoDataFrame 
    pkg : str
        Modflow package: RIV or DRN
    model_ds : xarray.DataSet
        DataSet containing model layer information

    Returns
    -------
    spd : list
        list containing stress period data: 
        - RIV: [(cellid), stage, cond, rbot]
        - DRN: [(cellid), stage, cond]
    """

    spd = []

    for cellid, row in tqdm(celldata.iterrows(),
                            total=celldata.index.size,
                            desc=f"Building {pkg} spd"):
        # rbot
        rbot = row["rbot"]

        # stage
        if model_ds.steady_state:
            stage = row["stage"]
            if np.isnan(stage):
                print(f"{cellid}: Cell skipped because stage is NaN")
                continue
            if stage < rbot:
                stage = rbot
        else:
            stage = row["name_largest"]
            if stage is None:
                stage = row["stage"].mean()
            elif isinstance(stage, str):
                stage = (stage.replace(".", "_")
                         .replace(" ", "_")
                         .replace("/", "_"))
            elif np.isnan(stage):
                continue

        # conductance
        cond = row["cond"]

        # check value
        if np.isnan(cond):
            print(f"{cellid}: Conductance is NaN! Info: area={row.area:.2f} "
                  f"len={row.len_estimate:.2f}, BL={row['rbot']}")
            continue
        if cond < 0:
            print(f"{cellid}, Conductance is < 0!, area={row.area:.2f}, "
                  f"len={row.len_estimate:.2f}, BL={row['rbot']}")
            continue

        lays, conds = distribute_cond_over_lays(cond,
                                                cellid,
                                                rbot,
                                                model_ds.top,
                                                model_ds.bot,
                                                model_ds.idomain,
                                                model_ds.kh,
                                                stage)
        # write SPD
        for lay, cond in zip(lays, conds):
            cid = (lay,) + cellid
            if pkg == "RIV":
                spd.append([cid, stage, cond, rbot])
            elif pkg == "DRN":
                spd.append([cid, stage, cond])

    return spd


def create_timeseries(celldata, t_start, t_end):
    """Create timeseries based on summer/winter levels between t_start and 
    t_end.

    Parameters
    ----------
    celldata : geopandas.GeoDataFrame
        DataFrame containing summer/winter levels
    t_start : str or pd.Timestamp
        start time
    t_end : str or pd.Timestamp]
        end time

    Returns
    -------
    tseries_list : list
        list containing timeseries
    """
    tseries_list = []

    # for peilvak in sfw.src_id_wla.unique():
    for peilvak in tqdm(celldata.name_largest.unique(),
                        desc="Generating timeseries"):
        if peilvak is None:
            continue
        if isinstance(peilvak, float):
            if np.isnan(peilvak):
                continue
        peilen = celldata.loc[celldata["name_largest"] == peilvak,
                              ["ZP", "WP"]].iloc[0]

        if peilen.isna().any():
            continue

        dt = pd.date_range(t_start, t_end, freq="MS")
        dt_apr_oct = [i for i in dt if i.month in [4, 10]]
        doffset = pd.tseries.offsets.DateOffset(months=6)
        dt_apr_oct.insert(0, dt_apr_oct[0] - doffset)
        dt_apr_oct.append(dt_apr_oct[-1] + doffset)
        dt = pd.DatetimeIndex(dt_apr_oct)
        dtnum = ((dt - t_start).days).to_numpy()
        dtnum[dtnum < 0] = 0.0
        ts = pd.Series(index=dtnum, dtype=float)
        ts.where(dt.month == 4, peilen.ZP, inplace=True)
        ts.where(dt.month == 10, peilen.WP, inplace=True)
        ts.name = (peilvak.replace(".", "_")
                   .replace(" ", "_")
                   .replace("/", "_"))
        tseries_list.append(ts)
    return tseries_list


def get_surfacewater_params(group, method, c0=1.0, cid=None, model_ds=None):

    if method == "area_weighted":
        zp = agg_area_weighted(group, "ZP")
        wp = agg_area_weighted(group, "WP")
        stage = np.mean([zp, wp])
        # cond
        cond = group.area.sum() / c0
        # rbot
        rbot = group["BL"].min()

    elif method == "max_area":
        # stage
        zp = agg_max_area(group, "ZP")
        wp = agg_max_area(group, "WP")
        stage = np.mean([zp, wp])
        # cond
        cond = group.area.sum() / c0
        # rbot
        rbot = group["BL"].min()

    elif method == "de_lange":
        c1 = 0.0
        N = 1e-3
        # stage
        zp = agg_area_weighted(group, "ZP")
        wp = agg_area_weighted(group, "WP")
        stage = np.mean([wp, zp])
        # cond
        _, _, cond = agg_de_lange(group, cid, model_ds, c1=c1, c0=c0, N=N)
        # rbot
        rbot = group["BL"].min()

    elif method == "individual":
        # stage = group["src_id_wla"]
        stage = group.loc[:, ["ZP", "WP"]].mean(axis=1)
        cond = group.area / c0
        rbot = group["BL"]
    else:
        raise ValueError(f"Method '{method}' not recognized!")

    return stage, cond, rbot


def get_subsurface_params(model_ds, cid):
    r, c = cid

    A = model_ds.delr * model_ds.delc  # cell area
    laytop = model_ds['top'].isel(x=c, y=r).data
    laybot = model_ds['bot'].isel(x=c, y=r).data
    kv = model_ds['kv'].isel(x=c, y=r).data
    kh = model_ds['kh'].isel(x=c, y=r).data
    # first_active_lay = model_ds["first_active_layer"].isel(x=c, y=r).data
    thickness = model_ds["thickness"].isel(x=c, y=r).data
    # idomain = model_ds["idomain"].isel(x=c, y=r).data
    return A, laytop, laybot, kh, kv, thickness


def agg_max_area(gdf, col):
    return gdf.loc[gdf.area.idxmax(), col]


def agg_area_weighted(gdf, col):
    nanmask = gdf[col].isna()
    aw = ((gdf.area * gdf[col]).sum(skipna=True) /
          gdf.loc[~nanmask].area.sum())
    return aw


def agg_de_lange(group, cid, model_ds, c1=0.0, c0=1.0, N=1e-3,
                 crad_positive=True):

    (A, laytop, laybot, kh, kv, thickness) = get_subsurface_params(model_ds,
                                                                   cid)

    rbot = group["BL"].min()

    # select active layers
    active = thickness > 0
    laybot = laybot[active]
    kh = kh[active]
    kv = kv[active]
    thickness = thickness[active]

    # layer thickn.
    H0 = laytop - laybot[laybot < rbot][0]
    ilay = 0
    rlay = np.where(laybot < rbot)[0][0]

    # equivalent hydraulic conductivities
    H = thickness[ilay:rlay + 1]
    kv = kv[ilay:rlay + 1]
    kh = kh[ilay:rlay + 1]
    kveq = np.sum(H) / np.sum(H / kv)
    kheq = np.sum(H * kh) / np.sum(H)

    # length
    len_est = estimate_polygon_length(group)
    li = len_est.sum()
    # correction if group contains multiple shapes
    # but covers whole cell
    if group.area.sum() == A:
        li = A / np.max([model_ds.delr, model_ds.delc])

    # width
    B = group.area.sum(skipna=True) / li

    # mean water level
    p = group.loc[group.area.idxmax(), ["ZP", "WP"]].mean()  # waterlevel

    # calculate params
    pstar, cstar, cond = de_lange_eqns(
        A, H0, kveq, kheq, c1, li, B, c0, p, N, crad_positive=crad_positive)

    return pstar, cstar, cond


def de_lange_eqns(A, H0, kv, kh, c1, li, Bin, c0, p, N, crad_positive=True):
    """Calculates the conductance according to De Lange

    Parameters
    ----------
    A : float
        celoppervlak (m2)
    H0 : float
        doorstroomde dikte (m)
    kv : float
        verticale doorlotendheid (m/d)
    kh : float
        horizontale doorlatendheid (m/d)
    c1 : float
        deklaagweerstand (d)
    li : float
        lengte van de waterlopen (m)
    Bin : float
        bodembreedte (m)
    c0 : float
        slootbodemweerstand (d)
    p : float
        water peil
    N : float
        grondwateraanvulling
    crad_positive: bool, optional
        whether to allow negative crad values. If True, crad will be set to 0
        if it is negative.

    Returns
    -------
    float
        Conductance (m2/d)

    """
    if li > 1e-3 and Bin > 1e-3 and A > 1e-3:
        Bcor = max(Bin, 1e-3)  # has no effect
        L = A / li - Bcor
        y = c1 + H0 / kv

        labdaL = np.sqrt(y * kh * H0)
        if L > 1e-3:
            xL = L / (2 * labdaL)
            FL = xL * coth(xL)
        else:
            FL = 0.0

        labdaB = np.sqrt(y * kh * H0 * c0 / (y + c0))
        xB = Bcor / (2 * labdaB)
        FB = xB * coth(xB)

        CL = (c0 + y) * FL + (c0 * L / Bcor) * FB
        if CL == 0.0:
            CB = 1.0
        else:
            CB = (c1 + c0 + H0 / kv) / (CL - c0 * L / Bcor) * CL

        # volgens Kees Maas mag deze ook < 0 zijn...
        # er miste ook een correctie in de log voor anisotropie
        # Crad = max(0., L / (np.pi * np.sqrt(kv * kh))
        #            * np.log(4 * H0 / (np.pi * Bcor)))
        crad = radial_resistance(L, Bcor, H0, kh, kv)
        if crad_positive:
            crad = max([0.0, crad])

        # Conductance
        pSl = Bcor * li / A
        if pSl >= 1.0 - 1e-10:
            Wp = 1 / (pSl / CB) + crad - c1
        else:
            Wp = 1 / ((1. - pSl) / CL + pSl / CB) + crad - c1
        cond = A / Wp

        # cstar, pstar
        cLstar = CL + crad

        pstar = p + N * (cLstar - y) * (y + c0) * L / (Bcor * cLstar + L * y)
        cstar = cLstar * (c0 + y) * (Bcor + L) / (Bcor * cLstar + L * y)

        return pstar, cstar, cond
    else:
        return 0., 0., 0.


def radial_resistance(L, B, H, kh, kv):
    return (L / (np.pi * np.sqrt(kh * kv)) *
            np.log(4 * H * np.sqrt(kh) / (np.pi * B * np.sqrt(kv))))


def coth(x):
    return 1.0 / np.tanh(x)


def estimate_polygon_length(gdf):
    # estimate length from polygon (for shapefactor > 4)
    shape_factor = gdf.length / np.sqrt(gdf.area)

    len_est1 = (gdf.length - np.sqrt(gdf.length**2 - 16 * gdf.area)) / 4
    len_est2 = (gdf.length + np.sqrt(gdf.length**2 - 16 * gdf.area)) / 4
    len_est = pd.concat([len_est1, len_est2], axis=1).max(axis=1)

    # estimate length from minimum rotated rectangle (for shapefactor < 4)
    min_rect = gdf.geometry.apply(lambda g: g.minimum_rotated_rectangle)
    xy = min_rect.apply(lambda g: np.sqrt(
        (np.array(g.exterior.xy[0]) - np.array(g.exterior.xy[0][0]))**2 +
        (np.array(g.exterior.xy[1]) - np.array(g.exterior.xy[1][0]))**2))
    len_est3 = xy.apply(lambda a: np.partition(a.flatten(), -2)[-2])

    # update length estimate where shape factor is lower than 4
    len_est.loc[shape_factor < 4] = len_est3.loc[shape_factor < 4]

    return len_est


def get_bc_type(df):

    maskriv_inf = ((df.SUB == 1) &
                   (df.has_slope == 0) &
                   (df.src_wla != "rws_krw") &
                   (df.CAT != 3))
    #    (~df.src_id_wla.isna()) &
    #    (~np.isnan(df.BL)))

    maskriv_drn = ((df.SUB == 0) &
                   (df.has_slope == 0) &
                   (df.src_wla != "rws_krw") &
                   (df.CAT != 3))
    #    (~df.src_id_wla.isna()) &
    #    (~np.isnan(df.BL)))

    maskdry = ((df.has_slope == 0) &
               (df.CAT == 3) &
               (df.src_wla != "rws_krw"))
    #    (~df.src_id_wla.isna()) &
    #    (~np.isnan(df.BL)))

    maskslope = (df.has_slope == 1)
    #  (~df.loc[:, ["ZP", "WP"]].isna().any(axis=1)))

    maskrws = (df.src_wla == "rws_krw")

    maskall = maskriv_inf | maskriv_drn | maskdry | maskslope | maskrws

    bc = pd.Series(index=df.index, dtype=object)
    bc.loc[maskriv_inf] = "riv_inf"
    bc.loc[maskriv_drn] = "riv_drn"
    bc.loc[maskdry] = "riv_drn"
    bc.loc[maskslope] = "riv_slp"
    bc.loc[maskrws] = "ghb"
    bc.loc[~maskall] = "none"

    return bc

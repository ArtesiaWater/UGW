# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 10:59:35 2020

NOTES
-----
for now the time factor of the bathymetry is ignored. For each tile the 
latest available bathymetry is used. 


"""
import os

import numpy as np
import requests
import xarray as xr

from . import mgrid, util


def get_modelgrid_bathymetry(model_ds, extent, datadir,
                             gridtype='structured',
                             gridprops=None,
                             cachedir=None, use_cache=False,
                             verbose=False):
    """ get bathymetry of the Northsea from the jarkus dataset.

    Parameters
    ----------
    model_ds : xarray.Dataset
        dataset with model data where bathymetry is added to
    extent : list, tuple or np.array
        model_extent (xmin, xmax, ymin, ymax)
    cachedir : str, optional
        directory to store cached values, if None a temporary directory is
        used. default is None
    use_cache : bool, optional
        if True the cached jarkus data is used. The default is False.
    verbose : bool, optional
        print additional information to the screen. The default is False.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with bathymetry
    """

    model_ds = util.get_cache_netcdf(use_cache, cachedir, 'bathymetry_model_ds.nc',
                                     bathymetry_to_model_dataset, model_ds,
                                     verbose=verbose,
                                     extent=extent, datadir=datadir,
                                     gridprops=gridprops,
                                     check_time=False,
                                     )
    return model_ds


def bathymetry_to_model_dataset(model_ds, extent, datadir,
                                gridprops=None):
    """ get bathymetry of the Northsea from the jarkus dataset.

    Parameters
    ----------
    model_ds : xarray.Dataset
        dataset with model data where bathymetry is added to
    extent : list, tuple or np.array
        model_extent (xmin, xmax, ymin, ymax)
    datadir : str, optional
        directory where jarkus data is stored, if None a temporary directory is
        used. default is None
    use_cache : bool, optional
        if True the cached jarkus data is used. The default is False.
    verbose : bool, optional
        print additional information to the screen. The default is False.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with bathymetry

    Notes
    -----
    Now the nan values in the original bathymetry are filled and than the 
    data is resampled to the modelgrid. Maybe we can speed up things by
    changing the order in which operations are executed.
    """
    pre_processed_jarkus = os.path.join(datadir, 'bathymetry',
                                        'jarkus_nhflopy.nc')

    if os.path.isfile(pre_processed_jarkus):
        ds_bathymetry_raw = xr.open_dataset(pre_processed_jarkus)
    else:
        ds_bathymetry_raw = get_dataset_jarkus(extent, datadir, use_cache=True)

    da_bathymetry_raw = ds_bathymetry_raw['z']

    # fill nan values in bathymetry
    da_bathymetry_filled = mgrid.fillnan_dataarray_structured_grid(
        da_bathymetry_raw)

    # bathymetrie mag nooit groter zijn dan NAP 0.0
    da_bathymetry_filled = xr.where(
        da_bathymetry_filled > 0, 0, da_bathymetry_filled)

    # bathymetry projected on model grid

    if model_ds.gridtype == 'structured':
        da_bathymetry = mgrid.resample_dataarray_to_structured_grid(da_bathymetry_filled,
                                                                    extent=extent,
                                                                    delr=model_ds.delr,
                                                                    delc=model_ds.delc,
                                                                    xmid=model_ds.x.data,
                                                                    ymid=model_ds.y.data[::-1])[0]
    elif model_ds.gridtype == 'unstructured':
        da_bathymetry = mgrid.resample_dataarray_to_unstructured_grid(da_bathymetry_filled,
                                                                      gridprops=gridprops)[0]

    model_ds_out = util.get_model_ds_empty(model_ds)

    model_ds_out['bathymetry'] = xr.where(
        model_ds['sea'], da_bathymetry, np.nan)

    return model_ds_out


def get_jarkus_tilenames(extent):
    """ Find all Jarkus tilenames within a certain extent

    Parameters
    ----------
    extent : list, tuple or np.array
        extent (xmin, xmax, ymin, ymax) of the desired grid. Should be RD-new
        coördinates (EPSG:28992)

    Returns
    -------
    netcdf_urls : list of str
        list of the urls of all netcdf files of the tiles with Jarkus data.

    """
    ds_jarkus_catalog = xr.open_dataset(
        'http://opendap.deltares.nl/thredds/dodsC/opendap/rijkswaterstaat/jarkus/grids/catalog.nc')
    ew_x = ds_jarkus_catalog['projectionCoverage_x']
    sn_y = ds_jarkus_catalog['projectionCoverage_y']

    mask_ew = (ew_x[:, 1] > extent[0]) & (ew_x[:, 0] < extent[1])
    mask_sn = (sn_y[:, 1] > extent[2]) & (sn_y[:, 0] < extent[3])

    indices_tiles = np.where(mask_ew & mask_sn)[0]
    all_netcdf_tilenames = get_netcdf_tiles()

    netcdf_tile_names = [all_netcdf_tilenames[i] for i in indices_tiles]

    return netcdf_tile_names


def get_netcdf_tiles():
    """ Find all Jarkus netcdf tile names. 


    Returns
    -------
    netcdf_urls : list of str
        list of the urls of all netcdf files of the tiles with Jarkus data.

    Notes
    -----
    This function would be redundant if the jarkus catalog 
    (http://opendap.deltares.nl/thredds/dodsC/opendap/rijkswaterstaat/jarkus/grids/catalog.nc)
    had a proper way of displaying the url's of each tile. It seems like an
    attempt was made to do this because here is a data variable
    named 'urlPath' in the catalog. However the dataarray of 'urlPath' has the
    same string for each tile.
    """
    url = 'http://opendap.deltares.nl/thredds/dodsC/opendap/rijkswaterstaat/jarkus/grids/catalog.nc.ascii'
    req = requests.get(url)
    s = req.content.decode('ascii')
    start = s.find('urlPath', s.find('urlPath') + 1)
    end = s.find('projectionCoverage_x', s.find('projectionCoverage_x') + 1)
    netcdf_urls = list(eval(s[start + 12:end - 2]))
    return netcdf_urls


def get_dataset_jarkus(extent, datadir, use_cache=True, verbose=True):
    """ Get bathymetry from Jarkus within a certain extent. If all_tiles is 
    False, the following steps are used:
    1. find Jarkus tiles within the extent
    2. download netcdf files of Jarkus tiles
    3. read Jarkus tiles and combine the 'z' parameter of the last time step
    of each tile, to a dataarray.

    Parameters
    ----------
    extent : list, tuple or np.array
        extent (xmin, xmax, ymin, ymax) of the desired grid. Should be RD-new
        coördinates (EPSG:28992)
    datadir : str
        directory where the Jarkus nc files are stored.
    use_cache : bool, optional
        if True the cached resampled dataset is used. 
        The default is True.
    verbose : bool, optional
        print additional information to the screen. The default is False.

    Returns
    -------
    z : xr.DataSet
        dataset containing bathymetry data

    """

    extent = [int(x) for x in extent]
    if use_cache:
        fname_bathymetry = 'bathymetry_{}_{}_{}_{}.nc'.format(*extent)
        path_bathymetry = os.path.join(datadir, fname_bathymetry)
        if os.path.exists(path_bathymetry):
            z_dataset = xr.open_dataset(path_bathymetry)
            return z_dataset

    netcdf_tile_names = get_jarkus_tilenames(extent)
    tiles = [xr.open_dataset(name) for name in netcdf_tile_names]
    # only use the last timestep
    tiles = [tile.isel(time=-1) for tile in tiles]
    z_dataset = xr.combine_by_coords(tiles, combine_attrs='drop')
    if use_cache:
        z_dataset.to_netcdf(path_bathymetry)

    return z_dataset

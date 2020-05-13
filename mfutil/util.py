# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:11:03 2020

@author: oebbe
"""

import os
import tempfile

import geopandas as gpd
import numpy as np
import xarray as xr
from shapely.geometry import box

import art_tools
import mgrid


def get_model_ds_empty(model_ds):
    """ get a copy of a model dataset with only grid and time information.

    Parameters
    ----------
    model_ds : xr.Dataset
        dataset with at least the variables layer, x, y and time

    Returns
    -------
    model_ds_out : xr.Dataset
        dataset with only model grid and time information

    """
    if model_ds.gridtype == 'structured':
        model_ds_out = model_ds[['layer', 'x', 'y', 'time']].copy()
        return model_ds_out
    elif model_ds.gridtype == 'unstructured':
        model_ds_out = model_ds[['cid', 'layer', 'x', 'y', 'time']].copy()
        return model_ds_out
    else:
        raise ValueError('no gridtype defined cannot compare model datasets')


def check_model_ds(model_ds, model_ds2):
    """ check if two model datasets have the same grid and time discretization.
    e.g. the same dimensions and coordinates.


    Parameters
    ----------
    model_ds : xr.Dataset
        dataset with model grid and time discretisation
    model_ds2 : xr.Dataset
        dataset with model grid and time discretisation. This is typically
        a cached dataset.

    Raises
    ------
    ValueError
        if the gridtype of model_ds is not structured or unstructured.

    Returns
    -------
    bool
        True if the two datasets have the same grid and time discretization.

    """

    if model_ds.gridtype == 'structured':
        try:
            check_x = (model_ds['x'] == model_ds2['x']).all()
            check_y = (model_ds['y'] == model_ds2['y']).all()
            check_layer = (model_ds['layer'] == model_ds2['layer']).all()
            check_time = (model_ds['time'] == model_ds2['time']).all()
            if check_x and check_y and check_layer and check_time:
                return True
        except KeyError:
            return False
    elif model_ds.gridtype == 'unstructured':
        try:
            check_cid = (model_ds['cid'] == model_ds2['cid']).all()
            check_x = (model_ds['x'] == model_ds2['x']).all()
            check_y = (model_ds['y'] == model_ds2['y']).all()
            check_layer = (model_ds['layer'] == model_ds2['layer']).all()
            check_time = (model_ds['time'] == model_ds2['time']).all()
            if check_cid and check_x and check_y and check_layer and check_time:
                return True
            else:
                return False
        except KeyError:
            return False
    else:
        raise ValueError('not gridtype defined cannot compare model datasets')


def get_cache_netcdf(use_cache, cachedir, cache_name, get_dataset_func,
                     model_ds=None, verbose=False, **get_kwargs):
    """ function is used to create, read or modify cached variables of a 
    model dataset.

    following steps are done:
        1. Read cached dataset and merge this with the current model_ds if all 
        of the following conditions are met:
            a. use_cache = True
            b. dataset exists in cachedir
            c. the grid and time discretisation of the cached dataset equals
            the grid and time discretisation of the model dataset
        2. if the conditions in step 1 are not met the get_dataset_func is
        called (with the **get_kwargs arguments).
        3. the dataset from step 2 is written to the cachedir.
        4. the dataset from step 2 is merged with the current model_ds.


    Parameters
    ----------
    use_cache : bool
        if True an attempt is made to use the cached dataset.
    cachedir : str
        directory to store cached values, if None a temporary directory is
        used. default is None
    cache_name : str
        named of the cached netcdf file with the dataset.
    get_dataset_func : function
        this function is called to obtain a new dataset (and not use the 
        cached dataset).
    model_ds : xr.Dataset
        dataset where the cached or new dataset is added to.
    verbose : bool, optional
        print additional information. default is False
    **get_kwargs : 
        keyword arguments are used when calling the get_dataset_func.

    Returns
    -------
    model_ds
        dataset with the cached or new dataset.

    """

    if cachedir is None:
        cachedir = tempfile.gettempdir()

    fname_model_ds = os.path.join(cachedir, cache_name)

    if use_cache:
        if os.path.exists(fname_model_ds):
            if verbose:
                print(f'found cached {cache_name}, loading cached dataset')

            cache_model_ds = xr.open_dataset(fname_model_ds)
            if model_ds is None:
                return cache_model_ds

            elif check_model_ds(model_ds, cache_model_ds):
                if verbose:
                    print('cached data has same grid as current model\n')
                model_ds.update(cache_model_ds)
                cache_model_ds.close()
                return model_ds
            else:
                cache_model_ds.close()
                if verbose:
                    print('cached data grid differs from current model')
    if verbose:
        print(f'creating and caching dataset {cache_name}\n')
    if model_ds is None:
        ds = get_dataset_func(**get_kwargs)
        ds.to_netcdf(fname_model_ds)
        return ds
    else:
        ds = get_dataset_func(model_ds, **get_kwargs)
        ds.to_netcdf(fname_model_ds)
        model_ds.update(ds)
        return model_ds


def get_ahn_dataset(model_ds, extent, resolution, use_cache=True, cachedir=None,
                    fname_netcdf='ahn_model_ds.nc', verbose=False,
                    gridprops=None):
    """ get an xarray dataset from the ahn values within an extent.

    Parameters
    ----------
    extent : list, tuple or np.array
        extent (xmin, xmax, ymin, ymax) of the desired grid.
    resolution : int or float
        resolution of the ahn grid.
    use_cache : bool, optional
        if True the cached resampled regis dataset is used. 
        The default is False.

    Returns
    -------
    ahn_ds : xarray.Dataset
        dataset with ahn data (not yet corresponding to the modelgrid)

    Notes
    -----
    1. The resolution is used to obtain the ahn from the wcs server. Not sure
    what kind of interpolation is used to resample the original grid.

    2. Still a dependency on art_tools to get the ahn.

    3. The ahn raster is now cached in a tempdir. Should be changed to the
    cachedir of the model I think.

    """
    ahn_ds = get_cache_netcdf(use_cache, cachedir, fname_netcdf,
                              get_ahn_at_grid, model_ds,
                              verbose=verbose, extent=extent,
                              resolution=resolution, gridprops=gridprops)

    return ahn_ds


def get_ahn_at_grid(model_ds, extent, resolution, gridprops=None):
    """ Get a model dataset with ahn variable.


    Parameters
    ----------
    model_ds : xr.Dataset
        dataset with the grid information.
    extent : list, tuple or np.array
        extent (xmin, xmax, ymin, ymax) of the desired grid.
    resolution : int or float
        resolution of the ahn grid.
    gridprops : dictionary, optional
        dictionary with grid properties output from gridgen. Only used if
        gridtype = 'unstructured'

    Returns
    -------
    model_ds_out : xr.Dataset
        dataset with the ahn variable.

    """

    fname_ahn = art_tools.get_ahn_within_extent(extent=extent, res=resolution,
                                                return_fname=True,
                                                cache=True)

    ahn_ds_raw = xr.open_rasterio(fname_ahn)
    ahn_ds_raw = ahn_ds_raw.rename({'band': 'layer'})
    ahn_ds_raw = xr.where(ahn_ds_raw > 1e35, np.nan, ahn_ds_raw)

    if model_ds.gridtype == 'structured':
        ahn_ds = mgrid.resample_dataarray_to_structured_grid(ahn_ds_raw,
                                                             extent=extent,
                                                             delr=model_ds.delr,
                                                             delc=model_ds.delc,
                                                             xmid=model_ds.x.data,
                                                             ymid=model_ds.y.data[::-1])
    elif model_ds.gridtype == 'unstructured':
        xyi, cid = mgrid.get_xyi_cid(gridprops)
        ahn_ds = mgrid.resample_dataarray_to_unstructured_grid(ahn_ds_raw,
                                                               gridprops,
                                                               xyi, cid)

    model_ds_out = get_model_ds_empty(model_ds)
    model_ds_out['ahn'] = ahn_ds[0]

    return model_ds_out


def get_regis_dataset(regis_ds_raw=None, extent=None,
                      delr=None, delc=None,
                      gridtype='structured',
                      gridprops=None,
                      interp_method="linear",
                      cachedir=None,
                      fname_netcdf='regis.nc',
                      use_cache=False,
                      verbose=False):
    """Get a regis dataset that is resampled to the modelgrid

    Parameters
    ----------
    regis_ds_raw : xarray.Dataset, optional
        dataset with raw regis data (netcdf from url). The default is None.
    extent : list, tuple or np.array
        extent (xmin, xmax, ymin, ymax) of the desired grid.
    delr : int or float
        cell size along rows of the desired grid (dx).
    delc : int or float
        cell size along columns of the desired grid (dy).
    gridtype : str, optional
        type of grid, options are 'structured' and 'unstructured'. 
        The default is 'structured'.
    gridprops : dictionary, optional
        dictionary with grid properties output from gridgen. Only used if
        gridtype = 'unstructured'
    interp_method : str, optional
        interpolation method, default is 'linear'
    cachedir : str, optional
        directory to store cached values, if None a temporary directory is
        used. default is None
    fname_netcdf : str
        name of the netcdf file that is stored in the cachedir.
    use_cache : bool, optional
        if True the cached resampled regis dataset is used. 
        The default is False.

    Returns
    -------
    regis_ds : xarray.Dataset
        dataset with regis data corresponding to the modelgrid

    """

    regis_ds = get_cache_netcdf(use_cache, cachedir, fname_netcdf,
                                get_regis_ds_from_raw_ds,
                                verbose=verbose,
                                regis_ds_raw=regis_ds_raw, extent=extent,
                                delr=delr, delc=delc, gridtype=gridtype, 
                                kind=interp_method)

    return regis_ds


def get_regis_ds_from_raw_ds(regis_ds_raw=None,
                             extent=None, delr=None, delc=None,
                             gridtype='structured',
                             gridprops=None, verbose=False,
                             kind='linear',
                             ):
    """ project regis data on to the modelgrid.


    Parameters
    ----------
    regis_ds_raw : xr.Dataset, optional
        regis dataset. If the gridtype is structured this should be the
        original regis netcdf dataset. If gridtype is unstructured this should
        be a regis dataset that is already projected on a structured modelgrid. 
        The default is None.
    extent : list, tuple or np.array
        extent (xmin, xmax, ymin, ymax) of the desired grid.
    delr : int or float
        cell size along rows of the desired grid (dx).
    delc : int or float
        cell size along columns of the desired grid (dy).
    gridtype : str, optional
        type of grid, options are 'structured' and 'unstructured'. 
        The default is 'structured'.
    gridprops : dictionary, optional
        dictionary with grid properties output from gridgen. Only used if
        gridtype = 'unstructured'
    verbose : bool, optional
        print additional information. default is False
    kind : str, optional
        kind of interpolation to use. default is 'linear'

    Returns
    -------
    regis_ds : xr.dataset
        regis dataset projected onto the modelgrid.

    """

    if gridtype == 'structured':
        if verbose:
            print(f'resample regis data to structured modelgrid')
        regis_ds = mgrid.resample_dataset_to_structured_grid(regis_ds_raw, extent,
                                                             delr, delc, kind=kind)
        regis_ds.attrs['extent'] = extent
        regis_ds.attrs['delr'] = delr
        regis_ds.attrs['delc'] = delc
        regis_ds.attrs['gridtype'] = gridtype
    elif gridtype == 'unstructured':
        if verbose:
            print(f'resample regis data to unstructured modelgrid')
        regis_ds = mgrid.resample_dataset_to_unstructured_grid(
            regis_ds_raw, gridprops)
        regis_ds['x'] = xr.DataArray([r[1] for r in gridprops['cell2d']],
                                     dims=('cid'),
                                     coords={'cid': regis_ds.cid.data})

        regis_ds['y'] = xr.DataArray([r[2] for r in gridprops['cell2d']],
                                     dims=('cid'),
                                     coords={'cid': regis_ds.cid.data})
        regis_ds.attrs['gridtype'] = gridtype
        regis_ds.attrs['extent'] = regis_ds_raw.extent

    return regis_ds


def find_most_recent_file(folder, name, extension='.pklz'):
    """ find the most recent file in a folder. File must startwith name and
    end width extension. If you want to look for the most recent folder use
    extension = ''.

    Parameters
    ----------
    folder : str
        path of folder to look for files
    name : str
        find only files that start with this name
    extension : str
        find only files with this extension

    Returns
    -------
    newest_file : str
        name of the most recent file
    """

    i = 0
    for file in os.listdir(folder):
        if file.startswith(name) and file.endswith(extension):
            if i == 0:
                newest_file = os.path.join(folder, file)
                time_prev_file = os.stat(newest_file).st_mtime
            else:
                check_file = os.path.join(folder, file)
                if os.stat(check_file).st_mtime > time_prev_file:
                    newest_file = check_file
                    time_prev_file = os.stat(check_file).st_mtime
            i += 1

    if i == 0:
        return None

    return newest_file


def gdf_within_extent(gdf, extent):
    """ select only parts of the geodataframe within the extent.
    Only works for polygon features.

    Parameters
    ----------
    gdf : geopandas GeoDataFrame
        dataframe with polygon features.
    extent : list or tuple
        extent to slice gdf, (xmin, xmax, ymin, ymax).

    Returns
    -------
    gdf : geopandas GeoDataFrame
        dataframe with only polygon features within the extent.

    """

    bbox = (extent[0], extent[2], extent[1], extent[3])
    geom_extent = box(*tuple(bbox))
    gdf_extent = gpd.GeoDataFrame(['extent'], geometry=[geom_extent],
                                  crs=gdf.crs)
    gdf = gpd.overlay(gdf, gdf_extent)

    return gdf

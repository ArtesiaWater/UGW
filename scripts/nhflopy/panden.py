# -*- coding: utf-8 -*-
"""
Created on Wed May 27 15:47:18 2020

@author: oebbe
"""
import os
import flopy
import pandas as pd
import geopandas as gpd
import numpy as np
import datetime as dt
import xarray as xr
import warnings

from . import mgrid
from . import util

def get_panden_shp(datadir, fname='Panden_ICAS_IKIEF.shp'):
    """ read shapefiles with locations of the infiltration ponds
    

    Parameters
    ----------
    datadir : str
        directory where the shapefile is located in the panden directory.
    fname : str, optional
        Name of the shapefile with the polygons of the infiltration ponds.
    
    Returns
    -------
    panden_shp : GeoDataFrame
       locations of the infiltration ponds.

    """
    panden_shp = gpd.read_file(os.path.join(datadir, 'panden', fname))
    panden_shp = panden_shp.explode().reset_index(drop=True)
    assert np.all(panden_shp.geom_type=='Polygon'),'For now only polygons are used'

    panden_shp.loc[panden_shp.ICAS==0, 'peil'] = 'pp_pand9'
    panden_shp.loc[panden_shp.ICAS!=0, 'peil'] = 'pp_mkom'
    panden_shp.loc[:, 'bweerstand'] = 1
    
    return panden_shp

def get_peilen(datadir, fname='PREP2XL_v122_Pandpeilen_kopie.xlsx'):
    """ read excel file with time series of the stage in the infiltration 
    ponds
    

    Parameters
    ----------
    datadir : str
        directory where the excel file is located in the panden directory.
    fname : str, optional
        Name of the excel file with the stages
        
    Returns
    -------
    peilen : DataFrame
        time series of infiltration ponds stage.

    """
    
    peilen = pd.read_excel(os.path.join(datadir, 'panden', fname),
                           skiprows=list(range(7))+[8,9], engine='openpyxl',
                           usecols=[0,8,10], index_col=0)
    peilen.rename(columns={'Pandniveau pand 9': 'pp_pand9',
                           'Pandniveau middenkom':'pp_mkom'}, inplace=True)
    
    return peilen

def get_panden(model_ds, panden, peilen, modelgrid, cachedir=None, 
               use_cache=False,
               verbose=False):
    """ Get datasets from dataframes with panden en peilen.

    Parameters
    ----------
    model_ds : xr.DataSet
        dataset containing relevant model information
    panden : GeoDataFrame
        GeoDataframe with locations of the infiltration ponds.
    peilen :  pandas DataFrame
        Dataframe with the ponds levels over time.
    modelgrid : flopy grid
        model grid.
    cachedir : str, optional
        directory to store cached values, if None a temporary directory is
        used. default is None
    use_cache : bool, optional
        if True the cached riv_panden data is used. The default is False.
    verbose : bool, optional
        print additional information to the screen. The default is False.

    Returns
    -------
    model_ds : xr.DataSet
        dataset with spatial model data including the riv rasters and 
        timeseries

    """
    model_ds = util.get_cache_netcdf(use_cache, cachedir, 'panden_model_ds.nc',
                                     add_panden_to_model_ds, 
                                     model_ds,
                                     panden=panden, 
                                     peilen=peilen, 
                                     modelgrid=modelgrid, 
                                     verbose=verbose)

    return model_ds
    
def add_panden_to_model_ds(model_ds, panden, peilen, modelgrid, 
                           name='pand_riv',
                           diepte_pand_gem = 2.0, verbose=False):
    """
    

    Parameters
    ----------
    model_ds : xr.DataSet
        dataset containing relevant model information
    panden : GeoDataFrame
        GeoDataframe with locations of the infiltration ponds.
    peilen :  pandas DataFrame
        Dataframe with the ponds levels over time.
    modelgrid : flopy grid
        model grid.
    name : str, optional
        name, used as variable name in model_ds. The default is 'panden'.
    diepte_pand_gem : float, optional
        depth of the ponds. Used to add as river bottom to the model. 
        The default is 2.0.
    verbose : bool, optional
        print additional information to the screen. The default is False.


    Returns
    -------
    model_ds_out : xr.DataSet
        dataset with the riv rasters and timeseries


    """
    
    if model_ds.gridtype != 'unstructured':
        warnings.warn(
            'this function is not yet tested for structured grids')

    area = xr.zeros_like(model_ds['top'])
    cond = xr.zeros_like(model_ds['top'])
    bot = xr.zeros_like(model_ds['top'])
    peil = xr.full_like(model_ds['top'], '' , dtype="S13")
    peil = peil.astype(str)
    
    for i, row in panden.iterrows():
        area_pol = mgrid.polygon_to_area(modelgrid, row['geometry'],
                                         xr.ones_like(model_ds['top']),
                                         model_ds.gridtype)
        cond = xr.where(area_pol != 0, cond+area_pol / row['bweerstand'], cond)
        area = xr.where(area_pol != 0, area+area_pol, area)
        bot = xr.where(area_pol != 0, peilen[row['peil']].mean() - diepte_pand_gem, bot)
        peil = xr.where(area_pol != 0, row['peil'], peil)
    
    model_ds_out = util.get_model_ds_empty(model_ds)
    model_ds_out[f'{name}_area'] = area
    model_ds_out[f'{name}_cond'] = cond
    model_ds_out[f'{name}_bot'] = bot

    # get riv peil data array 
    if (model_ds.gridtype == 'structured') and model_ds.steady_state:
        empty_time_array = np.zeros((model_ds.dims['y'],
                                     model_ds.dims['x']))*np.nan
        model_ds_out[f'{name}_peil'] = xr.DataArray(empty_time_array,
                                          dims=('y', 'x'),
                                          coords={'x': model_ds.x,
                                                  'y': model_ds.y})

    elif (model_ds.gridtype == 'structured') and (not model_ds.steady_state):
        empty_time_array = np.zeros((model_ds.dims['y'],
                                     model_ds.dims['x'],
                                     model_ds.dims['time']))*np.nan
        model_ds_out[f'{name}_peil'] = xr.DataArray(empty_time_array,
                                          dims=('y', 'x', 'time'),
                                          coords={'time': model_ds.time,
                                                  'x': model_ds.x,
                                                  'y': model_ds_out.y})

    elif (model_ds.gridtype == 'unstructured') and model_ds.steady_state:
        empty_time_array = np.zeros((model_ds.dims['cid']))*np.nan
        model_ds_out[f'{name}_peil'] = xr.DataArray(empty_time_array,
                                          dims=('cid'),
                                          coords={'cid': model_ds.cid})
    elif (model_ds.gridtype == 'unstructured') and (not model_ds.steady_state):
        empty_time_array = np.zeros((model_ds.dims['cid'],
                                     model_ds.dims['time']))*np.nan
        model_ds_out[f'{name}_peil'] = xr.DataArray(empty_time_array,
                                          dims=('cid', 'time'),
                                          coords={'time': model_ds.time,
                                                  'cid': model_ds.cid})

    # start en eindtijd  
    start_ts = pd.Timestamp(model_ds.time.data[0])
    end_ts = pd.Timestamp(model_ds.time.data[-1])
    if model_ds.steady_state or model_ds.steady_start:
        start = dt.datetime(start_ts.year - 1, start_ts.month, start_ts.day)
        end = end_ts + pd.Timedelta(1, unit='D')
    else:
        start = start_ts - pd.Timedelta(1, unit='D')
        end = end_ts + pd.Timedelta(1, unit='D')
    
    # voeg peil toe aan model dataset
    for col in peilen.columns:
        pandpeil_ts = peilen[col]
        pandpeil_average = pandpeil_ts.mean()
        if pandpeil_ts.index[-1] < end:
            if verbose:
                print(f'no pandpeil available after {pandpeil_ts.index[-1]} using average pandpeil NAP+{pandpeil_average:.2f} after this date')
            pandpeil_ts.loc[end] = pandpeil_average
            pandpeil_ts.loc[pandpeil_ts.index[-1]+dt.timedelta(days=1)] = pandpeil_average
            pandpeil_ts.sort_index(inplace=True)
        elif pandpeil_ts.index[0] > start:
            if verbose:
                print(f'no pandpeil available before {pandpeil_ts.index[0]} using average pandpeil NAP+{pandpeil_average:.2f} before this date')
            pandpeil_ts.loc[start] = pandpeil_average
            pandpeil_ts.loc[pandpeil_ts.index[0]-dt.timedelta(days=1)] = pandpeil_average
            pandpeil_ts.sort_index(inplace=True)
            
        if model_ds.steady_state:
            pandpeil_ts = pd.Series(index=[0], data=[pandpeil_average])
            model_ds_out[f'{name}_peil'] = xr.where(peil==col, pandpeil_average, model_ds_out[f'{name}_peil'])
        else:
            if model_ds.steady_start:
                pandpeil_ts = pandpeil_ts.reindex(model_ds.time.data, method='bfill')
                pandpeil_ts.loc[model_ds.time.data[0]] = pandpeil_average
            else:
                pandpeil_ts = pandpeil_ts.reindex(model_ds.time.data, method='bfill')

            # add data to model_ds_out
            if model_ds.gridtype == 'structured':
                allrow, allcol = np.where(peil==col)
                for row, col in zip(allrow, allcol):
                    model_ds_out[f'{name}_peil'].data[row, col, :] = pandpeil_ts.values
            elif model_ds.gridtype == 'unstructured':
                model_ds_out[f'{name}_peil'].loc[peil==col, :] = pandpeil_ts.values

    return model_ds_out

def model_dataset_to_panden_riv(model_ds, gwf, name='pand_riv'):
    """ add the infiltration ponds to the model via the riv package.
    

    Parameters
    ----------
    model_ds : xr.DataSet
        dataset containing relevant model information
    gwf : flopy.mf6.modflow.mfgwf.ModflowGwf
        groundwater flow model.
    name : str, optional
        name, used as variable name in model_ds. The default is 'panden'.
    
    Raises
    ------
    NotImplementedError
        function not yet implemented for structured grids.

    Returns
    -------
    riv_panden : flopy.mf6.ModflowGwfriv
        riv package.

    """
    if model_ds.gridtype!='unstructured':
        raise NotImplementedError('this function is not yet available for structured grids')
        
    if model_ds.steady_state:
        riv_spd = mgrid.data_array_unstructured_to_rec_list(model_ds, 
                                                            model_ds[f'{name}_area']>0,
                                                            col1=f'{name}_peil',
                                                            col2=f'{name}_cond',
                                                            col3=f'{name}_bot',
                                                            first_active_layer=True,
                                                            only_active_cells=False)
        
        riv_panden = flopy.mf6.ModflowGwfriv(gwf, 
                                             filename=f'{gwf.name}_panden.riv',
                                             pname=f'riv_{name}',
                                             maxbound=len(riv_spd),
                                             print_input=True,
                                             stress_period_data={0: riv_spd})
    
    # transient pandpeilen
    
    
    # vind unieke reeksen met pandpeilen
    pp_unique_arr = np.unique(model_ds[f'{name}_peil'].data, axis=0)
    pandpeilen = pp_unique_arr[~np.isnan(pp_unique_arr).all(axis=1)]
    
    # maak raster met unieke naam per unieke reeks
    empty_str_array = np.zeros_like(model_ds['idomain'][0], dtype="S13")
    model_ds[f'{name}_name'] = xr.DataArray(empty_str_array,
                                            dims=('cid'),
                                            coords={'cid': model_ds.cid})
    model_ds[f'{name}_name'] = model_ds[f'{name}_name'].astype(str)
    pp_unique_dic = {}
    for i, pandpeil in enumerate(pandpeilen):
        model_ds[f'{name}_name'][(model_ds[f'{name}_peil'].data == pandpeil).all(
            axis=1)] = f'{name}_{i}'
        
        pp_unique_dic[f'{name}_{i}'] = pandpeil
    
    # create riv package
    mask = model_ds[f'{name}_name'] != ''
    riv_spd = mgrid.data_array_unstructured_to_rec_list(model_ds, mask,
                                                        col1=f'{name}_name',
                                                        col2=f'{name}_cond',
                                                        col3=f'{name}_bot',
                                                        first_active_layer=True,
                                                        only_active_cells=False)
    
    riv_panden = flopy.mf6.ModflowGwfriv(gwf, 
                                         filename=f'{gwf.name}_panden.riv',
                                         pname=f'riv_{name}',
                                         maxbound=len(riv_spd),
                                         print_input=True,
                                         stress_period_data={0: riv_spd})
    
    # create timeseries packages
    for i, key in enumerate(pp_unique_dic.keys()):
        riv_peil_val = list(pp_unique_dic[key]) + [0.0]
        time_steps_riv = list(model_ds['time_steps'].data) + [model_ds.nper]
        riv_peil = list(zip(time_steps_riv, riv_peil_val))
        if i == 0:
            riv_panden.ts.initialize(filename=f'{key}.ts',
                                     timeseries=riv_peil,
                                     time_series_namerecord=key,
                                     interpolation_methodrecord='stepwise')
        else:
            riv_panden.ts.append_package(filename=f'{key}.ts',
                                         timeseries=riv_peil,
                                         time_series_namerecord=key,
                                         interpolation_methodrecord='stepwise')    
            
    return riv_panden
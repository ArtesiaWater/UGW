# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 11:16:45 2020

@author: oebbe
"""

import os

import pandas as pd
import xarray as xr
import numpy as np
from . import util, mgrid

def get_layer_models(datadir, extent, delr, delc,
                     regis=True, geotop=True, pwn_model=False,
                     cachedir=None,
                     fname_netcdf='combined_layer_ds.nc',
                     use_cache=False, verbose=False):
    
    combined_ds = util.get_cache_netcdf(use_cache, cachedir, fname_netcdf,
                                        get_combined_layer_models,
                                        verbose=verbose,
                                        datadir=datadir, extent=extent,
                                        delr=delr, delc=delc,
                                        regis=True, geotop=True,
                                        pwn_model=pwn_model)

    return combined_ds


def get_combined_layer_models(datadir, extent, delr, delc, 
                              regis=True, geotop=True, pwn_model=False,
                              cachedir=None, use_cache=False, verbose=True):
    """ combine layer models into a single layer model. 
    
    Possibilities so far inlcude:
        - regis -> full model based on regis
        - regis and geotop -> holoceen of REGIS is filled with geotop
        - regis and pwn -> use pwn model where available, below and on the 
        sides use regis
        - regis, geotop and pwn -> use pwn model where available, below and on 
        the sides use the regis/geotop combination.
    

    Parameters
    ----------
    datadir : str, optional
        datadirectory met koppeltabel. The default is None.
    extent : list, tuple or np.array
        desired model extent (xmin, xmax, ymin, ymax)
    delr : int or float,
        cell size along rows, equal to dx
    delc : int or float,
        cell size along columns, equal to dy
    regis : bool, optional
        True if part of the layer model should be REGIS. The default is True.
    geotop : bool, optional
        True if part of the layer model should be geotop. The default is True.
    pwn_model : bool, optional
        True if part of the layer model should be pwn. The default is False.
    cachedir : str
        directory to store cached values, if None a temporary directory is
        used. default is None
    use_cache : bool, optional
        if True the cached resampled regis dataset is used.
        The default is False.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        if an invalid combination of layers is used.

    Returns
    -------
    combined_ds : xarray dataset
        combination of layer models.

    """
    
    if regis:
        # get local regis dataset
        regis_path = os.path.join(datadir,'regis_nhflopy.nc')
        regis_ds_raw = xr.open_dataset(regis_path).sel(x=slice(extent[0], extent[1]),
                                                       y=slice(extent[2], extent[3]))
        
        # convert regis dataset to grid
        regis_ds = util.get_ml_layer_dataset_struc(raw_ds=regis_ds_raw, extent=extent, 
                                                   delr=delr, delc=delc, 
                                                   cachedir=cachedir, 
                                                   fname_netcdf='regis_ds.nc',
                                                   use_cache=use_cache, 
                                                   verbose=verbose)
    else:
        raise ValueError('layer models without REGIS not supported')
    
    if geotop:
        geotop_path = os.path.join(datadir,'geotop_nhflopy.nc')
        litho_translate_df = pd.read_csv(os.path.join(datadir, 'geotop','litho_eenheden.csv'), index_col=0)
        geo_eenheid_translate_df = pd.read_csv(os.path.join(datadir, 'geotop','geo_eenheden.csv'), index_col=0,
                                               keep_default_na=False)
        
        geotop_ds_raw = get_geotop_raw(geotop_path=geotop_path, 
                                        regis_ds=regis_ds_raw, 
                                        regis_layer='HLc',
                                        litho_translate_df=litho_translate_df,
                                        geo_eenheid_translate_df=geo_eenheid_translate_df,
                                        cachedir=cachedir, 
                                        fname_netcdf='geotop_raw.nc',
                                        use_cache=use_cache,
                                        verbose=verbose)
        
        geotop_ds = util.get_ml_layer_dataset_struc(raw_ds=geotop_ds_raw, 
                                                    extent=extent, 
                                                    delr=delr, delc=delc, 
                                                    cachedir=cachedir, 
                                                    fname_netcdf='geotop_ds.nc',
                                                    use_cache=use_cache, 
                                                    verbose=verbose)
    
    if regis and geotop:
        regis_geotop_ds = get_regis_geotop_ds(regis_ds, geotop_ds,
                                              cachedir=cachedir, 
                                              fname_netcdf='regis_geotop_ds.nc',
                                              use_cache=use_cache, 
                                              verbose=verbose)
        
    if pwn_model:
        model_ws_pwn = r'c:\Users\oebbe\02_python\NHFLO\NHFLO\work\pwn_model'

        # load pwn model dataset
        model_ds_pwn = xr.open_dataset(os.path.join(model_ws_pwn, 'full_model_ds.nc'))
        model_ds_pwn.attrs['model_name'] = 'pwn_model'
        
        # find overlap between models -> extent_pwn
        extent_pwn = get_pwn_extent(regis_ds, model_ds_pwn.extent, 
                                    verbose=verbose)
        
        # convert pwn dataset to grid
        pwn_ds = util.get_ml_layer_dataset_struc(raw_ds=model_ds_pwn[['bot','top','kh','kv']], 
                                                 extent=list(extent_pwn), 
                                                 delr=delr, delc=delc, 
                                                 cachedir=cachedir, 
                                                 fname_netcdf='pwn_ds.nc',
                                                 use_cache=use_cache, 
                                                 verbose=verbose)
        
    if regis and geotop and pwn_model:
        
        
        fname_koppeltabel = os.path.join(datadir,'pwn_modellagen', 
                                         'combine_regis_gtop_pwn.xlsx')
        df_koppeltabel = pd.read_excel(fname_koppeltabel, skiprows=1, 
                                       index_col=4)
        
        
        
        combined_ds = get_pwn_regis_ds(pwn_ds, regis_geotop_ds, datadir,
                                       cachedir=cachedir,
                                       fname_netcdf='pwn_regis_geotop_ds.nc',
                                       df_koppeltabel=df_koppeltabel,
                                       use_cache=use_cache,
                                       verbose=verbose)
        
        
        
        
    elif regis and pwn_model:
        
        fname_koppeltabel = os.path.join(datadir,'pwn_modellagen', 
                                         'combine_regis_pwn.xlsx')
        df_koppeltabel = pd.read_excel(fname_koppeltabel, skiprows=1, 
                                       index_col=4)
        
        combined_ds = get_pwn_regis_ds(pwn_ds, regis_ds, datadir,
                                        cachedir=cachedir,
                                        fname_netcdf='pwn_regis_ds.nc',
                                        df_koppeltabel=df_koppeltabel,
                                        use_cache=False,
                                        verbose=verbose)
    
    elif regis and geotop:
        combined_ds = regis_geotop_ds
        
    elif regis:
        combined_ds = regis_ds
        
        
    else:
        raise ValueError('combination of model layers not supported')

    
    return combined_ds


def get_geotop_raw(geotop_path, regis_ds=None, regis_layer=None,
                   litho_translate_df=None, 
                   geo_eenheid_translate_df=None,
                   cachedir=None,
                   fname_netcdf='geotop_raw.nc',
                   use_cache=False,
                   verbose=False):
    """ get geotop raw layer model
    """
    
    geotop_ds_raw = util.get_cache_netcdf(use_cache, cachedir, fname_netcdf,
                                          convert_geotop_to_ml_layers,
                                          verbose=verbose,
                                          geotop_path=geotop_path, 
                                          regis_ds=regis_ds,
                                          regis_layer=regis_layer,
                                          litho_translate_df=litho_translate_df,
                                          geo_eenheid_translate_df=geo_eenheid_translate_df,
                                          check_time=False)
    
    return geotop_ds_raw


def convert_geotop_to_ml_layers(geotop_path, regis_ds=None, regis_layer=None,
                                litho_translate_df=None, 
                                geo_eenheid_translate_df=None,
                                verbose=False):
    """ does the following steps to obtain model layers based on geotop:
        1. get geotop from url (or local file)
        2. slice by regis layer (if not None)
        3. compute kh from lithoklasse
        4. create a layer model based on geo-eenheden
    
    Parameters
    ----------
    geotop_path: str
        path to geotop netcdf, can be weblocation or preprocessed netcdf file
    regis_ds: xarray.DataSet
        regis dataset used to cut geotop to the same x and y coördinates    
    regis_layer: str, optional
        layer of regis dataset that will be filled with geotop 
    litho_translate_df: pandas.DataFrame
        horizontal conductance (kh)
    geo_eenheid_translate_df: pandas.DataFrame
        dictionary to translate geo_eenheid to a geo name    
    verbose : bool, optional
        print additional information. default is False
    
    Returns
    -------
    geotop_ds_raw: xarray.DataSet
        geotop dataset with added horizontal conductance

    """
    
    # stap 1 and 2
    if (regis_ds is not None) and (regis_layer is not None):
        if verbose:
            print(f'slice geotop with regis layer {regis_layer}')
        top_rl = regis_ds.top.sel(layer=regis_layer)
        bot_rl = regis_ds.bottom.sel(layer=regis_layer)
    
        geotop_ds_raw = xr.open_dataset(geotop_path).sel(x=regis_ds.x,
                                                         y=regis_ds.y,
                                                         z=slice(np.floor(bot_rl.min().data),
                                                                 np.ceil(top_rl.max().data)))
    else:
        geotop_ds_raw = xr.open_dataset(geotop_path)
        
    # stap 3 maak kh matrix a.d.v. lithoklasse
    if verbose:
        print('create kh matrix from lithoklasse and csv file')
    kh_from_litho = xr.zeros_like(geotop_ds_raw.lithok)
    for i, row in litho_translate_df.iterrows():
        kh_from_litho = xr.where(geotop_ds_raw.lithok==i, 
                                 row['hor_conductance_default'], 
                                 kh_from_litho)
    geotop_ds_raw['kh_from_litho'] = kh_from_litho
    
    # stap 4 maak een laag per geo-eenheid
    geotop_ds_mod = get_top_bot_from_geo_eenheid(geotop_ds_raw, geo_eenheid_translate_df,
                                                 verbose=verbose)
        
    return geotop_ds_mod


def get_top_bot_from_geo_eenheid(geotop_ds_raw, geo_eenheid_translate_df,
                                 verbose=False):
    """ get top, bottom and kh of each geo-eenheid in geotop dataset.
    
    Parameters
    ----------
    geotop_ds_raw: xr.DataSet
        geotop dataset with added horizontal conductance
    geo_eenheid_translate_df: pandas.DataFrame
        dictionary to translate geo_eenheid to a geo name
    verbose : bool, optional
        print additional information. default is False
        
    Returns
    -------
    geotop_ds_mod: xr.DataSet
        geotop dataset with top, bot, kh and kv per geo_eenheid
    
    Note
    ----
    the 'geo_eenheid' >6000 are 'stroombanen' these are difficult to add because
    they can occur above and below any other 'geo_eenheid' therefore they are
    added to the geo_eenheid below the stroombaan.
    
    """
    

    # vindt alle geo-eenheden in model_extent
    geo_eenheden = np.unique(geotop_ds_raw.strat.data)
    geo_eenheden = geo_eenheden[np.isfinite(geo_eenheden)]
    stroombaan_eenheden = geo_eenheden[geo_eenheden<5999]
    geo_eenheden = geo_eenheden[geo_eenheden<5999]

    # geo eenheid 2000 zit boven 1130
    if (2000. in geo_eenheden) and (1130. in geo_eenheden):
        geo_eenheden[(geo_eenheden==2000.) + (geo_eenheden==1130.)]= [2000., 1130.]

    geo_names = [geo_eenheid_translate_df.loc[geo_eenh, 'Code (lagenmodel en boringen)'] for geo_eenh in geo_eenheden]

    # fill top and bot
    top = np.ones((geotop_ds_raw.y.shape[0], geotop_ds_raw.x.shape[0], len(geo_names))) * np.nan
    bot = np.ones((geotop_ds_raw.y.shape[0], geotop_ds_raw.x.shape[0], len(geo_names))) * np.nan
    lay = 0
    if verbose:
        print('creating top and bot per geo eenheid')
    for geo_eenheid in geo_eenheden:
        if verbose:
            print(geo_eenheid)
        
        mask = geotop_ds_raw.strat == geo_eenheid
        geo_z = xr.where(mask, geotop_ds_raw.z, np.nan)
       
        top[:,:,lay] = geo_z.max(dim='z').T+0.5
        bot[:,:,lay] = geo_z.min(dim='z').T
        
        lay+=1
        
    geotop_ds_mod = add_stroombanen_and_get_kh(geotop_ds_raw, top, bot, 
                                               geo_names,
                                               verbose=verbose)
    
    geotop_ds_mod.attrs['stroombanen'] = stroombaan_eenheden
    
    return geotop_ds_mod

def add_stroombanen_and_get_kh(geotop_ds_raw, top, bot, geo_names, verbose=False):
    """ add stroombanen to tops and bots of geo_eenheden, also computes kh per 
    geo_eenheid. Kh is computed by taking the average of all kh's of a geo_eenheid
    within a cell (e.g. if one geo_eenheid has a thickness of 1,5m in a certain
    cell the kh of the cell is calculated as the mean of the 3 cells in geotop)
    
    Parameters
    ----------
    geotop_ds_raw: xr.DataSet
        geotop dataset with added horizontal conductance
    top: np.array
        raster with top of each geo_eenheid, shape(nlay,nrow,ncol)
    bot: np.array
        raster with bottom of each geo_eenheid, shape(nlay,nrow,ncol)
    geo_names: list of str
        names of each geo_eenheid
    verbose : bool, optional
        print additional information. default is False
        
    Returns
    -------
    geotop_ds_mod: xr.DataSet
        geotop dataset with top, bot, kh and kv per geo_eenheid

    """
    kh = np.ones((geotop_ds_raw.y.shape[0], geotop_ds_raw.x.shape[0], len(geo_names))) * np.nan
    thickness = np.ones((geotop_ds_raw.y.shape[0], geotop_ds_raw.x.shape[0], len(geo_names))) * np.nan
    z = xr.ones_like(geotop_ds_raw.lithok)*geotop_ds_raw.z
    if verbose:
        print(f'adding stroombanen to top and bot of each layer')
        print(f'get kh for each layer')
        
    for lay in range(top.shape[2]):
        if verbose:
            print(geo_names[lay])
        if lay==0:
            top[:,:,0] = np.nanmax(top, axis=2)
        else:
            top[:,:,lay] = bot[:,:,lay-1]
        bot[:,:,lay] = np.where(np.isnan(bot[:,:,lay]), top[:,:,lay], bot[:,:,lay])
        thickness[:,:,lay] = (top[:,:,lay] - bot[:,:,lay])
        
        #check which geotop voxels are within the range of the layer
        bool_z = xr.zeros_like(z)
        for i in range(z.z.shape[0]):
            bool_z[:,:,i] = np.where((z[:,:,i]>=bot[:,:,lay].T) * (z[:,:,i]<top[:,:,lay].T), True, False)
            
        kh_geo = xr.where(bool_z, geotop_ds_raw['kh_from_litho'], np.nan)
        kh[:,:,lay] = kh_geo.mean(dim='z').T
    
    
    da_top = xr.DataArray(data=top, dims=('y', 'x', 'layer'),
                          coords={'y': geotop_ds_raw.y,'x': geotop_ds_raw.x,
                                  'layer': geo_names})
    da_bot = xr.DataArray(data=bot, dims=('y', 'x', 'layer'),
                          coords={'y': geotop_ds_raw.y,'x': geotop_ds_raw.x,
                                  'layer': geo_names})
    da_kh = xr.DataArray(data=kh, dims=('y', 'x', 'layer'),
                         coords={'y': geotop_ds_raw.y,'x': geotop_ds_raw.x,
                                 'layer': geo_names})
    da_thick = xr.DataArray(data=thickness, dims=('y', 'x', 'layer'),
                            coords={'y': geotop_ds_raw.y,'x': geotop_ds_raw.x,
                                    'layer': geo_names})
    
    geotop_ds_mod = xr.Dataset()   
    
    geotop_ds_mod['top'] = da_top
    geotop_ds_mod['bottom'] = da_bot
    geotop_ds_mod['kh']  = da_kh
    geotop_ds_mod['kv']  = geotop_ds_mod['kh'] * .25
    geotop_ds_mod['thickness'] = da_thick
    
    return geotop_ds_mod


def get_regis_geotop_ds(regis_ds, geotop_ds, cachedir=None,
                        fname_netcdf='regis_geotop.nc',
                        use_cache=False,
                        verbose=False):
    """ Get a combination of regis and geotop for the model layers.
    
    Parameters
    ----------
    regis_ds: xarray.DataSet
        regis dataset
    geotop_ds: xarray.DataSet
        regis dataset
    verbose : bool, optional
        print additional information. default is False
    
    Returns
    -------
    regis_geotop_ds: xr.DataSet
        combined dataset  
    """
    
    
    
    regis_geotop_ds = util.get_cache_netcdf(use_cache, cachedir, fname_netcdf,
                                            add_geotop_to_regis_hlc,
                                            verbose=verbose,
                                            regis_ds=regis_ds,
                                            geotop_ds=geotop_ds,
                                            check_time=False)
    
    return regis_geotop_ds

def add_geotop_to_regis_hlc(regis_ds, geotop_ds, float_correctie=0.001,
                            verbose=False):
    """ Combine geotop and regis in such a way that the holoceen in Regis is
    replaced by the geo_eenheden of geotop.
    
    Parameters
    ----------
    regis_ds: xarray.DataSet
        regis dataset
    geotop_ds: xarray.DataSet
        regis dataset
    float_correctie: float
        due to floating point precision some floating point numbers that are
        the same are not recognised as the same. Therefore this correction is
        used.
    verbose : bool, optional
        print additional information. default is False
    
    Returns
    -------
    regis_geotop_ds: xr.DataSet
        combined dataset  
    
    
    """
    regis_geotop_ds = xr.Dataset()
    
    new_layers = np.append(geotop_ds.layer.data, regis_ds.layer.data[1:].astype('<U8')).astype('O')
    
    top = xr.DataArray(dims=('layer', 'y', 'x'),
                       coords={'y': geotop_ds.y,'x': geotop_ds.x,
                               'layer': new_layers})
    bot = xr.DataArray(dims=('layer', 'y', 'x'),
                       coords={'y': geotop_ds.y,'x': geotop_ds.x,
                               'layer': new_layers})
    kh = xr.DataArray(dims=('layer', 'y', 'x'),
                      coords={'y': geotop_ds.y,'x': geotop_ds.x,
                              'layer': new_layers})
    kv = xr.DataArray(dims=('layer', 'y', 'x'),
                      coords={'y': geotop_ds.y,'x': geotop_ds.x,
                              'layer': new_layers})
    
    # haal overlap tussen geotop en regis weg
    if verbose:
        print('cut geotop layer based on regis holoceen')
    for lay in range(geotop_ds.dims['layer']):
        # Alle geotop cellen die onder de onderkant van het holoceen liggen worden inactief
        mask1 = geotop_ds['top'][lay]<=(regis_ds['bottom'][0] - float_correctie)
        geotop_ds['top'][lay] = xr.where(mask1, np.nan, geotop_ds['top'][lay])
        geotop_ds['bottom'][lay] = xr.where(mask1, np.nan, geotop_ds['bottom'][lay])
        
        # Alle geotop cellen waarvan de bodem onder de onderkant van het holoceen ligt, krijgen als bodem de onderkant van het holoceen
        mask2 = geotop_ds['bottom'][lay]<regis_ds['bottom'][0]
        geotop_ds['bottom'][lay] = xr.where(mask2 * (~mask1), regis_ds['bottom'][0], geotop_ds['bottom'][lay])
        
        # Alle geotop cellen die boven de bovenkant van het holoceen liggen worden inactief
        mask3 = geotop_ds['bottom'][lay]>=(regis_ds['top'][0] - float_correctie)
        geotop_ds['top'][lay] = xr.where(mask3, np.nan, geotop_ds['top'][lay])
        geotop_ds['bottom'][lay] = xr.where(mask3, np.nan, geotop_ds['bottom'][lay])
        
        # Alle geotop cellen waarvan de top boven de top van het holoceen ligt, krijgen als top het holoceen van regis
        mask4 = geotop_ds['top'][lay]>=regis_ds['top'][0]
        geotop_ds['top'][lay] = xr.where(mask4 * (~mask3), regis_ds['top'][0], geotop_ds['top'][lay])
        
        # overal waar holoceen inactief is, wordt geotop ook inactief
        mask5 = regis_ds['bottom'][0].isnull()
        geotop_ds['top'][lay] = xr.where(mask5, np.nan, geotop_ds['top'][lay])
        geotop_ds['bottom'][lay] = xr.where(mask5, np.nan, geotop_ds['bottom'][lay])
        if verbose:
            if (mask2 * (~mask1)).sum()>0:
                print(f'regis holoceen snijdt door laag {geotop_ds.layer[lay].values}')
        
    
    top[:len(geotop_ds.layer),:,:] = geotop_ds['top'].data
    top[len(geotop_ds.layer):,:,:] = regis_ds['top'].data[1:]
    
    bot[:len(geotop_ds.layer),:,:] = geotop_ds['bottom'].data
    bot[len(geotop_ds.layer):,:,:] = regis_ds['bottom'].data[1:]  
    
    kh[:len(geotop_ds.layer),:,:] = geotop_ds['kh'].data
    kh[len(geotop_ds.layer):,:,:] = regis_ds['kh'].data[1:]  
    
    kv[:len(geotop_ds.layer),:,:] = geotop_ds['kv'].data
    kv[len(geotop_ds.layer):,:,:] = regis_ds['kv'].data[1:]  
        
    regis_geotop_ds['top'] = top
    regis_geotop_ds['bottom'] = bot
    regis_geotop_ds['kh']  = kh
    regis_geotop_ds['kv']  = kv
    
    _ = [regis_geotop_ds.attrs.update({key: item})
                 for key, item in regis_ds.attrs.items()]
    
    #maak bottom nan waar de laagdikte 0 is
    regis_geotop_ds['bottom'] = xr.where((regis_geotop_ds['top']-regis_geotop_ds['bottom'])<float_correctie, 
                                          np.nan,
                                          regis_geotop_ds['bottom'])
    
    return regis_geotop_ds

def get_pwn_regis_ds(pwn_ds, regis_ds, datadir=None,
                     df_koppeltabel=None,
                     cachedir=None,
                     fname_netcdf='pwn_regis_ds.nc',
                     use_cache=False,
                     verbose=False):
    """
    

    Parameters
    ----------
    pwn_ds : xr.DataSet
        lagenmodel van pwn.
    regis_ds : xr.DataSet
        lagenmodel regis.
    cachedir : TYPE, optional
        DESCRIPTION. The default is None.
    fname_netcdf : TYPE, optional
        DESCRIPTION. The default is 'pwn_regis_ds.nc'.
    use_cache : TYPE, optional
        DESCRIPTION. The default is False.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    pwn_regis_ds : xr.DataSet
        combined model

    """
    
    pwn_regis_ds = util.get_cache_netcdf(use_cache, cachedir, fname_netcdf,
                                         create_pwn_regis_ds,
                                         verbose=verbose,
                                         check_time=False, 
                                         pwn_ds=pwn_ds, regis_ds=regis_ds,
                                         datadir=datadir,
                                         df_koppeltabel=df_koppeltabel)
    
    return pwn_regis_ds

def create_pwn_regis_ds(pwn_ds, regis_ds, datadir=None,
                        df_koppeltabel=None, verbose=True):
    """ Create a new layer model based on regis and pwn models.
    

    Parameters
    ----------
    pwn_ds : xr.DataSet
        lagenmodel van pwn.
    regis_ds : xr.DataSet
        lagenmodel regis.
    datadir : str, optional
        datadirectory met koppeltabel. The default is None.
    df_koppeltabel : pandas DataFrame, optional
        dataframe van koppeltabel. The default is None.
    verbose : TYPE
        DESCRIPTION.

    Raises
    ------
    NotImplementedError
        some combinations of regis and pwn are not implemented yet.

    Returns
    -------
    pwn_regis_ds : xr.DataSet
        combined model

    """
    
    if util.compare_model_extents(regis_ds.extent, pwn_ds.extent, verbose)==1:
        pwn_regis_ds = add_regis_to_bottom_of_pwn(pwn_ds, regis_ds, verbose)
    else:
        pwn_regis_ds = combine_layer_models_regis_pwn(pwn_ds, regis_ds, 
                                                      datadir, df_koppeltabel,
                                                      verbose)
        
    return pwn_regis_ds

def combine_layer_models_regis_pwn(pwn_ds, regis_ds, datadir=None,
                                   df_koppeltabel=None, verbose=False):
    """ combine model layers from regis and pwn using a 'koppeltabel'
    

    Parameters
    ----------
    pwn_ds : xr.DataSet
        lagenmodel van pwn.
    regis_ds : xr.DataSet
        lagenmodel regis.
    datadir : str, optional
        datadirectory met koppeltabel. The default is None.
    df_koppeltabel : pandas DataFrame, optional
        dataframe van koppeltabel. The default is None.
    verbose : bool, optional
        

    Raises
    ------
    ValueError
        invalid values in koppeltabel.

    Returns
    -------
    pwn_regis_ds : xr.DataSet
        combined model
    """

    if df_koppeltabel is None:
        fname_koppeltabel = os.path.join(datadir,'pwn_modellagen', 
                                         'combine_regis_pwn.xlsx')
        df_koppeltabel = pd.read_excel(fname_koppeltabel, skiprows=1, 
                                       index_col=4)

    pwn_regis_ds = xr.Dataset(coords={'x': regis_ds.x.data,
                                      'y': regis_ds.y.data,
                                      'layer': df_koppeltabel.index.values})
    
    _ = [pwn_regis_ds.attrs.update({key: item})
         for key, item in regis_ds.attrs.items()]
    
    empty_da = xr.DataArray(dims=('layer', 'y', 'x'),
                            coords={'x': regis_ds.x.data,
                                    'y': regis_ds.y.data,
                                    'layer': df_koppeltabel.index.values})
    
    
    bot = empty_da.copy()
    top = empty_da.copy()
    kh = empty_da.copy()
    kv = empty_da.copy()
    
    y_mask = [True  if y in pwn_ds.y else False for y in regis_ds.y.values]
    x_mask = [True  if x in pwn_ds.x else False for x in regis_ds.x.values]
    
    column_reg_mod = df_koppeltabel.columns[1]
    for i, lay in enumerate(df_koppeltabel.index):
        regis_lay = df_koppeltabel.loc[lay, column_reg_mod]
        pwn_lay   = df_koppeltabel.loc[lay, 'pwn_lay']
        
        if verbose:
            print(f'combine regis layer {regis_lay} with pwn layer {pwn_lay}')
            
        
        if regis_lay==pwn_lay:
            raise ValueError(f'invalid values encountered, regis layer is {regis_lay} and pwn layer is {pwn_lay}')
        
        if isinstance(regis_lay, str):
            bot[i] = regis_ds.bottom.sel(layer=regis_lay)
            top[i] = regis_ds.top.sel(layer=regis_lay)
            kh[i] = regis_ds.kh.sel(layer=regis_lay)
            kv[i] = regis_ds.kv.sel(layer=regis_lay)
        elif np.isnan(regis_lay):
            bot[i] = np.nan
            top[i] = np.nan
            kh[i] = np.nan
            kv[i] = np.nan
        else:
            raise ValueError('invalid value encountered in regis_lay_nam')
        
        
        if isinstance(pwn_lay,int):
            # brand pwn modellaag in regis laag
            bot[i,y_mask,x_mask] = pwn_ds.bot.sel(layer=pwn_lay-1)
            top[i,y_mask,x_mask] = pwn_ds.top.sel(layer=pwn_lay-1)
            kh[i,y_mask,x_mask] = pwn_ds.kh.sel(layer=pwn_lay-1)
            kv[i,y_mask,x_mask] = pwn_ds.kv.sel(layer=pwn_lay-1)
            pwn_final_lay = pwn_lay
        elif pwn_lay =='REGIS':
            # plak REGIS model onder pwn model
            regis_bot = bot[i,y_mask, x_mask]
            regis_top = top[i,y_mask,x_mask]
            regis_kh  = kh[i,y_mask,x_mask]
            regis_kv  = kv[i,y_mask,x_mask]
            pwn_bot = pwn_ds.bot.sel(layer=pwn_final_lay-1)
            
            bot[i,y_mask, x_mask] = xr.where(regis_bot<pwn_bot, regis_bot, np.nan)
            top[i,y_mask, x_mask] = xr.where(regis_top<pwn_bot, regis_top, pwn_bot)
            kh[i,y_mask,x_mask] = xr.where(regis_bot<pwn_bot, regis_kh, np.nan)
            kv[i,y_mask,x_mask] = xr.where(regis_bot<pwn_bot, regis_kv, np.nan)
        elif np.isnan(pwn_lay):
            # maak laag met idomain -1 waar wel regis maar geen pwn laag zit
            bot[i,y_mask,x_mask] = np.nan
            top[i,y_mask,x_mask] = np.nan
            kh[i,y_mask,x_mask] = np.nan
            kv[i,y_mask,x_mask] = np.nan
        else:
            raise ValueError('invalid value encountered in pwn_lay')
    
    pwn_regis_ds['bottom'] = bot
    pwn_regis_ds['top'] = top
    pwn_regis_ds['kh'] = kh
    pwn_regis_ds['kv'] = kv
    
    return pwn_regis_ds

    
def add_regis_to_bottom_of_pwn(pwn_ds, regis_ds, verbose=False):
    """ extend the pwn model by using the regis model for the layers below
    the pwn model.
    

    Parameters
    ----------
    pwn_ds : xr.DataSet
        lagenmodel van pwn.
    regis_ds : xr.DataSet
        lagenmodel regis.

    Returns
    -------
    pwn_regis_ds : xr.DataSet
        combined model

    """
    lay_count = len(pwn_ds.layer)
    new_bot = pwn_ds['bot'].data.copy()
    new_top = pwn_ds['top'].data.copy()
    new_kh  = pwn_ds['kh'].data.copy()
    new_kv  = pwn_ds['kv'].data.copy()
    new_layer = pwn_ds['bot'].layer.data.copy()
    lname_regis = []
    for lay in regis_ds.layer:
        mask_lay = pwn_ds['bot'][-1].data>regis_ds['bottom'].sel(layer=lay).data
        if mask_lay.any():
            bot_lay = np.where(mask_lay,
                               regis_ds['bottom'].sel(layer=lay).data,
                               np.nan)
            kh_lay = np.where(mask_lay,
                              regis_ds['kh'].sel(layer=lay).data,
                              np.nan)
            kv_lay = np.where(mask_lay,
                              regis_ds['kv'].sel(layer=lay).data,
                              np.nan)
            top_lay = new_bot[lay_count-1]
            new_bot = np.concatenate((new_bot, np.array([bot_lay])))
            new_kh  = np.concatenate((new_kh, np.array([kh_lay])))
            new_kv  = np.concatenate((new_kv, np.array([kv_lay])))
            new_top = np.concatenate((new_top, np.array([top_lay])))
            lname_regis.append(str(lay.values))
            new_layer = np.append(new_layer, lay_count)
            lay_count+=1
            if verbose:
                print(f'adding regis layer {str(lay.values)}  to pwn_model layers' )
         
            
    pwn_regis_ds = xr.Dataset(coords={'x':pwn_ds.x, 
                                      'y':pwn_ds.y, 
                                      'layer':new_layer})
        
    pwn_regis_ds['bottom'] = xr.DataArray(new_bot, dims=('layer', 'y', 'x'),
                                       coords={'layer':new_layer,
                                               'x':pwn_ds.x, 
                                               'y':pwn_ds.y})
    pwn_regis_ds['top'] = xr.DataArray(new_top, dims=('layer', 'y', 'x'),
                                       coords={'layer':new_layer,
                                               'x':pwn_ds.x, 
                                               'y':pwn_ds.y})
    pwn_regis_ds['kh'] = xr.DataArray(new_kh, dims=('layer', 'y', 'x'),
                                      coords={'layer':new_layer,
                                               'x':pwn_ds.x, 
                                               'y':pwn_ds.y})
    pwn_regis_ds['kv'] = xr.DataArray(new_kv, dims=('layer', 'y', 'x'),
                                      coords={'layer':new_layer,
                                               'x':pwn_ds.x, 
                                               'y':pwn_ds.y})
    
    lname_pwn = [f'pwn_lay_{i+1}' for i in range(len(pwn_ds.layer))]           
    lnames_pwn_regis = lname_pwn + lname_regis
    
    pwn_regis_ds['lnames'] = xr.DataArray(lnames_pwn_regis, dims=('layer'),
                                          coords={'layer':new_layer})           
    
    _ = [pwn_regis_ds.attrs.update({key: item}) for key, item in regis_ds.attrs.items()]
    
    return pwn_regis_ds

def get_pwn_extent(regis_ds, pwn_extent_original,
                   verbose=False):
    """ get the extent of the part of the pwn model that is inside the 
    regis model.
    

    Parameters
    ----------
    regis_ds : xarray dataset
        discretization of regis model.
    pwn_extent_original : list, tuple or numpy array
        extent of the original pwn model [xmin, xmax, ymin, ymax]

    Returns
    -------
    extent_pwn : list, tuple or numpy array
        extent of the part of the pwn model that is inside the regis model

    """
    
    model_layer_combi_type = util.compare_model_extents(regis_ds.extent, 
                                                        pwn_extent_original, 
                                                        verbose=verbose)
    
    delr = regis_ds.delr
    delc = regis_ds.delc
    
    x = regis_ds.x.values
    y = regis_ds.y.values
    
    if model_layer_combi_type==1:
        new_pwn_extent = regis_ds.extent.copy()
        
    else:
        xmin = x[x>=(pwn_extent_original[0]+0.5*delr)].min() - 0.5*delr
        xmax = x[x<=(pwn_extent_original[1]-0.5*delr)].max() + 0.5*delr
        ymin = y[y>=(pwn_extent_original[2]+0.5*delc)].min() - 0.5*delc
        ymax = y[y<=(pwn_extent_original[3]-0.5*delc)].max() + 0.5*delc
        new_pwn_extent = [xmin, xmax, ymin, ymax]
        
    return new_pwn_extent

def add_kh_kv_from_ml_layer_to_dataset(ml_layer_ds, model_ds, anisotropy,
                                       fill_value_kh, fill_value_kv,
                                       verbose=False):
    """ add kh and kv from a model layer dataset to THE model dataset.

    Supports structured and unstructured grids.

    Parameters
    ----------
    ml_layer_ds : xarray.Dataset
        dataset with model layer data with kh and kv
    model_ds : xarray.Dataset
        dataset with model data where kh and kv are added to
    anisotropy : int or float
        factor to calculate kv from kh or the other way around
    fill_value_kh : int or float, optional
        use this value for kh if there is no data in regis. The default is 1.0.
    fill_value_kv : int or float, optional
        use this value for kv if there is no data in regis. The default is 1.0.
    verbose : bool, optional
        print additional information. default is False

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data with new kh and kv

    Notes
    -----
    some model dataset, such as regis, also have 'c' and 'kd' values. These 
    are ignored at the moment
    """
    model_ds.attrs['anisotropy'] = anisotropy
    model_ds.attrs['fill_value_kh'] = fill_value_kh
    model_ds.attrs['fill_value_kv'] = fill_value_kv
    kh_arr = ml_layer_ds['kh'].data
    kv_arr = ml_layer_ds['kv'].data

    if verbose:
        print('add kh and kv from model layer dataset to modflow model')

    kh, kv = get_kh_kv(kh_arr, kv_arr, anisotropy,
                       fill_value_kh=fill_value_kh,
                       fill_value_kv=fill_value_kv,
                       verbose=verbose)

    model_ds['kh'] = xr.ones_like(model_ds['idomain']) * kh

    model_ds['kv'] = xr.ones_like(model_ds['idomain']) * kv

    return model_ds


def get_kh_kv(kh_in, kv_in, anisotropy,
              fill_value_kh=1.0, fill_value_kv=1.0,
              verbose=False):
    """maak kh en kv rasters voor flopy vanuit een regis raster met nan
    waardes.

    vul kh raster door:
    1. pak kh uit regis, tenzij nan dan:
    2. pak kv uit regis vermenigvuldig met anisotropy, tenzij nan dan:
    3. pak fill_value_kh

    vul kv raster door:
    1. pak kv uit regis, tenzij nan dan:
    2. pak kh uit regis deel door anisotropy, tenzij nan dan:
    3. pak fill_value_kv

    Supports structured and unstructured grids.

    Parameters
    ----------
    kh_in : np.ndarray
        kh from regis with nan values shape(nlay, nrow, ncol) or 
        shape(nlay, len(cid))
    kv_in : np.ndarray
        kv from regis with nan values shape(nlay, nrow, ncol) or 
        shape(nlay, len(cid))
    anisotropy : int or float
        factor to calculate kv from kh or the other way around
    fill_value_kh : int or float, optional
        use this value for kh if there is no data in regis. The default is 1.0.
    fill_value_kv : int or float, optional
        use this value for kv if there is no data in regis. The default is 1.0.
    verbose : bool, optional
        print additional information. default is False

    Returns
    -------
    kh_out : np.ndarray
        kh without nan values (nlay, nrow, ncol) or shape(nlay, len(cid))
    kv_out : np.ndarray
        kv without nan values (nlay, nrow, ncol) or shape(nlay, len(cid))
    """
    kh_out = np.zeros_like(kh_in)
    for i, kh_lay in enumerate(kh_in):
        kh_new = kh_lay.copy()
        kv_new = kv_in[i].copy()
        if ~np.all(np.isnan(kh_new)):
            if verbose:
                print(f'layer {i} has a kh')
            kh_out[i] = np.where(np.isnan(kh_new), kv_new * anisotropy, kh_new)
            kh_out[i] = np.where(np.isnan(kh_out[i]), fill_value_kh, kh_out[i])
        elif ~np.all(np.isnan(kv_new)):
            if verbose:
                print(f'layer {i} has a kv')
            kh_out[i] = np.where(
                np.isnan(kv_new), fill_value_kh, kv_new * anisotropy)
        else:
            if verbose:
                print(f'kv and kh both undefined in layer {i}')
            kh_out[i] = fill_value_kh

    kv_out = np.zeros_like(kv_in)
    for i, kv_lay in enumerate(kv_in):
        kv_new = kv_lay.copy()
        kh_new = kh_in[i].copy()
        if ~np.all(np.isnan(kv_new)):
            if verbose:
                print(f'layer {i} has a kv')
            kv_out[i] = np.where(np.isnan(kv_new), kh_new / anisotropy, kv_new)
            kv_out[i] = np.where(np.isnan(kv_out[i]), fill_value_kv, kv_out[i])
        elif ~np.all(np.isnan(kh_new)):
            if verbose:
                print(f'layer {i} has a kh')
            kv_out[i] = np.where(
                np.isnan(kh_new), fill_value_kv, kh_new / anisotropy)
        else:
            if verbose:
                print(f'kv and kh both undefined in layer {i}')
            kv_out[i] = fill_value_kv

    return kh_out, kv_out


def add_top_bot_to_model_ds(ml_layer_ds, model_ds,
                            nodata=None,
                            max_per_nan_bot=50,
                            gridtype='structured',
                            verbose=False):
    """ add top and bot from a model layer dataset to THE model dataset.

    Supports structured and unstructured grids.

    Parameters
    ----------
    ml_layer_ds : xarray.Dataset
        dataset with model layer data with a top and bottom
    model_ds : xarray.Dataset
        dataset with model data where top and bot are added to
    nodata : int, optional
        if the first_active_layer data array in model_ds has this value,
        it means this cell is inactive in all layers. If nodata is None the
        nodata value in model_ds is used.
        the default is None
    max_per_nan_bot : int or float, optional
        if the percentage of cells that have nan values in all layers is
        higher than this an error is raised. The default is 50.
    gridtype : str, optional
        type of grid, options are 'structured' and 'unstructured'. 
        The default is 'structured'.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data including top and bottom

    """
    if nodata is None:
        nodata = model_ds.attrs['nodata']

    if verbose:
        print('using top and bottom from model layers dataset for modflow model')
        print('replace nan values for inactive layers with dummy value')

    if gridtype == 'structured':

        model_ds = add_top_bot_structured(ml_layer_ds, model_ds,
                                          nodata=nodata,
                                          max_per_nan_bot=max_per_nan_bot)

    elif gridtype == 'unstructured':
        model_ds = add_top_bot_unstructured(ml_layer_ds, model_ds,
                                            nodata=nodata,
                                            max_per_nan_bot=max_per_nan_bot)

    return model_ds


def add_top_bot_unstructured(ml_layer_ds, model_ds, nodata=-999,
                             max_per_nan_bot=50):
    """ voeg top en bottom vanuit model layer dataset toe aan de model dataset

    Deze functie is bedoeld voor unstructured arrays in modflow 6

    stappen:
    1. Zorg dat de onderste laag altijd een bodemhoogte heeft, als de bodem
    van alle bovenliggende lagen nan is, pak dan 0.
    2. Zorg dat de top van de bovenste laag altijd een waarde heeft, als de
    top van alle onderligende lagen nan is, pak dan 0.
    3. Vul de nan waarden in alle andere lagen door:
        a. pak bodem uit regis, tenzij nan dan:
        b. gebruik bodem van de laag erboven (of de top voor de bovenste laag)

    Supports only unstructured grids.

    Parameters
    ----------
    ml_layer_ds : xarray.Dataset
        dataset with model layer data with a top and bottom
    model_ds : xarray.Dataset
        dataset with model data where top and bottom are added to
    nodata : int, optional
        if the first_active_layer data array in model_ds has this value,
        it means this cell is inactive in all layers 
    max_per_nan_bot : int or float, optional
        if the percentage of cells that have nan values in all layers. 
        The default is 50.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data including top and bottom
    """
    # step 1:
    # set nan-value in bottom array
    # set to zero if value is nan in all layers
    # set to minimum value of all layers if there is any value in any layer
    active_domain = model_ds['first_active_layer'].data != nodata

    lowest_bottom = ml_layer_ds['bottom'].data[-1].copy()
    if np.any(active_domain == False):
        percentage = 100 * (active_domain == False).sum() / \
            (active_domain.shape[0])
        if percentage > max_per_nan_bot:
            print(f"{percentage:0.1f}% cells with nan values in every layer"
                  " check your extent")
            raise MemoryError('if you want to'
                              ' ignore this error set max_per_nan_bot higher')

        # set bottom to zero if bottom in a cell is nan in all layers
        lowest_bottom = np.where(active_domain, lowest_bottom, 0)

    if np.any(np.isnan(lowest_bottom)):
        # set bottom in a cell to lowest bottom of all layers
        i_nan = np.where(np.isnan(lowest_bottom))
        for i in i_nan:
            val = np.nanmin(ml_layer_ds['bottom'].data[:, i])
            lowest_bottom[i] = val
            if np.isnan(val):
                raise ValueError(
                    'this should never happen please contact Artesia')

    # step 2: get highest top values of all layers without nan values
    highest_top = ml_layer_ds['top'].data[0].copy()
    if np.any(np.isnan(highest_top)):
        highest_top = np.where(active_domain, highest_top, 0)

    if np.any(np.isnan(highest_top)):
        i_nan = np.where(np.isnan(highest_top))
        for i in i_nan:
            val = np.nanmax(ml_layer_ds['top'].data[:, i])
            highest_top[i] = val
            if np.isnan(val):
                raise ValueError(
                    'this should never happen please contact Artesia')

    # step 3: fill nans in all layers
    nlay = model_ds.dims['layer']
    top_bot_raw = np.ones((nlay + 1, model_ds.dims['cid']))
    top_bot_raw[0] = highest_top
    top_bot_raw[1:-1] = ml_layer_ds['bottom'].data[:-1].copy()
    top_bot_raw[-1] = lowest_bottom
    top_bot = np.ones_like(top_bot_raw)
    for i_from_bot, blay in enumerate(top_bot_raw[::-1]):
        i_from_top = nlay - i_from_bot
        new_lay = blay.copy()
        if np.any(np.isnan(new_lay)):
            lay_from_bot = i_from_bot
            lay_from_top = nlay - lay_from_bot
            while np.any(np.isnan(new_lay)):
                new_lay = np.where(np.isnan(new_lay),
                                   top_bot_raw[lay_from_top],
                                   new_lay)
                lay_from_bot += 1
                lay_from_top = nlay - lay_from_bot

        top_bot[i_from_top] = new_lay

    model_ds['bot'] = xr.DataArray(top_bot[1:], dims=('layer', 'cid'),
                                   coords={'cid': model_ds.cid.data,
                                           'layer': model_ds.layer.data})
    model_ds['top'] = xr.DataArray(top_bot[0], dims=('cid'),
                                   coords={'cid': model_ds.cid.data})

    return model_ds


def add_top_bot_structured(ml_layer_ds, model_ds, nodata=-999,
                           max_per_nan_bot=50):
    """ voeg top en bottom vanuit een model layer dataset toe aan DE model 
    dataset

    Deze functie is bedoeld voor structured arrays in modflow 6

    stappen:
    1. Zorg dat de onderste laag altijd een bodemhoogte heeft, als de bodem
    van alle bovenliggende lagen nan is, pak dan 0.
    2. Zorg dat de top van de bovenste laag altijd een waarde heeft, als de
    top van alle onderligende lagen nan is, pak dan 0.
    3. Vul de nan waarden in alle andere lagen door:
        a. pak bodem uit de model layer dataset, tenzij nan dan:
        b. gebruik bodem van de laag erboven (of de top voor de bovenste laag)

    Supports only structured grids.

    Parameters
    ----------
    ml_layer_ds : xarray.Dataset
        dataset with model layer data with a top and bottom
    model_ds : xarray.Dataset
        dataset with model data where top and bottom are added to
    nodata : int, optional
        if the first_active_layer data array in model_ds has this value,
        it means this cell is inactive in all layers 
    max_per_nan_bot : int or float, optional
        if the percentage of cells that have nan values in all layers. 
        The default is 50.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data including top and bottom

    """

    active_domain = model_ds['first_active_layer'].data != nodata

    # step 1:
    # set nan-value in bottom array
    # set to zero if value is nan in all layers
    # set to minimum value of all layers if there is any value in any layer
    lowest_bottom = ml_layer_ds['bottom'].data[-1].copy()
    if np.any(active_domain == False):
        percentage = 100 * (active_domain == False).sum() / \
            (active_domain.shape[0] * active_domain.shape[1])
        if percentage > max_per_nan_bot:
            print(f"{percentage:0.1f}% cells with nan values in every layer"
                  " check your extent")
            raise MemoryError('if you want to'
                              ' ignore this error set max_per_nan_bot higher')

        # set bottom to zero if bottom in a cell is nan in all layers
        lowest_bottom = np.where(active_domain, lowest_bottom, 0)

    if np.any(np.isnan(lowest_bottom)):
        # set bottom in a cell to lowest bottom of all layers
        rc_nan = np.where(np.isnan(lowest_bottom))
        for row, col in zip(rc_nan[0], rc_nan[1]):
            val = np.nanmin(ml_layer_ds['bottom'].data[:, row, col])
            lowest_bottom[row, col] = val
            if np.isnan(val):
                raise ValueError(
                    'this should never happen please contact Onno')

    # step 2: get highest top values of all layers without nan values
    highest_top = ml_layer_ds['top'].data[0].copy()
    if np.any(np.isnan(highest_top)):
        # set top to zero if top in a cell is nan in all layers
        highest_top = np.where(active_domain, highest_top, 0)

    if np.any(np.isnan(highest_top)):
        # set top in a cell to highest top of all layers
        rc_nan = np.where(np.isnan(highest_top))
        for row, col in zip(rc_nan[0], rc_nan[1]):
            val = np.nanmax(ml_layer_ds['top'].data[:, row, col])
            highest_top[row, col] = val
            if np.isnan(val):
                raise ValueError(
                    'this should never happen please contact Onno')

    # step 3: fill nans in all layers
    nlay = model_ds.dims['layer']
    nrow = model_ds.dims['y']
    ncol = model_ds.dims['x']
    top_bot_raw = np.ones((nlay + 1, nrow, ncol))
    top_bot_raw[0] = highest_top
    top_bot_raw[1:-1] = ml_layer_ds['bottom'].data[:-1].copy()
    top_bot_raw[-1] = lowest_bottom
    top_bot = np.ones_like(top_bot_raw)
    for i_from_bot, blay in enumerate(top_bot_raw[::-1]):
        i_from_top = nlay - i_from_bot
        new_lay = blay.copy()
        if np.any(np.isnan(new_lay)):
            lay_from_bot = i_from_bot
            lay_from_top = nlay - lay_from_bot
            while np.any(np.isnan(new_lay)):
                new_lay = np.where(np.isnan(new_lay),
                                   top_bot_raw[lay_from_top],
                                   new_lay)
                lay_from_bot += 1
                lay_from_top = nlay - lay_from_bot

        top_bot[i_from_top] = new_lay

    model_ds['bot'] = xr.DataArray(top_bot[1:], dims=('layer', 'y', 'x'),
                                   coords={'x': model_ds.x.data,
                                           'y': model_ds.y.data,
                                           'layer': model_ds.layer.data})

    model_ds['top'] = xr.DataArray(top_bot[0], dims=('y', 'x'),
                                   coords={'x': model_ds.x.data,
                                           'y': model_ds.y.data})

    return model_ds


def fill_top_bot_kh_kv_at_mask(model_ds, fill_mask,
                               gridtype='structured',
                               gridprops=None):
    """ fill values in top, bot, kh and kv where:
        1. the cell is True in mask
        2. the cell thickness is greater than 0

    fill values:
        top: 0
        bot: minimum of bottom_filled or top
        kh: kh_filled if thickness is greater than 0
        kv: kv_filled if thickness is greater than 0

    Parameters
    ----------
    model_ds : xr.DataSet
        model dataset, should contain 'first_active_layer'
    fill_mask : xr.DataArray
        1 where a cell should be replaced by masked value.
    gridtype : str, optional
        type of grid.        
    gridprops : dictionary, optional
        dictionary with grid properties output from gridgen. Default is None

    Returns
    -------
    model_ds : xr.DataSet
        model dataset with adjusted data variables:
            'top', 'bot', 'kh', 'kv'
    """

    # zee cellen hebben altijd een top gelijk aan 0
    model_ds['top'] = xr.where(fill_mask, 0, model_ds['top'])

    if gridtype == 'structured':
        fill_function = mgrid.fillnan_dataarray_structured_grid
        fill_function_kwargs = {}
    elif gridtype == 'unstructured':
        fill_function = mgrid.fillnan_dataarray_unstructured_grid
        fill_function_kwargs = {'gridprops': gridprops}

    for lay in range(model_ds.dims['layer']):
        bottom_nan = xr.where(fill_mask, np.nan, model_ds['bot'][lay])
        bottom_filled = fill_function(bottom_nan, **fill_function_kwargs)[0]

        kh_nan = xr.where(fill_mask, np.nan, model_ds['kh'][lay])
        kh_filled = fill_function(kh_nan, **fill_function_kwargs)[0]

        kv_nan = xr.where(fill_mask, np.nan, model_ds['kv'][lay])
        kv_filled = fill_function(kv_nan, **fill_function_kwargs)[0]

        if lay == 0:
            # top ligt onder bottom_filled -> laagdikte wordt 0
            # top ligt boven bottom_filled -> laagdikte o.b.v. bottom_filled
            mask_top = model_ds['top'] < bottom_filled
            model_ds['bot'][lay] = xr.where(fill_mask * mask_top,
                                            model_ds['top'],
                                            bottom_filled)
            model_ds['kh'][lay] = xr.where(fill_mask * mask_top,
                                           model_ds['kh'][lay],
                                           kh_filled)
            model_ds['kv'][lay] = xr.where(fill_mask * mask_top,
                                           model_ds['kv'][lay],
                                           kv_filled)

        else:
            # top ligt onder bottom_filled -> laagdikte wordt 0
            # top ligt boven bottom_filled -> laagdikte o.b.v. bottom_filled
            mask_top = model_ds['bot'][lay - 1] < bottom_filled
            model_ds['bot'][lay] = xr.where(fill_mask * mask_top,
                                            model_ds['bot'][lay - 1],
                                            bottom_filled)
            model_ds['kh'][lay] = xr.where(fill_mask * mask_top,
                                           model_ds['kh'][lay],
                                           kh_filled)
            model_ds['kv'][lay] = xr.where(fill_mask * mask_top,
                                           model_ds['kv'][lay],
                                           kv_filled)
       

    return model_ds


def add_bathymetry_to_top_bot_kh_kv(model_ds, bathymetry,
                                    fill_mask,
                                    kh_sea=10,
                                    kv_sea=10):
    """ add bathymetry to the top and bot of each layer for all cells with
    fill_mask.

    Parameters
    ----------
    model_ds : xarray.Dataset
        dataset with model data, should 
    fill_mask : xr.DataArray
        cell value is 1 if you want to add bathymetry

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data where the top, bot, kh and kv are changed

    """
    model_ds['top'] = xr.where(fill_mask,
                               0.0,
                               model_ds['top'])

    lay = 0
    model_ds['bot'][lay] = xr.where(fill_mask,
                                    bathymetry,
                                    model_ds['bot'][lay])

    model_ds['kh'][lay] = xr.where(fill_mask,
                                   kh_sea,
                                   model_ds['kh'][lay])

    model_ds['kv'][lay] = xr.where(fill_mask,
                                   kv_sea,
                                   model_ds['kv'][lay])

    # reset bot for all layers based on bathymetrie
    for lay in range(1, model_ds.dims['layer']):
        model_ds['bot'][lay] = np.where(model_ds['bot'][lay] > model_ds['bot'][lay - 1],
                                        model_ds['bot'][lay - 1],
                                        model_ds['bot'][lay])

    return model_ds

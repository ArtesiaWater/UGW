# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:14:50 2020

@author: oebbe
"""

import os

import xarray as xr
import numpy as np

import rioxarray


def data_array_3d_to_gis(da, gisdir, model_name, da_name=None, 
                         nan_value=None):
    
    if da_name is None:
        da_name = da.name
        
    for i, lay in enumerate(da.layer[:]):
        da_exp = da.sel(layer=lay)
        data_array_2d_to_gis(da_exp, gisdir, model_name, da_name=da_name,
                             lay_name=str(lay.values), lay_no=i+1,
                             nan_value=nan_value)
    
    return 1

def data_array_2d_to_gis(da, gisdir, model_name, da_name=None,
                         lay_name='', lay_no=1, nan_value=None):
    
    if da_name is None:
        da_name = da.name
        
    if nan_value is not None:
        da = xr.where(da==nan_value, np.nan, da)
        
    if da.isnull().all():
        print(f'{da_name} layer {lay_name} contains only nan values -> no export was made')
    else:
        da = da.rio.set_crs(28992)
        fname = os.path.join(gisdir,model_name,f'{da_name}_{lay_no:03d}_{lay_name}.tif')
        da.rio.to_raster(fname)
    
    return 1 
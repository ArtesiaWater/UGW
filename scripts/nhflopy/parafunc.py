# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 13:45:24 2020

@author: oebbe
"""
import os

import xarray as xr
import numpy as np

from scipy import interpolate
from pyevtk.hl import gridToVTK





def model_to_vtk(model_ds, paradir, fname, 
                 layer=14,
                 datavar='kh', nan_factor=0.01):
    """
    
    Parameters
    ----------
    model_ds : xarray DataSet
        model dataset within the extent that is exported to vtk.
    paradir : str
        directory to save the vtk file.
    fname : str
        filename of the vtk file.
    layer : int
        number of the regis layer you want to convert to vtk
    datavar : str
        name of the data variabel you want to convert to vtk
    nan_factor : float
        if the value in a cell is more than (0.01) part influenced by 
        nan values make it a nan value

    Returns
    -------
    1 if succesfull.

    """
    X = np.append(model_ds.x.values-(model_ds.delr/2), model_ds.x.values[-1]+(model_ds.delr/2))
    Y = np.append(model_ds.y.values[::-1]-(model_ds.delc/2), model_ds.y.values[::-1][-1]+(model_ds.delc/2))
    Z = np.linspace(-1,1,2)
    
    x,y,z = np.meshgrid(X,Y,Z, indexing='ij')
    x = x.astype(float)
    y = y.astype(float)
    if layer==0:
        top = project_to_grid(model_ds['top'], X, Y, nan_factor)
    else:
        top = project_to_grid(model_ds['bot'][layer-1], X, Y, nan_factor)
        
    bot = project_to_grid(model_ds['bot'][layer], X, Y, nan_factor)
    z[:,:,0] = (top.T*100)
    z[:,:,1] = (bot.T*100)
    
    arr = model_ds[datavar][layer:layer+1].values
    arr[0] = arr[0][::-1]
    arr = arr.T
    
    
    gridToVTK(os.path.join(paradir, fname), x, y, z, cellData={datavar:arr})
    
    return 1


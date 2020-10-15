# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:17:28 2020

@author: danie
"""

import rasterio
import numpy as np
from rasterio.warp import Resampling

raster_in = r'd:\repositories\UGW\data\sources\Rijkswaterstaat\ijg2020_v1101\w001001.adf'
raster_uit = r'd:\repositories\UGW\data\sources\Bathymetry\Rijkswaterstaat\ijg2020_v1101.tif'
shape_uit = r'd:\repositories\UGW\data\sources\Bathymetry\Rijkswaterstaat\ijg2020_v1101.shp'


with rasterio.open(raster_in) as src:
    profile = src.profile
    src_nodata = profile['nodata']
    profile['driver'] = 'GTiff'
    profile['compress'] ='lzw'
    profile['predictor'] = 2
    profile['count'] = 1
    
    with rasterio.open(raster_uit,'w', **profile) as dst:
        print('write data (blocks in memory)')
        for block_index, window in src.block_windows(1):
            data = src.read(window=window)[0]
            dst.write(data, window=window,indexes=1)
            dst.scales  = [0.01,]
            
        #%%
        print('create overviews')
        cell_size = profile['transform'][0]
        factors = [int(size/cell_size) for size in [5] if size > cell_size]
        dst.build_overviews(factors, Resampling.average)
        dst.update_tags(ns='rio_overview', resampling='average')
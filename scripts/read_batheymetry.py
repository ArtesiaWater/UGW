# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:30:15 2020

@author: danie
"""

import rasterio
from rasterio.enums import Resampling
import numpy as np

bag_file = r'..\data\sources\Rijkswaterstaat\MN_zuid_NAP.bag'
tif_file = r'..\data\MN_zuid_NAP.tif'
nodata = 32767
dtype = rasterio.int16
scale = 0.001

print('read bag-file')
with rasterio.open(bag_file) as src:
    profile = src.profile
    data = src.read(1)

data = np.where(data == profile['nodata'], nodata, data / scale).astype(dtype)

#%%
profile['driver'] = 'GTiff'
profile['predictor'] = 2
profile['dtype'] = dtype
profile['nodata'] = nodata
profile['count'] = 1

print('write tif-file')
with rasterio.open(tif_file,'w', **profile) as dst:
    cell_size = profile['transform'][0]
    dst.scales  = [scale,]
    dst.write(data,1)
    factors = [int(size/cell_size) for size in [5] if size > cell_size]
    dst.build_overviews(factors, Resampling.average)
    dst.update_tags(ns='rio_overview', resampling='average')
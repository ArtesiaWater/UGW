# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:30:15 2020

@author: danie
"""
import glob
import os
import rasterio
from rasterio.enums import Resampling
import numpy as np

bag_dir = r'..\data\sources\Rijkswaterstaat\BAG'
tif_dir = r'..\data\sources\Bathymetry\Rijkswaterstaat'
nodata = 32767
dtype = rasterio.int16
scale = 0.001

block_writing = True #use if larger than memory rasters
over_write = True

#%%
for bag_file in glob.glob(os.path.join(bag_dir, '*.bag')):
    file_name = os.path.splitext(os.path.basename(bag_file))[0]
    tif_file = os.path.join(tif_dir,'{}.tif'.format(file_name))
    if over_write or not os.path.exists(tif_file):
        print('converting {}'.format(bag_file))
        with rasterio.open(bag_file) as src:
            profile = src.profile
            src_nodata = profile['nodata']
            profile['driver'] = 'GTiff'
            profile['compress'] ='lzw'
            profile['predictor'] = 2
            profile['dtype'] = dtype
            profile['nodata'] = nodata
            profile['count'] = 1
            with rasterio.open(tif_file,'w', **profile) as dst:
                if block_writing:
                    print('write data (blocks in memory)')
                    for block_index, window in src.block_windows(1):
                        data = src.read(window=window)[0]
                        data = np.where(data == src_nodata, nodata, data / scale).astype(dtype)
                        dst.write(data, window=window,indexes=1)
                        
                else:
                    print('write data (all in memory)')
                    data = src.read(1)
                    data = np.where(data == src_nodata, nodata, data / scale).astype(dtype)
                    dst.write(data,1)
                
                dst.scales  = [scale,]
                
                #%%
                print('create overviews')
                cell_size = profile['transform'][0]
                factors = [int(size/cell_size) for size in [5] if size > cell_size]
                dst.build_overviews(factors, Resampling.average)
                dst.update_tags(ns='rio_overview', resampling='average')
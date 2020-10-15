# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:34:48 2020

@author: danie
"""

import rasterio.features
import rasterio
from rasterio.enums import Resampling
import glob
import numpy as np
import os
import geopandas as gpd
from shapely.ops import unary_union
from shapely.geometry import Polygon, MultiPolygon

input_dir = r'..\data\sources\Rijnland\bathymetrie'
output_dir = r'..\data\sources\Bathymetry\Rijnland'
shape_file = os.path.join(os.path.abspath(input_dir),'Rijnand_meren_plassen_diepte.shp')
tif_files = glob.glob(os.path.join(os.path.abspath(input_dir),'*.tif'))

tif_specs = {'Diepte_Langeraarse_Plassen.tif':{'peil':-1.57,
                                               'name':'Langeraarse Plassen'},
            'Diepte_Westeinder_Plassen_1m.tif':{'peil':-0.61,
                                                'name':None}}   

scale = 0.01
dtype = rasterio.int16
nodata = 32767

#%%
def match_file(x):
    result = [key for key, value in tif_specs.items() if value['name'] == x]
    if len(result) != 0:
        if result[0] != None:
            return os.path.join(output_dir,result[0])
        else:
            return np.NaN     

#%% tif-files
for tif_file in tif_files:
    print('converting {}'.format(tif_file))
    file_name = os.path.basename(tif_file)
    level = tif_specs[file_name]['peil']
    
    with rasterio.open(tif_file) as src:
        data = src.read(1)
        profile = src.profile
              
        data = np.where(data == profile['nodata'], nodata, (level - data)/scale).astype(dtype)
        
        profile['compress'] ='lzw'
        profile['predictor'] = 2
        profile['dtype'] = dtype
        profile['nodata'] = nodata
        
        with rasterio.open(os.path.join(output_dir,file_name), 'w', **profile) as dst:
            dst.write(data,1)
            dst.scales  = [scale,]
            
            print('create overviews')
            cell_size = profile['transform'][0]
            factors = [int(size/cell_size) for size in [5] if size > cell_size]
            dst.build_overviews(factors, Resampling.average)
            dst.update_tags(ns='rio_overview', resampling='average')

#%% shape-file
gdf = gpd.read_file(shape_file)

gdf['BL_AVG'] = gdf.apply((lambda x: x['ZP'] - x['diepte_gem'] if x['diepte_gem'] != 0 else np.NaN),axis=1)
gdf['BL_MIN'] = gdf.apply((lambda x: x['ZP'] - x['diepte_max'] if x['diepte_max'] != 0 else np.NaN),axis=1)
gdf['FILE'] = gdf.apply((lambda x: match_file(x['NAAM'])),axis=1)
drop_cols = [col for col in gdf.columns if not col in ['BL_AVG','BL_MIN','FILE','geometry']]
gdf = gdf.drop(drop_cols,axis=1)
gdf.to_file(os.path.join(output_dir,'bathymetry.shp'))
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:34:48 2020

@author: danie
"""

import json
import rasterio.features
import rasterio
from rasterio.warp import Resampling
import glob
import numpy as np
import os
import geopandas as gpd
from shapely.ops import unary_union
from shapely.geometry import Polygon, MultiPolygon

mask_dict = {'FILE':[],
             'BL_AVG':[],
             'BL_MIN':[],
             'geometry':[]}
cell_size = 5

tif_dirs = [r'..\data\sources\Bathymetry\Rijkswaterstaat',
            r'..\data\sources\Bathymetry\Waternet',
            r'..\data\sources\Bathymetry\Rijnland']

tif_dirs = [r'..\data\sources\Bathymetry\ijsselmeer']

mask_shape = r'..\data\sources\Bathymetry\bathymetry_masks.shp'
mask_shape = r'..\data\sources\Bathymetry\bathymetry_masks2.shp'

tif_files = [x for sublist in [glob.glob(os.path.join(tif_dir, '*.tif')) 
                               for tif_dir in tif_dirs] for x in sublist]

mask_files = [x for sublist in [glob.glob(os.path.join(tif_dir, '*.shp')) 
                               for tif_dir in tif_dirs] for x in sublist]

#%% add properties of bathymetry masks
for mask in mask_files:
    m_dict = json.loads(gpd.read_file(mask).to_json())
    for feature in m_dict['features']:
        for item in list(mask_dict.keys()):
            if item == 'geometry':
                mask_dict[item].append(gpd.GeoDataFrame.from_features([feature]).geometry.values[0])
            else:   
                mask_dict[item].append(feature['properties'][item])

skip_files = [file for file in mask_dict['FILE'] if not file == None]
      
#%%
for tif_file in tif_files:
    if not tif_file in skip_files:
        print('adding {}'.format(tif_file))
        mask_dict['FILE'].append(tif_file)
        
        with rasterio.open(tif_file) as src:
            mask = None
            upscale_factor = cell_size/src.profile['transform'][0]
            #data = src.read(1) # first band
            data = src.read(
                        out_shape=(
                            src.profile['count'],
                            int(src.profile['height'] / upscale_factor),
                            int(src.profile['width'] / upscale_factor)
                            ),
                        resampling=Resampling.average
                    )[0]
            
            xmin, ymin, xmax, ymax = src.bounds
            transform = rasterio.transform.from_bounds(xmin,ymin,xmax,ymax, data.shape[1], data.shape[0])
            
            data = np.where(data == src.profile['nodata'], 0, 1)
            results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                rasterio.features.shapes(data, mask=data.astype(np.uint8), transform=transform)))

            mask_dict['BL_MIN'].append(np.NaN)
            mask_dict['BL_AVG'].append(np.NaN)
            
            data = None
            
            geoms = list(results)
            
            gdf = gpd.GeoDataFrame.from_features(geoms)
            gdf = gdf.simplify(20, preserve_topology=True)
            
            geometry = unary_union(gdf.geometry.values)
            mask_dict['geometry'].append(geometry)
        
#%% write shape-file
mask_gdf = gpd.GeoDataFrame(data=mask_dict)
mask_gdf.crs = 'epsg:28992'

schema = {'geometry': 'Polygon',
          'properties': {'FILE':'str',
                         'BL_AVG':'float: 5.2',
                         'BL_MIN':'float: 5.2'}
          }

mask_gdf.to_file(mask_shape, schema=schema)

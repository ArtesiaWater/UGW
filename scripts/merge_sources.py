# -*- coding: utf-8 -*-

import pandas as pd
import geopandas as gpd
import pyproj
from owslib.wfs import WebFeatureService
import requests
import json


admins = json.loads(open(r'..\config\administrations.json','r').read())
sources = json.loads(open(r'..\config\sources.json','r').read())

#%% 
def set_crs(gdf):
    default_crs = sources['default_crs']
    if gdf.crs == None: 
        gdf.crs = default_crs
    else:
       if not pyproj.CRS(gdf.crs).equals(pyproj.CRS(sources['default_crs'])):
           gdf = gdf.to_crs(sources['default_crs'])
           
    return gdf

class WFS:
    
    def __init__(self,url,crs='epsg:28992'):
        self.wfs = WebFeatureService(url,version='2.0.0', timeout=300)
        self.layers = list(self.wfs.contents.keys())
        self.crs = 'epsg:28992'
        
    def get_features(self,layer,poly):
        if layer in self.layers:
            bbox = poly.bounds
            response = self.wfs.getfeature(typename=layer, 
                          bbox=bbox, 
                          outputFormat='json')
            geojson = json.loads(response.read())
            gdf = set_crs(gpd.GeoDataFrame.from_features(geojson['features']))
            gdf = gdf[gdf.within(poly)]
            
            return gdf
        else:
            print('{} not in layers({})'.format(layer,self.layers))
            
            return None

class HyDAMO:
    
    def __init__(self):
        self.wfs = WebFeatureService('https://data.nhi.nu/geoserver/ows',version='2.0.0', timeout=300)
        self.layers = list(self.wfs.contents.keys())
        self.crs = 'epsg:28992' 
   
    def get_features(self,layer,poly):
        if layer in self.layers:
            bbox = poly.bounds
            response = self.wfs.getfeature(typename=layer, 
                          bbox=bbox, 
                          outputFormat='json')
            geojson = json.loads(response.read())
            gdf = set_crs(gpd.GeoDataFrame.from_features(geojson['features']))
            gdf = gdf[gdf.within(poly)]
            
            return gdf
        else:
            print('{} not in layers({})'.format(layer,self.layers))
            
            return None
            

class ArcREST:
    
    def __init__(self,url,crs='epsg:28992'):
        epsg = crs[crs.find(':')+1:]
        self.maxRecordCount = requests.get('{}?f=pjson'.format(url)).json()['maxRecordCount']
        self.url = url
        self.epsg = epsg
        self.query = '{url}/query?where={min_objects}<=OBJECTID and {max_objects}>OBJECTID&outFields=*&outSR={epsg}&f=geojson'
        self.crs = crs
        
    def get_features(self):
        
        url = self.query.format(url=self.url,min_objects=0,max_objects=1000000,epsg=self.epsg)
        url = '{}&returnIdsOnly=true'.format(url)
        print(url)
        response = requests.get(url)
        object_ids = response.json()['properties']['objectIds']
        object_ids.sort()
        
        downloads = round(len(object_ids)/self.maxRecordCount + 0.5)
        gdf_list = []
        for download in range(downloads):
            min_objects = download * self.maxRecordCount
            max_objects = min_objects + self.maxRecordCount
            url = self.query.format(url=self.url,min_objects=min_objects,max_objects=max_objects,epsg=self.epsg)
            response = requests.post(url)
            gdf = gpd.GeoDataFrame.from_features(response.json()['features'])
            gdf.crs = self.crs
            gdf_list += [gdf]
            
        if len(gdf_list) > 1:
            gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True))
            
        else: gdf = gdf_list[0]
        
        return(gdf)
        
        
#%% read admin_boundaries and update settings (later a class)
admin_src = sources['boundaries']['administrations']
id_field = admin_src['id_field']
gdf = gpd.read_file(admin_src['path'])
gdf = set_crs(gdf)
{key: value.update({'boundary':gdf.loc[gdf[id_field] == value['id']].geometry.item()}) for key, value in admins.items()}

#%%
#init water_bodies
gdf_list = []

#add admin_id and geometry to water_bodies
for key, value in admins.items():
    print('waterlopen {}'.format(key))
    admin_id = value['id']
    src = sources['water_lines'][value['water_src']]
    src_type = src['type']
    poly = value['boundary']
          
    if src_type == 'wfs':
        url = src['url']
        layer = src['layer']
        wfs = WFS(url)
        gdf = wfs.get_features(layer,poly)
        
    elif src_type == 'arcrest':
        url = src['url']
        rest = ArcREST(url)
        gdf = rest.get_features()
           
    elif src_type =='file':
        gdf = set_crs(gpd.read_file(src['path']))
        
    else:
        print("type '{}' not supported".format(src_type))
    
    drop_cols = [col for col in list(gdf.columns) if col != 'geometry']
    gdf = gdf.drop(drop_cols, axis=1)
    gdf['admin_id'] = admin_id
    gdf['admin'] = key
    gdf_list += [gdf]

#merge all into one dataframe    
gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True))
gdf.crs = sources['default_crs']

gdf.to_file('..\data\waterlopen.shp')   
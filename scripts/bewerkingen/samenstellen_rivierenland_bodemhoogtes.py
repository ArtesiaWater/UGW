# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 11:38:59 2020

@author: danie
"""


import geopandas as gpd
import pandas as pd
import sys
sys.path.insert(1, r'../')
import services
import logging
import os

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

#%% paden
reaches_shp = r'..\..\data\sources\Rivierenland\reaches.shp'
admins_shp = r'..\..\data\sources\admins.shp'
watergangen_shp = r'..\..\data\sources\Rivierenland\waterlijnen_20200724.shp'
url = r'https://kaarten.wsrl.nl/arcgis/rest/services/Kaarten/Legger_Watersysteem/MapServer'
layer = 13

#%% inlezen reaches & admins
reaches_gdf = gpd.read_file(reaches_shp)
admins_gdf = gpd.read_file(admins_shp)

#%% downloaden legger
rest = services.ArcREST(url)
poly = admins_gdf[admins_gdf['name'] == 'Waterschap Rivierenland']['geometry'].values[0]
watergangen_gdf = rest.get_features(layer,poly,object_id='se_sdo_rowid')

#%% combineren reaches met legger

reaches_filtered_gdf = reaches_gdf[(reaches_gdf['BEDLEVELUP'] != 0) & (reaches_gdf['BEDLEVELDN'] != 0)]

reaches_filtered_gdf = reaches_filtered_gdf.drop(['ROUTENUM', 
                                                  'BEDWIDTHUP',
                                                  'BEDWIDTHDN', 
                                                  'MAXWIDTHUP', 
                                                  'MAXWIDTHDN', 
                                                  'TLUP', 
                                                  'TLDN', 
                                                  'TLWIDTHUP',
                                                  'TLWIDTHDN'], axis=1)


codes = np.unique(np.array([value[0:6] for value in reaches_filtered_gdf['ID'].values]))

watergangen_filtered_gdf = watergangen_gdf[~watergangen_gdf['se_sdo_rowid'].isin(codes)]
watergangen_filtered_gdf = watergangen_filtered_gdf.rename(columns={'se_sdo_rowid':'ID'})

#%%
watergangen_totaal = gpd.GeoDataFrame(pd.concat([reaches_filtered_gdf,watergangen_filtered_gdf],ignore_index=True))
#%%
watergangen_totaal.to_file(watergangen_shp)
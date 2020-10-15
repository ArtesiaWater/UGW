# -*- coding: utf-8 -*-
"""
Created on Tue May 26 10:07:37 2020

@author: danie
"""

import pandas as pd
import geopandas as gpd
import get_ahn
import numpy as np
import datetime
from rasterstats import zonal_stats
import rasterio

'''
ToDo:
    -   .str.contains(profiel_ident) werkend maken voor situaties waarin juist de omschrijving
        contained moet zijn in de profiel_ident. bijvoorbeeld L9(nv)L(nv)R (profiel ident) 
        is contained in de omschrijving L9. En P17L5R5 (profiel_ident) is contained 
        in de omschrijving P17
'''


legger_xlsx = r'..\..\data\sources\zuiderzeeland\leggertabellen_Zuiderzeeland.xlsx'
admins_shp = r'..\..\data\sources\admins.shp'
profiel_shp = r'..\..\data\sources\Zuiderzeeland\Profiel.shp'
watergangen_shp = r'..\..\data\sources\Zuiderzeeland\Watergangen.shp'
peilgebieden_shp = r'..\..\data\sources\Zuiderzeeland\Peilgebieden Besluit.shp'
profielen_result_shp = r'..\..\data\sources\Zuiderzeeland\Profiel_merged_{}.shp'.format(datetime.datetime.now().strftime('%Y%m%d'))
dem_tif = r'..\..\data\sources\Zuiderzeeland\zuiderzeeland_dtm.tif'
get_dem = False
cell_size = 5

#%% inlezen poly van waterschap

admins_gdf = gpd.read_file(admins_shp)
poly = admins_gdf.loc[admins_gdf['id'] == 37]['geometry'].values[0]

if get_dem:
    get_ahn.to_tif(poly,dem_tif,cell_size=cell_size)


#%% inlezen legger tabellen
vaarten_df = pd.read_excel(legger_xlsx,
                           header=[0, 1, 2],
                           sheet_name='Vaarten&Tochten')
cols = [' '.join([cel for cel in col if not 'Unnamed' in cel]) for col in list(vaarten_df.columns)]
vaarten_df.columns = cols


stedelijk_owt_df = pd.read_excel(legger_xlsx,
                           header=[0, 1, 2],
                           sheet_name='Stedelijk_OWT')
cols = [' '.join([cel for cel in col if not 'Unnamed' in cel]) for col in list(stedelijk_owt_df.columns)]
stedelijk_owt_df.columns = cols

stedelijk_bwt_df = pd.read_excel(legger_xlsx,
                           header=[0, 1, 2],
                           sheet_name='Stedelijk_BWT')
cols = [' '.join([cel for cel in col if not 'Unnamed' in cel]) for col in list(stedelijk_bwt_df.columns)]
stedelijk_bwt_df.columns = cols

sloten_df = pd.read_excel(legger_xlsx,
                           sheet_name='Sloten')

#%% inlezen profielen & watergangen
profiel_gdf = gpd.read_file(profiel_shp)
watergangen_gdf = gpd.read_file(watergangen_shp)

#%% inlezen peilvakken
peilgebieden_gdf = gpd.read_file(peilgebieden_shp)

#%% open data_array for zonal_stats
with rasterio.open(dem_tif) as src:
    scale = src.scales[0]
    raster_profile = src.profile
    raster_data = src.read(1) * scale
    affine = raster_profile['transform']


#%% toevoegen Naam, BL, WD, BB, CAT & Z aan profielen
def add_attributes(row):
    ident = row['OWA_OWA_ID']
    name = BL = WD = BB = CAT = Z = ZP = SUB = np.NaN
    
    #add name
    name_array = watergangen_gdf.loc[watergangen_gdf['GLOBALID'] == ident]['OWANAAM'].values
    if len(name_array) > 0:
        if not name_array[0] == None:
            name = name_array[0]
            
    #add zomerpeil
    peil_gdf = peilgebieden_gdf[peilgebieden_gdf.intersects(row['geometry'])]
    
    if len(peil_gdf) > 0:
        if not peil_gdf['GPGZMRPL'].values[0] in [0,-999]:
            ZP = peil_gdf['GPGZMRPL'].values[0]
        elif not peil_gdf['IWS_ONDERG'].values[0] in [0,-999]:
            ZP = peil_gdf['IWS_ONDERG'].values[0]
        elif not peil_gdf['IWS_BOVENG'].values[0] in [0,-999]:
            ZP = peil_gdf['IWS_BOVENG'].values[0] 
    
    Z = zonal_stats(row['geometry'].buffer(50),
                    raster_data,
                    affine=affine,
                    nodata=raster_profile['nodata'],
                    stats='percentile_50',raster_out=False)[0]['percentile_50']        
            
    profiel_ident = row['OSMOMSCH']
    if type(profiel_ident) == str:
        if not vaarten_df[vaarten_df['OMSCHRIJVING'].str.contains(profiel_ident,regex=False)].empty:
            profiel = vaarten_df[vaarten_df['OMSCHRIJVING'].str.contains(profiel_ident,regex=False)].to_dict('records')[0]
            BL = ZP - profiel['d Bodemdiepte']
            BB = profiel['b Bodembreedte']
            CAT = 1 #vaarten zijn altijd primair
            SUB = 1
            
        elif not sloten_df[sloten_df['OMSCHRIJVING'].str.contains(profiel_ident,regex=False)].empty:
            profiel = sloten_df[sloten_df['OMSCHRIJVING'].str.contains(profiel_ident,regex=False)].to_dict('records')[0]

            
            if profiel['Profielhoogte'] == profiel['Profielhoogte']:
                BL = Z - profiel['Profielhoogte']
                
            elif profiel_ident in ['K-01','K-02']:
                BL = ZP + 0.2 # uit tabellen 'De bodem van de kavelsloot moet ca 0-20cm boven tochtpeil liggen
            
            BB = profiel['Bodembreedte']
            
            if BL < ZP:
                SUB = 1
                CAT = 2
            else:
                SUB = 0
                CAT = 3 # droge sloot, want BL < ZP
            
        elif not stedelijk_owt_df[stedelijk_owt_df['OMSCHRIJVING'].str.contains(profiel_ident,regex=False)].empty:
            profiel = stedelijk_owt_df[stedelijk_owt_df['OMSCHRIJVING'].str.contains(profiel_ident,regex=False)].to_dict('records')[0]
            
            if profiel['d Bodemdiepte'] == profiel['d Bodemdiepte']:
                BL = ZP - profiel['d Bodemdiepte']
                
                if BL < ZP:
                    SUB = 1
                    CAT = 2
                else:
                    SUB = 0
                    CAT = 3 # droge sloot, want BL < ZP
            
    return name, BL, WD, BB, CAT, Z, ZP, SUB

profiel_gdf[['name', 'BL', 'WD', 'BB', 'CAT', 'Z', 'ZP', 'SUB']] = profiel_gdf.apply(add_attributes, 
                                                                   axis=1, 
                                                                   result_type="expand")




drop_cols = [col for col in profiel_gdf.columns if not col in ['name','Z','OSMOMSCH','OWA_OWA_ID', 'BL', 'WD', 'BB', 'CAT','Z','ZP', 'SUB', 'geometry']]
profiel_gdf = profiel_gdf.drop(drop_cols,axis=1)
profiel_gdf.to_file(profielen_result_shp)
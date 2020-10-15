# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:30:34 2020

@author: danie
"""
import pandas as pd
import fiona
import gdal
import geopandas as gpd
import services
import warnings
import requests
warnings.simplefilter(action='ignore', category=FutureWarning)

#brondata
waterschappen_gml = r'..\\data\\sources\\Waterschapsgrenzen_4258.gml'
landsgrens_gml = r'..\\data\\sources\\Landsgrens.gml'
rws_url = 'https://geoservices.rijkswaterstaat.nl/rws_legger_2.0'

#resultaat-bestand
admin_shp = r'..\\data\\sources\\admins.shp'


def wfs2gdf(url):
    response = requests.get(url)
    
    try:
        vsif = fiona.ogrext.buffer_to_virtual_file(bytes(response.content))
        vsiz = vsif+'.gml'
        gdal.Rename(vsif, vsiz)
        fc = fiona.Collection(vsiz)
        
        return gpd.GeoDataFrame.from_features([feature for feature in fc], crs='epsg:28992')
    except:
        print(response.text)

#%% inlezen landsgrens

url = ('https://geoservices.rijkswaterstaat.nl/rws_legger_2.0?&'
       'service=WFS&'
       'version=1.1.0&'
       'request=GetFeature&'
       'typeName={layer}').format(layer='basis_nl')


landsgrens_gdf = wfs2gdf(url)
landsgrens_poly = landsgrens_gdf.dissolve(by='NAAM').loc['Nederland'].geometry

#%% inlezen begrenzing rijksvaarweg
url = ('https://geoservices.rijkswaterstaat.nl/rws_legger_2.0?&'
       'service=WFS&'
       'version=1.1.0&'
       'request=GetFeature&'
       'typeName={layer}').format(layer='oppervlwaterlich')

rijksvaarweg_gdf = wfs2gdf(url)
rijksvaarweg_gdf['dissolve'] = 0
rijksvaarweg_poly = rijksvaarweg_gdf.dissolve(by='dissolve').geometry[0]

#%% inlezen begrenzing waterschappen
waterschappen_gdf = gpd.read_file(waterschappen_gml)
waterschappen_gdf = waterschappen_gdf.to_crs('epsg:28992')

#%% maken van een negatief van de rijkswaarwegen
rijksvaarweg_negatief_poly = landsgrens_poly.difference(rijksvaarweg_poly)

#%% knippen van waterschaps admins op rws_negatief
waterschappen_gdf = gpd.clip(waterschappen_gdf[waterschappen_gdf['geometry'].is_valid], rijksvaarweg_negatief_poly)

#waterschappen_gdf = gpd.clip(waterschappen_gdf, rijksvaarweg_negatief_poly)


#%% samenvoegen tot 1 admin dataframe
rws_gdf = gpd.GeoDataFrame(data={'name':['Rijkswaterstaat'],
                                  'id':[1],
                                  'geometry':[rijksvaarweg_poly]})

waterschappen_gdf.rename(columns={'waterschapsnaam':'name',
                                  'publicerende_waterbeheerder':'id'},
           inplace=True)

drop_cols = [col for col in waterschappen_gdf.columns if not col in ['name','id','geometry']]
waterschappen_gdf = waterschappen_gdf.drop(drop_cols,axis=1)

admins_gdf = gpd.GeoDataFrame(pd.concat([waterschappen_gdf,rws_gdf], ignore_index=True) )
admins_gdf.crs = 'epsg:28992'
admins_gdf.to_file(admin_shp)

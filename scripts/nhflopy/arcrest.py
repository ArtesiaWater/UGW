# -*- coding: utf-8 -*-

import geopandas as gpd
import pandas as pd
import requests


class ArcREST:
    def __init__(self, url, crs='epsg:28992', maxRecordCount=2000):
        epsg = crs[crs.find(':') + 1:]
        self.maxRecordCount = requests.get('{}?f=pjson'.format(url)).json()[
            'maxRecordCount']
        self.maxRecordCount = min(self.maxRecordCount, maxRecordCount)
        self.url = url
        self.epsg = epsg
        self.crs = crs

    def get_features(self, layer, poly, object_filter=''):
        ''' to download features from a layer for a shapely polygon 

            Parameters:
            layer: integer layer number
            poly: shapely Polygon object used as a boundary
        '''

        if hasattr(poly, 'bounds'):
            xmin, ymin, xmax, ymax = poly.bounds
        else:
            xmin, xmax, ymin, ymax = poly

        if not object_filter == '':
            object_filter = ' and {}'.format(object_filter)

        url = ('{url}/{layer}/query?'
               'where=OBJECTID>=0{object_filter}'
               '&f=json'
               '&geometry={xmin},{ymin},{xmax},{ymax}'
               '&geometryType=esriGeometryEnvelope'
               '&inSR={epsg}'
               '&returnIdsOnly=true').format(url=self.url,
                                             layer=layer,
                                             xmin=xmin,
                                             ymin=ymin,
                                             xmax=xmax,
                                             ymax=ymax,
                                             epsg=self.epsg,
                                             object_filter=object_filter)

        response = requests.get(url)
        if response.status_code == 200:
            if 'objectIds' in list(response.json().keys()):
                object_ids = response.json()['objectIds']
            else:
                object_ids = response.json()['properties']['objectIds']
            object_ids.sort()

            downloads = round(len(object_ids) / self.maxRecordCount + 0.5)
            gdf_list = []
            for download in range(downloads):
                min_object = download * self.maxRecordCount
                max_object = min(
                    min_object + self.maxRecordCount - 1, len(object_ids) - 1)
                url = ('{url}/{layer}/query?'
                       'where={min_objects}<=OBJECTID and {max_objects}>=OBJECTID{object_filter}'
                       '&outFields=*'
                       '&geometry={xmin},{ymin},{xmax},{ymax}'
                       '&geometryType=esriGeometryEnvelope'
                       '&inSR={epsg}'
                       '&outSR={epsg}&f=geojson').format(
                           url=self.url,
                           layer=layer,
                           min_objects=object_ids[min_object],
                           max_objects=object_ids[max_object],
                           object_filter=object_filter,
                           xmin=xmin,
                           ymin=ymin,
                           xmax=xmax,
                           ymax=ymax,
                           epsg=self.epsg)
                response = requests.post(url)
                features = response.json()['features']
                gdf = gpd.GeoDataFrame.from_features(features)
                gdf.crs = self.crs
                if hasattr(poly, 'bounds'):
                    gdf = gdf[gdf.intersects(poly)]
                gdf_list += [gdf]

            if len(gdf_list) > 1:
                gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True))
            else:
                gdf = gdf_list[0]

            return gdf

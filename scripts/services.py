import logging

import fiona
import gdal
import geopandas as gpd
import pandas as pd
import requests


class ArcREST:

    def __init__(self, url, output_format='geojson', crs='epsg:28992', maxRecordCount=2000):
        epsg = crs[crs.find(':') + 1:]
        self.maxRecordCount = min(requests.get('{}?f=pjson'.format(url)).json()['maxRecordCount'],
                                  maxRecordCount)
        self.url = url
        self.epsg = epsg
        self.crs = crs
        self.format = output_format

    def get_features(self, layer, poly, object_filter=''):
        ''' to download features from a layer for a shapely polygon 

            Parameters:
            layer: integer layer number
            poly: shapely Polygon object used as a boundary
        '''

        xmin, ymin, xmax, ymax = poly.bounds
        try:
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
                           '&outSR={epsg}&f={output_format}').format(url=self.url,
                                                                     layer=layer,
                                                                     min_objects=object_ids[min_object],
                                                                     max_objects=object_ids[max_object],
                                                                     object_filter=object_filter,
                                                                     xmin=xmin,
                                                                     ymin=ymin,
                                                                     xmax=xmax,
                                                                     ymax=ymax,
                                                                     epsg=self.epsg,
                                                                     output_format=self.format)
                    response = requests.post(url)
                    if self.format == 'json':
                        logging.warning(
                            'reading ESRI-json format (GeoJSON is preferred)')
                        vsif = fiona.ogrext.buffer_to_virtual_file(
                            bytes(response.content))
                        vsiz = vsif + '.json'
                        gdal.Rename(vsif, vsiz)
                        fc = fiona.Collection(vsiz)
                        gdf = gpd.GeoDataFrame.from_features(
                            [feature for feature in fc], crs='epsg:28992')
                        columns = [
                            col for col in gdf.columns if not col == 'geometry'] + ['geometry']
                        gdf = gdf[columns]
                    else:
                        features = response.json()['features']
                        gdf = gpd.GeoDataFrame.from_features(features)
                        gdf.crs = self.crs
                    gdf = gdf[gdf.intersects(poly)]
                    gdf_list += [gdf]

                if len(gdf_list) > 1:
                    gdf = gpd.GeoDataFrame(
                        pd.concat(gdf_list, ignore_index=True))

                else:
                    gdf = gdf_list[0]

                return(gdf)

        except:
            logging.error(
                'Processing data from the following url failed: {}'.format(url))


class WFS:

    def __init__(self, url, crs='epsg:28992'):
        self.wfs = WebFeatureService(url, version='2.0.0', timeout=300)
        self.layers = list(self.wfs.contents.keys())
        self.crs = 'epsg:28992'

    def get_features(self, layer, poly):
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
            logging.error('{} not in layers({})'.format(layer, self.layers))

            return None

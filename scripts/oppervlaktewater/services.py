# -*- coding: utf-8 -*-

import logging
import sys
import xml.etree.ElementTree as ET
from io import StringIO

import fiona
import gdal
import geopandas as gpd
import pandas as pd
import pyproj
import requests


def wfs2gdf(url, crs, output_type, timeout):
    response = requests.get(url, timeout=timeout)
    try:
        if output_type == 'application/json':
            features = response.json()['features']
            return gpd.GeoDataFrame.from_features(features)
        else:
            vsif = fiona.ogrext.buffer_to_virtual_file(bytes(response.content))
            vsiz = vsif + '.gml'
            gdal.Rename(vsif, vsiz)
            fc = fiona.Collection(vsiz)

            return gpd.GeoDataFrame.from_features([feature for feature in fc], crs=crs)
    except:
        print(response.text)

# %%


def set_crs(gdf, crs='epsg:28992'):
    if gdf.crs == None:
        gdf.crs = crs
    else:
        if hasattr(pyproj, 'CRS'):
            update_crs = not pyproj.CRS(gdf.crs).equals(pyproj.CRS(crs))
        else:
            update_crs = pyproj.Proj(gdf.crs).srs != pyproj.Proj(init=crs).srs
        if update_crs:
            gdf = gdf.to_crs({'init': crs})

    return gdf


class ArcREST:

    def __init__(self, url, output_format='geojson', crs='epsg:28992', maxRecordCount=2000):
        epsg = crs[crs.find(':') + 1:]

        properties = requests.get('{}?f=pjson'.format(url)).json()
        self.maxRecordCount = min(properties['maxRecordCount'], maxRecordCount)
        self.layers = properties['layers']
        self.url = url
        self.epsg = epsg
        self.crs = crs
        self.format = output_format

    def get_features(self, layer, poly=None, object_filter='', object_id=None):
        ''' to download features from a layer for a shapely polygon 

            Parameters:
            layer: integer layer number
            poly: shapely Polygon object used as a boundary
        '''

        if object_id == None:
            properties = requests.get(
                ('{url}/{layer}/?f=pjson').format(url=self.url, layer=layer)).json()

            if "uniqueIdField" in list(properties.keys()):
                object_id = properties["uniqueIdField"]["name"]
            else:
                if "fields" in list(properties.keys()):
                    field = [field['name'] for field in properties["fields"]
                             if field['name'].lower() == "objectid"]
                    if len(field) == 1:
                        object_id = field[0]
            if object_id == None:
                logging.error(
                    'Processing data from the following url failed: {url}/{layer}/?f=pjson'.format(url=self.url, layer=layer))
                logging.error(('ArcREST Layer has no Unique ID Field, script defaulted to {object_id}.'
                               'Please specify a correct object_id for this layer & adminstration').format(object_id=object_id))
                sys.exit()

        xmin, ymin, xmax, ymax = poly.bounds
        try:
            if not object_filter == '':
                object_filter = ' and {}'.format(object_filter)
            url = ('{url}/{layer}/query?'
                   'where={object_id}>=0{object_filter}'
                   '&geometry={xmin},{ymin},{xmax},{ymax}'
                   '&geometryType=esriGeometryEnvelope'
                   '&f=json'
                   '&inSR={epsg}'
                   '&returnIdsOnly=true').format(url=self.url,
                                                 layer=layer,
                                                 object_id=object_id,
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
                if object_ids == None:
                    return gpd.GeoDataFrame()
                else:
                    object_ids.sort()

                    downloads = round(len(object_ids) /
                                      self.maxRecordCount + 0.5)
                    gdf_list = []
                    for download in range(downloads):
                        min_object = download * self.maxRecordCount
                        max_object = min(
                            min_object + self.maxRecordCount - 1, len(object_ids) - 1)
                        url = ('{url}/{layer}/query?'
                               'where={min_objects}<={object_id} and {max_objects}>={object_id}{object_filter}'
                               '&outFields=*'
                               '&geometry={xmin},{ymin},{xmax},{ymax}'
                               '&geometryType=esriGeometryEnvelope'
                               '&inSR={epsg}'
                               '&outSR={epsg}&f={output_format}').format(url=self.url,
                                                                         layer=layer,
                                                                         object_id=object_id,
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
                        gdf = gpd.GeoDataFrame()
                        if len(response.json()['features']) > 0:
                            if self.format == 'json':
                                logging.warning(
                                    'reading ESRI-json format (GeoJSON is preferred)')
                                vsif = fiona.ogrext.buffer_to_virtual_file(
                                    bytes(response.content))
                                vsiz = vsif + '.json'
                                gdal.Rename(vsif, vsiz)
                                fc = fiona.Collection(vsiz)
                                gdf = gpd.GeoDataFrame.from_features(
                                    [feature for feature in fc], crs=self.crs)
                                columns = [
                                    col for col in gdf.columns if not col == 'geometry'] + ['geometry']
                                gdf = gdf[columns]
                            else:
                                features = response.json()['features']
                                gdf = gpd.GeoDataFrame.from_features(features)
                        else:
                            logging.warning(
                                'no features returned for url: {}'.format(url))

                        if len(gdf) > 0:
                            gdf.crs = self.crs
                            gdf = gdf[gdf.intersects(poly)]
                            gdf_list += [gdf]

                    if not gdf.empty:
                        if len(gdf_list) > 1:
                            gdf = gpd.GeoDataFrame(
                                pd.concat(gdf_list, ignore_index=True))

                        else:
                            gdf = gdf_list[0]

                    layer_name = [lay['name']
                                  for lay in self.layers if lay['id'] == layer][0]
                    gdf['layer_name'] = layer_name

                    return(gdf)

        except Exception as e:
            logging.error(
                'Processing data from the following url failed: {} with error {}'.format(url, e))


class WFS:

    def __init__(self, url, timeout=300):
        self.url = url
        self.timeout = timeout
        response = requests.get('{url}?service=WFS&request=GetCapabilities'.format(url=self.url),
                                timeout=self.timeout)

        self.capabilities_ns = dict([node for _,
                                     node in ET.iterparse(StringIO(response.text),
                                                          events=['start-ns'])])

        self.capabilities = ET.ElementTree(
            ET.fromstring(response.text)).getroot()

        self.version = self.capabilities.attrib['version']

        self.feature_type_list = self.capabilities.find('wfs:FeatureTypeList',
                                                        self.capabilities_ns)
        self.layers = [feature_type.find('wfs:Name', self.capabilities_ns).text
                       for feature_type
                       in self.feature_type_list.findall('wfs:FeatureType', self.capabilities_ns)]

        operations = self.capabilities.find('ows:OperationsMetadata', self.capabilities_ns).findall('ows:Operation',
                                                                                                    self.capabilities_ns)

        parameters = [operation for operation in operations if
                      operation.attrib['name'] == 'GetFeature'][0].findall('ows:Parameter',
                                                                           self.capabilities_ns)

        allowed_values = [parameter for parameter in parameters if
                          parameter.attrib['name'] == 'outputFormat'][0].find('ows:AllowedValues',
                                                                              self.capabilities_ns)

        if not allowed_values == None:
            self.output_formats = [value.text for
                                   value in allowed_values.findall('ows:Value',
                                                                   self.capabilities_ns)]
        else:
            self.output_formats = [[parameter for parameter in parameters if
                                    parameter.attrib['name'] == 'outputFormat'][0].find('ows:Value',
                                                                                        self.capabilities_ns).text]

        self.output_format = None
        if 'application/json' in self.output_formats:
            self.output_format = 'application/json'
        else:
            gml_formats = [
                output_format for output_format in self.output_formats if 'text/xml' in output_format]
            if len(gml_formats) > 0:
                self.output_format = gml_formats[0]

        if self.output_format == None:
            logging.warning(
                'no suitable output format available from {}'.format(self.output_formats))
        else:
            logging.info('output format set to {}'.format(self.output_format))

    def get_features(self, layer, crs, poly=None):
        if layer in self.layers:
            feature_type = self.feature_type_list[self.layers.index(layer)]

            src_crs = feature_type.find('wfs:DefaultSRS', self.capabilities_ns)
            if src_crs == None:
                src_crs = feature_type.find(
                    'wfs:DefaultCRS', self.capabilities_ns)

            src_crs = src_crs.text
            src_crs = src_crs[src_crs.find('EPSG'):].replace('::', ':')

            url = ('{url}?&'
                   'service=WFS&'
                   'version={version}&'
                   'request=GetFeature&'
                   'typeName={layer}&'
                   'outputFormat={output_format}').format(url=self.url,
                                                          version=self.version,
                                                          layer=layer,
                                                          output_format=self.output_format)

            if not poly == None:
                bbox = ",".join(map(str, list(poly.bounds)))
                url = '{url}&bbox={bbox}'.format(url=url, bbox=bbox)

            gdf = wfs2gdf(
                url, src_crs, output_type=self.output_format, timeout=self.timeout)
            gdf = set_crs(gdf, crs=crs)

            if not poly == None and len(gdf) > 0:
                gdf = gdf[gdf.intersects(poly)]

            return gdf
        else:
            logging.error('{} not in layers({})'.format(layer, self.layers))

            return None

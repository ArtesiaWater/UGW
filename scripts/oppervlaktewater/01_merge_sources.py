# -*- coding: utf-8 -*-

# public packages
import configparser
import json
import logging
import os
import shutil
import sys
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
from shapely.geometry import Polygon

import get_ahn
import get_bgt
# project scripts
import services

start = time.time()
bgt_service = 'website'

# %% alle paden van output bestanden
config_ini = configparser.ConfigParser()
config_ini.read(r'config.ini')
project = config_ini['general']['project']

# %%
extent_shp = 'extent.shp'
if 'extent' in config_ini['general'].keys():
    extent_shp = config_ini['general']['extent']

os.chdir('../../config')
project_shp = Path(f'../data/{project}/{extent_shp}')
input_dir = Path(f'../data/{project}/input')

admins = {'file_name': Path(r'../config/administrations.json')}
sources = {'file_name': Path(r'../config/sources.json')}
admins.update(json.loads(open(admins['file_name'], 'r').read()))
sources.update(json.loads(open(sources['file_name'], 'r').read()))

admins_shp = input_dir.joinpath('waterschappen.shp')

if input_dir.exists():
    shutil.rmtree(input_dir)
input_dir.mkdir(parents=True)

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

src_layers = ['water-level_areas', 'water_areas', 'water_lines']

get_dem = True
cell_size = 5

dem_tif = input_dir.joinpath(f'ahn3_{int(cell_size)}m_dtm.tif')

# %% alle functies die worden aangeroepen


def set_crs(gdf, crs=sources['default_crs']):
    '''sets the coordinate reference system to a gdf if not specified'''

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


def to_lineString(gdf):
    '''converts a mixed MultiLineString + LineString polyline gdf to single lines'''

    outgdf = gdf[gdf.geometry.type == 'LineString']

    for index, row in gdf[gdf.geometry.type == 'MultiLineString'].iterrows():
        multdf = gpd.GeoDataFrame(columns=gdf.columns)
        recs = len(row.geometry)
        multdf = multdf.append([row] * recs, ignore_index=True)
        for geom in range(recs):
            multdf.loc[geom, 'geometry'] = row.geometry[geom]
        outgdf = outgdf.append(multdf, ignore_index=True)

    return outgdf


def MultiPolgyon2Polygon(gdf):
    outgdf = gpd.GeoDataFrame(columns=gdf.columns)
    for idx, row in gdf.iterrows():
        multdf = gpd.GeoDataFrame(columns=gdf.columns)
        recs = len(row.geometry)
        multdf = multdf.append([row] * recs, ignore_index=True)
        for geom in range(recs):
            multdf.loc[geom, 'geometry'] = row.geometry[geom]
        outgdf = outgdf.append(multdf, ignore_index=True)
    return outgdf


def explode(gdf):
    ''' explodes Multipolygons to Polygons '''
    outgdf = gdf[gdf.geometry.type == "Polygon"]
    outgdf = outgdf.append(MultiPolgyon2Polygon(
        gdf[gdf.geometry.type == "MultiPolygon"]), ignore_index=True)
    for idx, row in gdf[gdf.geometry.type == "GeometryCollection"].iterrows():
        geoms = [geom for geom in row['geometry']
                 if geom.type in ['Polygon', 'MultiPolygon']]
        for geom in geoms:
            if geom.type == 'Polygon':
                outgdf.append(row, ignore_index=True)
                outgdf.tail(1)['geometry'] = geom
            elif geom.type == 'MultiPolygon':
                multdf = gpd.GeoDataFrame(columns=gdf.columns)
                recs = len(geom)
                multdf = multdf.append([row] * recs, ignore_index=True)
                for rec in range(recs):
                    multdf.loc[rec, 'geometry'] = geom[rec]
                outgdf = outgdf.append(multdf, ignore_index=True)

    return outgdf


admin_conversion = {value['bronhouder']:
                    {'name': key,
                     'id': admins['administrations'][key]}
                    for key, value in sources['water_areas'].items()
                    if value['source'] == 'bgt'}


def bgt_admin_conversion(row):
    ''' apply function to specify admin_id and admin_name from bronhouder or íntersection'''
    admin_id = admin_name = None
    if row['bronhouder'] in admin_conversion.keys():
        admin_id = admin_conversion[row['bronhouder']]['id']
        admin_name = admin_conversion[row['bronhouder']]['name']

    else:
        gdf = admins_gdf[admins_gdf.intersects(row['geometry'])]
        if not gdf.empty:
            admin_id = gdf['id'].values[0]
            admin_name = gdf['name'].values[0]

    return admin_id, admin_name


def get_layer(src_layer, admins_gdf):
    ''' get the layer (water_lines, water_areas or water-level_areas) in gdf '''
    gdf_list = []
    gdf_invalid_list = []
    bgt_names = []

    if src_layer == 'water_areas':
        bgt_gdf = admins_gdf[admins_gdf['name'].isin(
            [key for key, value in sources[src_layer].items() if value['source'] == 'bgt'])]
        bgt_names = list(bgt_gdf['name'].values)
        if len(bgt_gdf) > 0:
            logging.info('bgt for {}'.format(
                ", ".join(list(bgt_gdf['name'].values))))
            bgt_gdf['dissolve'] = 0
            logging.info('create exterior polygon for bgt-download')
            bgt_gdf = gpd.clip(bgt_gdf, project_mask)
            xmin = bgt_gdf['geometry'].bounds['minx'].values.min()
            xmax = bgt_gdf['geometry'].bounds['maxx'].values.max()
            ymin = bgt_gdf['geometry'].bounds['miny'].values.min()
            ymax = bgt_gdf['geometry'].bounds['maxy'].values.max()

            poly = Polygon([[xmin, ymax], [xmax, ymax], [
                           xmax, ymin], [xmin, ymin], [xmin, ymax]])

            bronhouders = [value['bronhouder'] for key, value in sources[src_layer].items(
            ) if key in bgt_gdf['name'].values]

            logging.info('starting bgt download')
            gdf = get_bgt.to_gdf(poly, layer="waterdeel",
                                 bronhouders=None,
                                 end_registration='now',
                                 log_level="INFO",
                                 service=bgt_service)

            valid_bronhouders = [bronhouder for bronhouder in np.unique(
                gdf['bronhouder'].values) if (bronhouder in bronhouders) | (not bronhouder[0] == 'W')]
            gdf = gdf[gdf['bronhouder'].isin(valid_bronhouders)]

            logging.info('clip bgt result to polygon exterior')
            gdf = gpd.clip(gdf, poly)

            logging.info('assigning admin names and ids')
            gdf[['admin_id', 'admin']] = gdf.apply(
                bgt_admin_conversion, axis=1, result_type="expand")

            drop_cols = [col for col in list(gdf.columns) if not col in [
                'gml_id', 'admin_id', 'admin', 'geometry']]
            gdf = gdf.drop(drop_cols, axis=1)
            gdf.rename(columns={'gml_id': 'src_id'},
                       inplace=True)
            gdf.insert(0, 'src', 'bgt')

            gdf_list += [gdf]

            admins_gdf = admins_gdf[admins_gdf['name'].isin(
                [key for key, value in sources[src_layer].items() if value['source'] != 'bgt'])]

    for index, row in admins_gdf.iterrows():
        admin = row['name']
        admin_id = row['id']
        if not admin in bgt_names:
            logging.info('{} {}'.format(src_layer, admin))

        if admin in list(sources[src_layer].keys()):
            src = sources[src_layer][admin]
            poly = row['geometry'].intersection(project_mask)

            if src['source'] in list(sources['services'].keys()):
                src_url = sources['services'][src['source']]
                url = src_url['url']
                if 'layer' in list(src.keys()):
                    layers = [src['layer']]
                elif 'layers' in list(src.keys()):
                    layers = src['layers']

                gdf_tmp_list = []
                for layer in layers:
                    if src_url['type'] == 'wfs':

                        wfs = services.WFS(url)
                        if layer in wfs.layers:
                            gdf = wfs.get_features(
                                layer, crs=sources['default_crs'], poly=poly)
                        else:
                            logging.error('{} not a layer in {}'.format(
                                layer, src['source']))

                    elif src_url['type'] == 'arcrest':
                        output_format = 'geojson'

                        if 'format' in list(src_url.keys()):
                            output_format = src_url['format']
                        object_filter = ''
                        if 'filter' in list(src.keys()):
                            object_filter = src['filter']
                        rest = services.ArcREST(
                            url, output_format=output_format)
                        if 'object_id' in list(src.keys()):
                            object_id = src['object_id']
                        else:
                            object_id = None

                        gdf = rest.get_features(
                            layer, poly, object_filter=object_filter, object_id=object_id)

                    elif src_url['type'] == 'gdb':
                        gdf = gpd.read_file(
                            url, driver='FileGDB', layer=layer, bbox=poly)

                    else:
                        logging.error(
                            "type '{}' not supported".format(src_url['type']))

                    gdf_tmp_list += [gdf]

                if len(gdf_tmp_list) == 1:
                    gdf = gdf_tmp_list[0]
                else:
                    gdf = gpd.GeoDataFrame(
                        pd.concat(gdf_tmp_list, ignore_index=True))

            elif src['source'] in list(sources['files'].keys()):
                src_path = sources['files'][src['source']]['path']
                gdf = gpd.read_file(src_path)
                gdf = set_crs(gdf)
                gdf = gdf[gdf.intersects(poly)]

            if len(gdf) > 0:
                if 'attributes' in list(src.keys()):

                    attributes = src['attributes']

                    merge_attributes = {
                        key: item for key, item in attributes.items() if type(item) == list}

                    if len(merge_attributes) > 0:
                        gdf.fillna(value=np.NaN, inplace=True)
                        for attrib, attrib_list in merge_attributes.items():
                            gdf[attrib] = np.NaN
                            for attrib_item in attrib_list:
                                gdf[attrib] = gdf.apply((lambda x: x[attrib_item] if
                                                         x[attrib_item] == x[attrib_item] else np.NaN),
                                                        axis=1)

                        attributes = {
                            key: item for key, item in attributes.items() if not type(item) == list}

                    drop_cols = [col for col in list(gdf.columns) if not col in list(
                        attributes.values()) + ['geometry']]
                    gdf = gdf.drop(drop_cols, axis=1)
                    gdf.rename(columns={value: key for key, value in attributes.items()},
                               inplace=True)
                else:
                    drop_cols = [col for col in list(
                        gdf.columns) if col != 'geometry']
                    gdf = gdf.drop(drop_cols, axis=1)

                if 'values' in list(src.keys()):
                    values = src['values']
                    for ident, val in values.items():
                        gdf[ident] = val

                gdf.insert(0, 'src', src['source'])
                gdf.insert(0, 'admin_id', admin_id)
                gdf.insert(0, 'admin', admin)
                gdf_invalid = gdf[~gdf['geometry'].is_valid]

                if (not len(gdf_invalid)) == 0 & all(geom_type in ['Polygon', 'MultiPolygon'] for geom_type in np.unique(gdf['geometry'].type)):
                    logging.warning(
                        'trying to fix geometry errors by 0m buffer')
                    gdf['geometry'] = gdf['geometry'].buffer(0)

                gdf_invalid_list += [gdf[~gdf['geometry'].is_valid]]
                gdf = gdf[gdf['geometry'].is_valid]
                gdf = gpd.clip(gdf, poly)
                gdf_list += [gdf]
            else:
                logging.warning(
                    'data-source does not contain features for this project-extent')

        else:
            logging.warning(
                'source for {} not specified in {}'.format(admin, src_layer))

    # merge all into one dataframe
    gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True))
    gdf.crs = sources['default_crs']
    gdf.fillna(value=np.NaN, inplace=True)

    # merge all invalid into a dataframe
    if len(gdf_invalid_list) > 0:
        gdf_invalid = gpd.GeoDataFrame(
            pd.concat(gdf_invalid_list, ignore_index=True))
        gdf_invalid.crs = sources['default_crs']
        gdf_invalid.fillna(value=np.NaN, inplace=True)
    else:
        gdf_invalid = gpd.GeoDataFrame()

    return gdf, gdf_invalid


def schema(gdf):
    ''' define schema for gdf '''
    cols = [col for col in gdf.columns if not col == 'geometry']
    schema = {'geometry': min(
        list(np.unique(gdf['geometry'].type.values)), key=len)}

    schema['properties'] = {}

    for col in cols:
        if str in [type(elem) for elem in list(gdf[col].values)]:
            schema['properties'].update({col: 'str'})
        else:
            schema['properties'].update({col: 'float'})
    return schema


def check_sources():
    ''' check if sources are available '''
    logging.info('checking config for errors')
    for index, row in admins_gdf.iterrows():
        admin = row['name']
        for src_layer in src_layers:
            if (not admin in list(sources[src_layer])) & (src_layer == 'water_areas'):
                logging.error('"{}" missing in layer "{}" in config-file "{}"'.format(
                    admin,
                    src_layer,
                    os.path.abspath(sources['file_name'])))
                sys.exit()


# %% aanmaken van een project-mask
logging.info('Dissolving project mask')
project_gdf = gpd.read_file(project_shp)
project_gdf['dissolve'] = 0
project_mask = project_gdf.dissolve(by='dissolve').geometry[0]

# %% clippen van de admin shape-file naar project-mask
admin_boundaries = admins['boundaries']
admins_gdf = gpd.read_file(admin_boundaries['path'])
admins_gdf = set_crs(admins_gdf)
admins_gdf = admins_gdf[admins_gdf.geometry.intersects(project_mask)]
id_field = admin_boundaries['id_field']
admins_gdf = admins_gdf[admins_gdf[id_field].isin(
    [value for key, value in admins['administrations'].items()])]
admins_gdf = admins_gdf.drop([col for col in admins_gdf.columns if not col in [
                             id_field, 'geometry']], axis=1)
admins_gdf.columns = ['id', 'geometry']
admins_gdf['name'] = admins_gdf['id'].apply(lambda x: list(
    admins['administrations'].keys())[list(admins['administrations'].values()).index(x)])
admins_gdf.to_file(admins_shp)

# %% controleren of bronnen beschikbaar zijn
check_sources()

# %% inlezen van alle bronnen
for src_layer in src_layers:
    gdf, gdf_invalid = get_layer(src_layer, admins_gdf)
    if 'GeometryCollection' in np.unique(gdf.geometry.type):
        gdf = explode(gdf)
    if not gdf.empty:
        gdf.to_file(input_dir.joinpath(
            '{}.shp'.format(src_layer)), schema=schema(gdf))
    if not gdf_invalid.empty:
        gdf_invalid.to_file(input_dir.joinpath(
            '{}_invalid.shp'.format(src_layer)), schema=schema(gdf))

# %% get sources: water areas (ahn)
if get_dem:
    get_ahn.to_tif(project_mask, dem_tif, cell_size=cell_size)


logging.info('done in {} seconds'.format(round(time.time() - start, 1)))

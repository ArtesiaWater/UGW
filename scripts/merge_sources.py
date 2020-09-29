# public packages
import json
import logging
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj

import get_ahn
import get_bgt

# project scripts
import services

'''
ToDo:
    - values uit get_layer slopen en in prepare_modflow zetten
    - explode naar prepare_modflow
'''

# %%
admins = json.loads(open('../config/administrations.json', 'r').read())
sources = json.loads(open('../config/sources.json', 'r').read())
project_shp = '../data/input/project.shp'

admins_shp = '../data/input/waterschappen.shp'
wl_areas_shp = '../data/input/peilvakken.shp'
w_areas_shp = '../data/input/watervakken.shp'
w_lines_shp = '../data/input/waterlijnen.shp'
dem_tif = '../data/input/ahn3_05m_dtm.tif'

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

get_dem = False

# %%


def set_crs(gdf, crs=sources['default_crs']):
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


def explode(gdf):
    outgdf = gdf[gdf.geometry.type == "Polygon"]
    for idx, row in gdf[gdf.geometry.type == "MultiPolygon"].iterrows():
        multdf = gpd.GeoDataFrame(columns=gdf.columns)
        recs = len(row.geometry)
        multdf = multdf.append([row] * recs, ignore_index=True)
        for geom in range(recs):
            multdf.loc[geom, 'geometry'] = row.geometry[geom]
        outgdf = outgdf.append(multdf, ignore_index=True)
    return outgdf


def get_layer(src_layer):
    gdf_list = []

    for index, row in admins_gdf.iterrows():
        admin = row['name']
        admin_id = row['id']
        logging.info('{} {}'.format(src_layer, admin))

        src = sources[src_layer][admin]
        poly = row['geometry'].intersection(project_mask)

        if src['source'] in list(sources['services'].keys()):
            src_url = sources['services'][src['source']]

            if src_url['type'] == 'wfs':
                url = src_url['url']
                layer = src['layer']
                wfs = services.WFS(url)
                gdf = wfs.get_features(layer, poly)

            elif src_url['type'] == 'arcrest':
                url = src_url['url']
                layer = src['layer']
                output_format = 'geojson'

                if 'format' in list(src_url.keys()):
                    output_format = src_url['format']
                object_filter = ''
                if 'filter' in list(src.keys()):
                    object_filter = src['filter']
                rest = services.ArcREST(url, output_format=output_format)
                gdf = rest.get_features(
                    layer, poly, object_filter=object_filter)

            else:
                logging.error(
                    "type '{}' not supported".format(src_url['type']))

        elif src['source'] in list(sources['files'].keys()):
            src_path = sources['files'][src['source']]['path']
            gdf = gpd.read_file(src_path)
            gdf = set_crs(gdf)
            gdf = gdf[gdf.intersects(poly)]

        elif src['source'] == 'bgt':
            gdf = get_bgt.to_gdf(poly, layer="waterdeel",
                                 bronhouders=src['bronhouders'],
                                 end_registration='now',
                                 log_level="INFO")

        if 'attributes' in list(src.keys()):
            attributes = src['attributes']
            drop_cols = [col for col in list(gdf.columns) if not col in list(
                attributes.values()) + ['geometry']]
            gdf = gdf.drop(drop_cols, axis=1)
            gdf.rename(columns={value: key for key, value in attributes.items()},
                       inplace=True)
        else:
            drop_cols = [col for col in list(gdf.columns) if col != 'geometry']
            gdf = gdf.drop(drop_cols, axis=1)

        if 'translations' in list(src.keys()):
            translations = src['translations']
            for ident, trans in translations.items():
                for org, new in trans.items():
                    gdf[ident][gdf[ident] == org] = new

        if 'values' in list(src.keys()):
            values = src['values']
            for ident, val in values.items():
                gdf[ident] = val

        gdf.insert(0, 'src', src['source'])
        gdf.insert(0, 'admin_id', admin_id)
        gdf.insert(0, 'admin', admin)
        gdf_list += [gdf]

    # merge all into one dataframe
    gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True))
    gdf.crs = sources['default_crs']
    gdf.fillna(value=np.NaN, inplace=True)

    return gdf


# %% read admin_boundaries and update settings (later a class)
logging.info('Dissolving project mask')
project_gdf = gpd.read_file(project_shp)
project_gdf['dissolve'] = 0
project_mask = project_gdf.dissolve(by='dissolve').geometry[0]

# %% create an admins geodataframe, clip by project_mask
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

# %% get sources: water-level areas (peilvakken)
wl_areas_gdf = get_layer('water-level_areas')
wl_areas_gdf.to_file(wl_areas_shp)

# %% get sources: water areas (watervakken)
w_areas_gdf = get_layer('water_areas')
w_areas_gdf = w_areas_gdf.explode()
w_areas_gdf.to_file(w_areas_shp)

# %% get sources: water areas (watervakken)
w_lines_gdf = get_layer('water_lines')
w_lines_gdf = w_lines_gdf.explode()
#w_lines_gdf = to_lineString(w_lines_gdf)
w_lines_gdf.to_file(w_lines_shp)

# %% get sources: water areas (ahn)
if get_dem:
    get_ahn.to_tif(project_mask, dem_tif)

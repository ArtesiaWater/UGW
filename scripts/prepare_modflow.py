# public packages
import json
import logging
import os
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import rasterio
from rasterstats import zonal_stats

'''
ToDo:
    - uitzoeken waarom sommige src_wl velden leeg blijven terwijl er wel de juiste src_id_wl staat
'''

warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)

# %%
admins = json.loads(open('../config/administrations.json', 'r').read())
sources = json.loads(open('../config/sources.json', 'r').read())
project_shp = '../data/input/project.shp'

admins_shp = '../data/input/waterschappen.shp'
wl_areas_shp = '../data/input/peilvakken.shp'
w_areas_shp = '../data/input/watervakken.shp'
w_lines_shp = '../data/input/waterlijnen.shp'
dem_tif = '../data/input/ahn3_05m_dtm.tif'

mod_areas_shp = '../data/modflow/waterareas.shp'
mask_shp = '../data/modflow/mask.shp'
w_lines_clipped_shp = '../data/modflow/waterlijnen_masked.shp'

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

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


def sjoin_lines(x):
    try:
        gdf = gpd.clip(w_lines_gdf, x['geometry'])
        gdf['length'] = gdf['geometry'].length
        src_id = gdf.loc[gdf['length'].idxmax()]['src_id']
        return src_id
    except:
        return np.NaN


# %% read admin_boundaries and update settings (later a class)
logging.info('reading input-data')
admins_gdf = gpd.read_file(admins_shp)
wl_areas_gdf = gpd.read_file(wl_areas_shp)
w_areas_gdf = gpd.read_file(w_areas_shp)
w_lines_gdf = gpd.read_file(w_lines_shp)
project_gdf = gpd.read_file(project_shp)

logging.info('Dissolving project mask')
project_gdf['dissolve'] = 0
project_mask = project_gdf.dissolve(by='dissolve').geometry[0]

verd_gdf = gpd.GeoDataFrame(columns=['admin',
                                     'admin_id',
                                     'Z',
                                     'ZP',
                                     'WP',
                                     'BL',
                                     'BB',
                                     'CAT',
                                     'src_wa',
                                     'src_id_wa',
                                     'src_wla',
                                     'src_id_wla',
                                     'src_wl',
                                     'src_id_wl',
                                     'geometry'])

# %% compute 4 ModFLOW: compute ZP/WP
logging.info('merge ZP/WP')
wl_areas_gdf.fillna(value=np.NaN, inplace=True)
wl_areas_gdf['ZP'] = wl_areas_gdf.apply(lambda x:
                                        x['ZP'] if x['ZP'] == x['ZP'] else
                                        (x['ZFPB'] + x['ZFPO']) / 2 if x['ZFPB'] == x['ZFPB'] and x['ZFPO'] == x['ZFPO'] else
                                        x['VP'] if x['VP'] == x['VP'] else
                                        (x['FPB'] + x['FPO']) / 2 if x['FPB'] == x['FPB'] and x['FPO'] == x['FPO'] else
                                        np.NaN,
                                        axis=1)

wl_areas_gdf['WP'] = wl_areas_gdf.apply(lambda x:
                                        x['WP'] if x['WP'] == x['WP'] else
                                        (x['WFPB'] + x['WFPO']) / 2 if x['WFPB'] == x['WFPB'] and x['WFPO'] == x['WFPO'] else
                                        x['VP'] if x['VP'] == x['VP'] else
                                        (x['FPB'] + x['FPO']) / 2 if x['FPB'] == x['FPB'] and x['FPO'] == x['FPO'] else
                                        np.NaN,
                                        axis=1)

drop_cols = [col for col in wl_areas_gdf.columns if not col in ['src',
                                                                'src_id',
                                                                'ZP',
                                                                'WP',
                                                                'geometry']]
wl_areas_gdf = wl_areas_gdf.drop(drop_cols, axis=1)

# %% compute 4 ModFLOW: water-level areas
logging.info('compute median surface level (Z)')
with rasterio.open(dem_tif) as src:
    scale = src.scales[0]
    raster_profile = src.profile
    raster_data = src.read(1)
    affine = raster_profile['transform']
    raster_stats = zonal_stats(wl_areas_gdf,
                               raster_data,
                               affine=affine,
                               nodata=raster_profile['nodata'],
                               stats='percentile_50', raster_out=True)

wl_areas_gdf.insert(2, 'Z', [stat['percentile_50']
                             * scale for stat in raster_stats])

# %% compute 4 ModFLOW: overlay water areas with water-level areas
logging.info('overlay water-level areas with water areas')
mod_areas_gdf = gpd.overlay(w_areas_gdf, wl_areas_gdf, how='union')
mod_areas_gdf = mod_areas_gdf[mod_areas_gdf['src_1'] == mod_areas_gdf['src_1']]
mod_areas_gdf = mod_areas_gdf[mod_areas_gdf['Z'] == mod_areas_gdf['Z']]

mod_areas_gdf.rename(columns={'src_1': 'src_wa',
                              'src_id_1': 'src_id_wa',
                              'src_2': 'src_wla',
                              'src_id_2': 'src_id_wla', },
                     inplace=True)

# %% compute 4 ModFLOW: sjoin water lines
logging.info('spatial-join water lines with water areas')
# add the src_id of water-lines 2 mod_areas_gdf
mod_areas_gdf['src_id_wl'] = mod_areas_gdf.apply(
    (lambda x: sjoin_lines(x)), axis=1)

# join the attributes of w_lines_gdf 2 mod_areas_gdf on src_id_wl
w_lines_df = pd.DataFrame(w_lines_gdf.drop(['geometry'], axis=1))
w_lines_df.rename(columns={'src_id': 'src_id_wl',
                           'src': 'src_wl'}, inplace=True)

# as we exploded multi-polygons, we need to remove duplicates
w_lines_df.drop_duplicates(subset='src_id_wl',
                           keep=False, inplace=True)

mod_areas_gdf = mod_areas_gdf.merge(w_lines_df,
                                    on='src_id_wl',
                                    how='left',
                                    suffixes=('_wa', '_wl'))

# some laundering on attribute-names
mod_areas_gdf.rename(columns={'admin_wa': 'admin',
                              'admin_id_wa': 'admin_id'},
                     inplace=True)

mod_areas_gdf.fillna(value=np.NaN, inplace=True)

# %% compute 4 ModFLOW: fill ZP/WP if not existing
logging.info('fill ZP/WP if not exiting (ZP = Z - 0.5, WP = Z - 0.7')
for level, correction in {'ZP': 0.5, 'WP': 0.7}.items():
    mod_areas_gdf['ZP'] = mod_areas_gdf.apply(lambda x: x[level] if x[level] == x[level] else
                                              x['Z'] - correction, axis=1)


# %% compute 4 ModFLOW: merge WD
logging.info('merge water depths (areas first)')
for attrib in ['WD', 'WD_wa', 'WD_wl']:
    if not attrib in list(mod_areas_gdf.columns):
        mod_areas_gdf[attrib] = np.NaN

mod_areas_gdf['WD'] = mod_areas_gdf.apply((lambda x:
                                           x['WD'] if x['WD'] == x['WD'] else
                                           x['WD_wa'] if x['WD_wa'] == x['WD_wa'] else
                                           x['WD_wl'] if x['WD_wl'] == x['WD_wl'] else
                                           np.NaN), axis=1)

mod_areas_gdf = mod_areas_gdf.drop(['WD_wa', 'WD_wl'], axis=1)

# %% compute 4 ModFLOW: merge BL
logging.info('merge water bottom levels (BL) (areas first)')
for attrib in ['BL', 'BL_wa', 'BL_wl']:
    if not attrib in list(mod_areas_gdf.columns):
        mod_areas_gdf[attrib] = np.NaN

mod_areas_gdf['BL'] = mod_areas_gdf.apply((lambda x:
                                           x['BL'] if x['BL'] == x['BL'] else
                                           x['BL_wa'] if x['BL_wa'] == x['BL_wa'] else
                                           x['BL_wl'] if x['BL_wl'] == x['BL_wl'] else
                                           np.NaN), axis=1)

# %% compute 4 ModFLOW: fill BL by WP & WD default = WP - 0.3
logging.info('fill bottom level if not exiting (WP - 0.3)')
mod_areas_gdf = mod_areas_gdf.drop(['BL_wa', 'BL_wl'], axis=1)

mod_areas_gdf['BL'] = mod_areas_gdf.apply(lambda x:
                                          x['BL'] if x['BL'] == x['BL'] else
                                          x['WP'] - x['WD'] if x['WD'] == x['WD'] else
                                          x['WP'] - 0.3,
                                          axis=1)

# %% compute 4 ModFLOW: merge CAT (default = 2)
logging.info('merge water depths (areas first)')
for attrib in ['CAT', 'CAT_wa', 'CAT_wl']:
    if not attrib in list(mod_areas_gdf.columns):
        mod_areas_gdf[attrib] = np.NaN

mod_areas_gdf['CAT'] = mod_areas_gdf.apply((lambda x:
                                            x['CAT'] if x['CAT'] == x['CAT'] else
                                            x['CAT_wa'] if x['CAT_wa'] == x['CAT_wa'] else
                                            x['CAT_wl'] if x['CAT_wl'] == x['CAT_wl'] else
                                            3), axis=1)

mod_areas_gdf = mod_areas_gdf.drop(['CAT_wa', 'CAT_wl'], axis=1)

# %% compute 4 ModFLOW: merge BB
logging.info('merge water bottom widths (BB) (areas first)')
for attrib in ['BB', 'BB_wa', 'BB_wl']:
    if not attrib in list(mod_areas_gdf.columns):
        mod_areas_gdf[attrib] = np.NaN

mod_areas_gdf['BB'] = mod_areas_gdf.apply((lambda x:
                                           x['BB'] if x['BB'] == x['BB'] else
                                           x['BB_wa'] if x['BB_wa'] == x['BB_wa'] else
                                           x['BB_wl'] if x['BB_wl'] == x['BB_wl'] else
                                           np.NaN), axis=1)


# %% compute 4 ModFLOW: cleaning and writing

drop_cols = [
    col for col in mod_areas_gdf.columns if not col in list(verd_gdf.columns)]

mod_areas_gdf = mod_areas_gdf.drop(drop_cols, axis=1)

mod_areas_gdf = mod_areas_gdf[list(verd_gdf.columns)]

mod_areas_gdf.to_file(mod_areas_shp)

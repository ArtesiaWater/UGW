# -*- coding: utf-8 -*-

# public packages
import configparser
import json
import logging
import os
import shutil
import time
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import rasterio
from rasterstats import zonal_stats
from tqdm import tqdm

warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)
tqdm.pandas()

start = time.time()


# %% alle paden van output bestanden
config_ini = configparser.ConfigParser()
config_ini.read(r'config.ini')
project = config_ini['general']['project']

extent_shp = 'extent.shp'
if 'extent' in config_ini['general'].keys():
    extent_shp = config_ini['general']['extent']

# %% paths
os.chdir('../../config')
admins = {'file_name': Path(r'../config/administrations.json')}
sources = {'file_name': Path(r'../config/sources.json')}
validation = {'file_name': Path(r'../config/validation.json')}

admins.update(json.loads(open(admins['file_name'], 'r').read()))
sources.update(json.loads(open(sources['file_name'], 'r').read()))
validation.update(json.loads(open(validation['file_name'], 'r').read()))

data_dir = Path(f'../data/{project}').absolute().resolve()
input_dir = data_dir.joinpath('input')
project_shp = input_dir.joinpath(extent_shp)

admins_shp = input_dir.joinpath('waterschappen.shp')
wl_areas_shp = input_dir.joinpath('water-level_areas.shp')
w_areas_shp = input_dir.joinpath('water_areas.shp')
w_lines_shp = input_dir.joinpath('water_lines.shp')
dem_tif = input_dir.joinpath('ahn3_5m_dtm.tif')
bathymetry_shp = r'../data/sources/Bathymetry/bathymetry_masks.shp'

modflow_dir = Path(f'../data/{project}/modflow').absolute().resolve()
mod_areas_shp = modflow_dir.joinpath('waterareas.shp')
ma_verdict_shp = modflow_dir.joinpath('waterareas_verdict.shp')
ma_verdict_json = modflow_dir.joinpath('waterareas_verdict.geojson')
wla_verdict_shp = modflow_dir.joinpath('water-level_areas_verdict.shp')
wla_verdict_json = modflow_dir.joinpath('water-level_areas_verdict.geojson')
w_lines_bl_shp = modflow_dir.joinpath('waterlines.shp')


if modflow_dir.exists():
    shutil.rmtree(modflow_dir)
modflow_dir.mkdir(parents=True)


logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

progress = 0
percentage = -1

for file in [ma_verdict_json, wla_verdict_json]:
    if file.exists():
        file.unlink()

# %% functions


def set_crs(gdf, crs=sources['default_crs']):
    if gdf.crs == None:
        gdf.crs = crs
    else:
        if hasattr(pyproj, 'CRS'):
            update_crs = not pyproj.CRS(gdf.crs).equals(pyproj.CRS(crs))
        else:
            update_crs = pyproj.Proj(gdf.crs).srs != pyproj.Proj(init=crs).srs
        if update_crs:
            gdf = gdf.to_crs(crs)

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


def merge_attrib(x, attribs, default, transformations=None, validation=None):
    cols = list(x.keys())
    value = default
    for attrib in attribs:
        if attrib in cols:
            if x[attrib] == x[attrib]:
                value = x[attrib]

                if not transformations == None:  # apply transformations
                    if x['admin'] in list(transformations.keys()):
                        if str(value) in list(transformations[x['admin']].keys()):
                            value = transformations[x['admin']][str(value)]

    return value


def sjoin(row):

    lines_idx = list(w_lines_sindex.intersection(row['geometry'].bounds))
    lines_gdf = w_lines_gdf.loc[lines_idx]
    lines_gdf = lines_gdf[lines_gdf.intersects(row['geometry'])]

    if len(lines_gdf) > 0:
        lines_gdf = lines_gdf[lines_gdf['admin'] == row['admin']]

    if len(lines_gdf) == 0:
        src_id_wl = np.NaN
        src_wl = np.NaN
        idx_wl = np.NaN

    elif len(lines_gdf) == 1:
        src_id_wl = 'NL.{:.0f}.{}'.format(
            row['admin_id'], lines_gdf['src_id'].values[0])
        src_wl = lines_gdf['src'].values[0]
        idx_wl = lines_gdf['src'].index[0]

    else:
        lines_gdf = gpd.clip(lines_gdf, row['geometry'])
        lines_gdf['length'] = lines_gdf['geometry'].length
        src_id_wl = 'NL.{:.0f}.{}'.format(row['admin_id'],
                                          lines_gdf.loc[lines_gdf['length'].idxmax()]['src_id'])
        src_wl = lines_gdf.loc[lines_gdf['length'].idxmax()]['src']
        idx_wl = lines_gdf.loc[lines_gdf['length'].idxmax()].name

    # add src_id_wla

    areas_idx = list(wl_areas_sindex.intersection(row['geometry'].bounds))
    areas_gdf = wl_areas_gdf.loc[areas_idx]
    areas_gdf = areas_gdf[areas_gdf.intersects(row['geometry'])]

    if len(areas_gdf) > 0:
        areas_gdf = areas_gdf[areas_gdf['admin'] == row['admin']]

    if len(areas_gdf) == 0:
        src_id_wla = np.NaN
        src_wla = np.NaN
        idx_wla = np.NaN

    elif len(areas_gdf) == 1:
        src_id_wla = 'NL.{:.0f}.{}'.format(
            row['admin_id'], areas_gdf['src_id'].values[0])
        src_wla = areas_gdf['src'].values[0]
        idx_wla = areas_gdf['src'].index[0]

    else:
        areas_gdf = gpd.clip(wl_areas_gdf, row['geometry'])
        areas_gdf['area'] = areas_gdf['geometry'].area
        src_id_wla = 'NL.{:.0f}.{}'.format(row['admin_id'],
                                           areas_gdf.loc[areas_gdf['area'].idxmax()]['src_id'])
        src_wla = areas_gdf.loc[areas_gdf['area'].idxmax()]['src']
        idx_wla = areas_gdf.loc[areas_gdf['area'].idxmax()].name

    if bathymetry_gdf[bathymetry_gdf.intersects(row['geometry'])].empty:
        has_bath = 0
    else:
        has_bath = 1

    if src_id_wl in w_lines_bl_gdf['src_id_wl'].values:
        has_slope = 1
    else:
        has_slope = 0

    return src_wl, src_id_wl, idx_wl, src_wla, src_id_wla, idx_wla, has_bath, has_slope


# %% read admin_boundaries and update settings (later a class)
logging.info('reading input-data')
admins_gdf = gpd.read_file(admins_shp)
wl_areas_gdf = gpd.read_file(wl_areas_shp)
w_areas_gdf = gpd.read_file(w_areas_shp)
w_lines_gdf = gpd.read_file(w_lines_shp)
project_gdf = gpd.read_file(project_shp)
bathymetry_gdf = gpd.read_file(bathymetry_shp)

logging.info('Dissolving project mask')
project_gdf['dissolve'] = 0
project_mask = project_gdf.dissolve(by='dissolve').geometry[0]

# %% build indeces
logging.info('indexing water-lines and water-level-areas')
istart = time.time()
w_lines_sindex = w_lines_gdf.sindex
wl_areas_sindex = wl_areas_gdf.sindex
logging.info('indeces build in {} seconds'.format(
    round(time.time() - istart, 1)))

# %%
logging.info('exploding water-areas')
estart = time.time()
w_areas_gdf = explode(w_areas_gdf)
#w_areas_gdf = explode(gpd.clip(w_areas_gdf,project_mask))
logging.info('exploded in {} seconds'.format(round(time.time() - estart, 1)))

# %% WL AREAS: compute Z
logging.info('compute median surface level on water-level areas (Z)')
with rasterio.open(dem_tif) as src:
    scale = src.scales[0]
    raster_profile = src.profile
    raster_data = src.read(1)
    affine = raster_profile['transform']
    raster_stats = zonal_stats(wl_areas_gdf,
                               raster_data,
                               affine=affine,
                               nodata=raster_profile['nodata'],
                               stats=['percentile_50', 'percentile_90'], raster_out=False)

Z = np.array([stat['percentile_50'] for stat in raster_stats])
Z[Z != None] *= scale
wl_areas_gdf.insert(2, 'Z50', Z)

Z = np.array([stat['percentile_90'] for stat in raster_stats])
Z[Z != None] *= scale
wl_areas_gdf.insert(2, 'Z90', Z)


# %% WL AREAS: validate and compute ZP/WP
wl_areas_gdf.fillna(value=np.NaN, inplace=True)


def combine_wl(gdf):

    gdf['SUB_wla'] = 0
    clean_cols = [col for col in gdf.columns if not col in [
        'admin', 'admin_id', 'geometry', 'src', 'src_id']]
    verdict_gdf = gdf.drop(clean_cols, axis=1)
    for col in clean_cols:
        verdict_gdf[col] = ''

    verdict_gdf.insert(0, 'mark', 10)

    exclude = ['RWS']

    for index, row in gdf.iterrows():
        #        has_zp = False
        #        has_wp = False

        params = ['VP', 'FPB', 'FPO', 'ZP',
                  'ZFPB', 'ZFPO', 'WP', 'WFPB', 'WFPO']
        for param in params:
            if not param in list(row.keys()):
                row[param] = np.NaN

        if not row['admin'] in exclude:
            row.fillna(value=np.NaN, inplace=True)
            mark = 10

            # compare all levels with Z
            level_cols = [col for col in gdf.columns if col.isupper() and not col in [
                'Z50', 'Z90', 'SUB_wla']]
            for col in level_cols:
                if row[col] == row[col]:
                    if not row[col] <= row['Z50']:
                        if row[col] <= row['Z90']:
                            verdict_gdf.loc[index, col] = '{}doubtfull:{}>Z50 (doubtfull);'.format(verdict_gdf.loc[index, col],
                                                                                                   col)
                            mark -= 1
                        else:
                            verdict_gdf.loc[index, col] = '{}doubtfull:{}>Z90 (deleted);'.format(verdict_gdf.loc[index, col],
                                                                                                 col)
                            row[col] = np.NaN
                            mark -= 2

                if row[col] == row[col]:
                    if not row[col] > row['Z50'] - 10:
                        verdict_gdf.loc[index, col] = '{}error:{}<<<Z (deleted);'.format(verdict_gdf.loc[index, col],
                                                                                         col)
                        row[col] = np.NaN
                        mark -= 1

            # make vp column
            if not row['VP'] == row['VP']:
                if (row['FPB'] == row['FPB']) and (row['FPO'] == row['FPO']):
                    row['VP'] = (row['FPB'] + row['FPO']) / 2

            # make zp column
            if not row['ZP'] == row['ZP']:
                if (row['ZFPB'] == row['ZFPB']) and (row['ZFPO'] == row['ZFPO']):
                    row['ZP'] = (row['ZFPB'] + row['ZFPO']) / 2
#                    has_zp = True
                elif row['VP'] == row['VP']:
                    row['ZP'] = row['VP']
                else:
                    verdict_gdf.loc[index, 'ZP'] = '{}error:no ZP (estimated Z - 0.5);'.format(
                        verdict_gdf.loc[index, 'ZP'])
                    row['ZP'] = row['Z50'] - 0.5
                    mark = 1
#            else:
#                has_zp = True

            # make wp column
            if not row['WP'] == row['WP']:
                if (row['WFPB'] == row['WFPB']) and (row['WFPO'] == row['WFPO']):
                    row['WP'] = (row['WFPB'] + row['WFPO']) / 2
#                    has_wp = True
                elif row['VP'] == row['VP']:
                    row['WP'] = row['VP']
                else:
                    verdict_gdf.loc[index, 'WP'] = '{}error:no WP (estimated Z - 0.7);'.format(
                        verdict_gdf.loc[index, 'WP'])
                    row['WP'] = row['Z50'] - 0.7
                    mark = 1
#            else:
#                has_wp = True

             # check if zp>wp
            if not row['WP'] <= row['ZP']:
                verdict_gdf.loc[index, 'WP'] = '{}doubtfull:WP>ZP;'.format(
                    verdict_gdf.loc[index, 'WP'])
                verdict_gdf.loc[index, 'ZP'] = '{}doubtfull:ZP<WP;'.format(
                    verdict_gdf.loc[index, 'ZP'])
                mark -= 1

            if not row['WP'] == row['ZP']:
                row['SUB_wla'] = 1
            verdict_gdf.loc[index, 'mark'] = max(mark, 1)
            gdf.loc[index] = row

    return gdf, verdict_gdf


# %%
logging.info('validate and combine ZP/WP in water-level areas')
wl_areas_gdf, wla_verdict_gdf = combine_wl(wl_areas_gdf)

wla_verdict_gdf.to_file(wla_verdict_shp)

# %% WL AREAS: drop columns
drop_cols = [col for col in wl_areas_gdf.columns if not col in ['src',
                                                                'src_id',
                                                                'Z50',
                                                                'ZP',
                                                                'WP',
                                                                'SUB_wla',
                                                                'geometry',
                                                                'admin',
                                                                'admin_id']]
wl_areas_gdf = wl_areas_gdf.drop(drop_cols, axis=1)
wl_areas_gdf.rename(columns={'Z50': 'Z'},
                    inplace=True)

# %% CREATE slope gdf
if ('BLU' in w_lines_gdf.columns) and ('BLD' in w_lines_gdf.columns):
    logging.info('createing slope water-lines')
    w_lines_gdf.fillna(value=np.NaN, inplace=True)
    w_lines_bl_gdf = w_lines_gdf[w_lines_gdf['BLU'].notnull()]
    w_lines_bl_gdf = w_lines_bl_gdf[w_lines_bl_gdf['BLD'].notnull()]
    drop_cols = [col for col in w_lines_bl_gdf.columns
                 if col.isupper() and not col in ['BLU', 'BLD', 'WD']]
    w_lines_bl_gdf = w_lines_bl_gdf.drop(drop_cols, axis=1)
    w_lines_bl_gdf.rename(columns={'src': 'src_wl',
                                   'src_id': 'src_id_wl'},
                          inplace=True)

    w_lines_bl_gdf['src_id_wl'] = w_lines_bl_gdf.apply((lambda x: 'NL.{:.0f}.{}'.format(x['admin_id'],
                                                                                        x['src_id_wl'])), axis=1)

    w_lines_bl_gdf.to_file(w_lines_bl_shp)
else:
    w_lines_bl_gdf = gpd.GeoDataFrame(columns=['src_id_wl'])

# %% JOIN: spatial join water-level areas and water lines with water areas
logging.info('spatial join water areas with water-level areas and water lines')
jstart = time.time()

mod_areas_gdf = w_areas_gdf

mod_areas_gdf[['src_wl',
               'src_id_wl',
               'idx_wl',
               'src_wla',
               'src_id_wla',
               'idx_wla',
               'has_bath',
               'has_slope']] = mod_areas_gdf.progress_apply(sjoin,
                                                            axis=1,
                                                            result_type="expand")

# %%
logging.info('join in {} seconds'.format(round(time.time() - jstart, 1)))
mod_areas_gdf.rename(columns={'src': 'src_wa',
                              'src_id': 'src_id_wa'},
                     inplace=True)
# %%
# join the attributes of w_lines_gdf 2 mod_areas_gdf on src_id_wl
w_lines_df = pd.DataFrame(w_lines_gdf.drop(['geometry'], axis=1))
w_lines_df = pd.DataFrame(w_lines_df.drop(
    ['admin', 'admin_id', 'src_id', 'src'], axis=1))


# %%
mod_areas_gdf = mod_areas_gdf.merge(w_lines_df,
                                    left_on='idx_wl',
                                    right_index=True,
                                    how='left',
                                    suffixes=('_wa', '_wl'))

# %%
wl_areas_df = pd.DataFrame(wl_areas_gdf.drop(['geometry'], axis=1))
wl_areas_df = pd.DataFrame(wl_areas_df.drop(
    ['admin', 'admin_id', 'src_id', 'src'], axis=1))

mod_areas_gdf = mod_areas_gdf.merge(wl_areas_df,
                                    left_on='idx_wla',
                                    right_index=True,
                                    how='left',
                                    suffixes=('_wa', '_wla'))

mod_areas_gdf.fillna(value=np.NaN, inplace=True)


# %% MERGE: names
mstart = time.time()
logging.info('merge names')
attrib = 'name'
attribs = ['{}_wl'.format(attrib), '{}_wa'.format(attrib), attrib]
mod_areas_gdf[attrib] = mod_areas_gdf.apply(
    (lambda x: merge_attrib(x, attribs, np.NaN)), axis=1)

# %% MERGE: water depths
logging.info('merge water depths (areas first)')
attrib = 'WD'
attribs = ['{}_wl'.format(attrib), '{}_wa'.format(attrib), attrib]
mod_areas_gdf[attrib] = mod_areas_gdf.apply(
    (lambda x: merge_attrib(x, attribs, np.NaN)), axis=1)

# %% MERGE: bottom level
logging.info('merge water bottom levels (BL) (areas first)')
attrib = 'BL'
attribs = ['{}_wl'.format(attrib), '{}_wa'.format(attrib), attrib]
mod_areas_gdf[attrib] = mod_areas_gdf.apply(
    (lambda x: merge_attrib(x, attribs, np.NaN)), axis=1)

# %% MERGE: supply.
logging.info('merge supply')
attrib = 'SUB'
default = validation[attrib]['default']
transformations = validation[attrib]['transformations']
attribs = ['{}_wla'.format(attrib), '{}_wa'.format(attrib), attrib]
mod_areas_gdf[attrib] = mod_areas_gdf.apply((lambda x: merge_attrib(x,
                                                                    attribs,
                                                                    default,
                                                                    transformations=transformations)), axis=1)

# %% MERGE: merge BB
logging.info('merge bottom widths (areas first)')
attrib = 'BB'
attribs = ['{}_wl'.format(attrib), '{}_wa'.format(attrib), attrib]
mod_areas_gdf[attrib] = mod_areas_gdf.apply(
    (lambda x: merge_attrib(x, attribs, np.NaN)), axis=1)

# %% MERGE: merge CAT (default = 2)
logging.info('merge category (areas first)')
attrib = 'CAT'
default = validation[attrib]['default']
transformations = validation[attrib]['transformations']
attribs = ['{}_wl'.format(attrib), '{}_wa'.format(attrib), attrib]
mod_areas_gdf[attrib] = mod_areas_gdf.apply((lambda x: merge_attrib(x,
                                                                    attribs,
                                                                    default,
                                                                    transformations=transformations)), axis=1)

logging.info('merging in {} seconds'.format(round(time.time() - mstart, 1)))
# %% VALIDATE:


def validate_wa(row):
    exclude = ['RWS']
    mark = 10
    CAT_vdct = BL_vdct = ZP_vdct = WP_vdct = BB_vdct = ''

    if not row['admin'] in exclude:
        if not row['CAT'] == row['CAT']:
            CAT_vdct = '{}not specified (estimated 2);'.format(CAT_vdct)
            row['CAT'] = 2
            mark -= 1

        not_specified = ' error: not specified'
        row['BL'] = float(row['BL'])
        if row['BL'] == row['BL']:
            not_specified = ''
            # check if BL < WP:
            if row['WP'] == row['WP']:
                if not row['BL'] <= row['WP']:
                    BL_vdct = '{}error: BL > WP (deleted);'.format(CAT_vdct)
                    row['BL'] = np.NaN
                    mark -= 5
            # check if BL < Z:
            if row['Z'] == row['Z']:
                if not row['BL'] <= row['Z']:
                    BL_vdct = '{}error: BL > Z (deleted);'.format(BL_vdct)
                    row['BL'] = np.NaN
                    mark -= 5

        if not row['BL'] == row['BL']:
            if (row['WD'] == row['WD']) & (row['WP'] == row['WP']):
                BL_vdct = '{}{} (estimated WP - WD);'.format(BL_vdct,
                                                             not_specified)
                row['BL'] = row['WP'] - row['WD']
                mark -= 2
            else:
                if row['WP'] == row['WP']:
                    if row['CAT'] == 1:
                        subtract = 1.5
                        mark -= 5
                    else:
                        subtract = 0.3
                        mark -= 4
                    BL_vdct = '{}{} (estimated WP - {});'.format(BL_vdct,
                                                                 not_specified,
                                                                 subtract)
                    row['BL'] = row['WP'] - subtract

                else:
                    BL_vdct = '{}{} cannot be determined;'.format(BL_vdct,
                                                                  not_specified)
                    mark = 1

        if not row['ZP'] == row['ZP']:
            ZP_vdct = '{}warning: not specified;'.format(ZP_vdct)

        if not row['WP'] == row['WP']:
            WP_vdct = '{}warning: not specified;'.format(WP_vdct)

    else:
        row['CAT'] = 1

    row['BB'] = float(row['BB'])
    if row['BB'] == row['BB']:
        if row['BB'] > row['geometry'].exterior.length / np.pi:
            BB_vdct = '{}error: BB unrealistic large (deleted);'.format(
                BB_vdct)
            row['BB'] = np.NaN
            mark -= 1
        elif row['BB'] < 0.5:
            BB_vdct = '{}error: BB < 0.5 (deleted);'.format(BB_vdct)
            row['BB'] = np.NaN
            mark -= 1

    return row['CAT'], np.round(row['BL'], 2), np.round(row['BB'], 2), max(mark, 1), CAT_vdct, BL_vdct, ZP_vdct, WP_vdct, BB_vdct


logging.info('validating water-areas')

vstart = time.time()
#mod_areas_gdf, ma_verdict_gdf = combine_wa(mod_areas_gdf)
mod_areas_gdf[['CAT',
               'BL',
               'BB',
               'mark',
               'CAT_vdct',
               'BL_vdct',
               'ZP_vdct',
               'WP_vdct',
               'BB_vdct']] = mod_areas_gdf.apply(validate_wa,
                                                 axis=1,
                                                 result_type="expand")
ma_verdict_gdf = mod_areas_gdf
drop_cols = [col for col in mod_areas_gdf.columns if col.isupper()]
ma_verdict_gdf = ma_verdict_gdf.drop(drop_cols, axis=1)

verdict_parameters = [col.replace('_vdct', '')
                      for col in ma_verdict_gdf.columns if '_vdct' in col]
ma_verdict_gdf.rename(columns={col: col.replace('_vdct', '')
                               for col in ma_verdict_gdf.columns
                               if '_vdct' in col},
                      inplace=True)

logging.info('validated in {} seconds'.format(round(time.time() - vstart, 1)))


# %% CLEAN and WRITE result
wstart = time.time()
logging.info('writing results-files')

schema = {'geometry': 'Polygon',
          'properties': {'admin': 'str',
                         'admin_id': 'int',
                         'name': 'str',
                         'Z': 'float: 5.2',
                         'ZP': 'float: 5.2',
                         'WP': 'float: 5.2',
                         'WD': 'float: 5.2',
                         'BL': 'float: 5.2',
                         'BB': 'float: 5.2',
                         'CAT': 'int',
                         'SUB': 'int',
                         'has_bath': 'int',
                         'has_slope': 'int',
                         'src_wa': 'str',
                         'src_id_wa': 'str',
                         'src_wla': 'str',
                         'src_id_wla': 'str',
                         'src_wl': 'str',
                         'src_id_wl': 'str'}
          }

attributes = list(schema['properties'].keys()) + ['geometry']
drop_cols = [col for col in mod_areas_gdf.columns if not col in attributes]
mod_areas_gdf = mod_areas_gdf.drop(drop_cols, axis=1)
mod_areas_gdf = mod_areas_gdf[attributes]
mod_areas_gdf.to_file(mod_areas_shp, schema=schema)

schema['properties'] = {**{key: value for key, value in schema['properties'].items()
                           if key.islower() and not 'src' in key},
                        'mark': 'int',
                        **{par: 'str' for par in verdict_parameters},
                        **{key: value for key, value in schema['properties'].items()
                            if 'src' in key}}

attributes = list(schema['properties'].keys()) + ['geometry']
drop_cols = [col for col in ma_verdict_gdf.columns if not col in attributes]
ma_verdict_gdf = ma_verdict_gdf.drop(drop_cols, axis=1)
ma_verdict_gdf = ma_verdict_gdf[attributes]
ma_verdict_gdf.to_file(ma_verdict_shp, schema=schema)

logging.info('written in {} seconds'.format(round(time.time() - wstart, 1)))
logging.info('done in {} seconds'.format(round(time.time() - start, 1)))

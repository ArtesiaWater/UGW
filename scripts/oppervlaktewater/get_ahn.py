# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:31:31 2019

@author: danie
"""

import time

import geopandas as gpd
import numpy as np
import rasterio
import requests
from owslib.wcs import WebCoverageService
from rasterio import Affine, features
from rasterio.enums import Resampling
from rasterio.io import MemoryFile
from rasterio.windows import Window
from shapely.geometry import Polygon

# %% script objects
max_size = 2000  # max pixels in x and y direction that can be loaded from NGR WCS
nodata = -32768  # nodata value to be used in int16
dtype = rasterio.int16  # dtype for raster
attempts = 10  # attempts it will try to succesfully download en process a data-array from NGR
timeout = 300

# %% functions


def get_layers(ahn_version='ahn3'):
    '''
    Returns a list of layers which can be used in function "to_tif"

    Parameters
    ----------    
    ahn_version: str, 
        specify 'ahn2' or 'ahn3'
    '''

    url = 'https://geodata.nationaalgeoregister.nl/{}/wcs?'.format(ahn_version)
    wcs = WebCoverageService(url)
    return list(wcs.contents.keys())


def to_tif(mask_poly, tif_file, layer='ahn3_05m_dtm', cell_size=0.5, src_nodata=None, overviews=[5, 25]):
    '''
    Download an ahn-layer clipped by a shapely-polygon

    Parameters
    ----------    
    mask_poly: shapely polygon geometry
        Polygon used as clipping mask
    tif_file: str, file object or pathlib.Path object
        Location of the tif-file to be stored
    layer: str, optional
        NGR WCS layer to be downloaded. By default 'ahn3_05m_dtm'
    cell_size: int,float
        Cell_size in which the layer will be downloaded and stored
    src_nodata: int, float, optional
        Over-write the nodata value returned by the WCS. Usefull @ ahn2, as 
        nodata is not provided in gtiff profile
    overviews: list, optional
        Specify a list of raster-overviews in m. E.g.: overviews=[5,25] and 
        cell_size=0.5 will create 2 overviews with a cell size of 5m and 25m. 
        With the same overviews and cell_size=5 only an overview of 25m will 
        be included
    '''

    url = 'https://geodata.nationaalgeoregister.nl/{}/wcs?'.format(
        layer[:layer.find('_')])
    wcs = WebCoverageService(url, version='1.0.0', timeout=timeout)

    # see if the mask_poly is within the boundary of the wcs layer
    bboxes = wcs.contents[layer].boundingboxes
    xmin, ymin, xmax, ymax = [item['bbox']
                              for item in bboxes if item['nativeSrs'] == 'EPSG:28992'][0]
    wcs_poly = Polygon([(xmin, ymin), (xmin, ymax),
                        (xmax, ymax), (xmax, ymin), (xmin, ymin)])
    if not wcs_poly.intersects(mask_poly):
        print('Bounds of poly-mask do not intersect with bounds of wcs-layer')

    else:
        bounds = list(mask_poly.bounds)
        # xmin, ymin rounddown 2 cellsize
        bounds[0], bounds[2] = [
            round(bounds[idx] / cell_size - cell_size) * cell_size for idx in [0, 2]]
        # xmax, ymax rounddown 2 cellsize
        bounds[1], bounds[3] = [
            round(bounds[idx] / cell_size + cell_size) * cell_size for idx in [1, 3]]

        profile = {'driver': 'GTiff',
                   'dtype': dtype,
                   'nodata': -32768,
                   'width': int((bounds[2] - bounds[0]) / cell_size),
                   'height': int((bounds[3] - bounds[1]) / cell_size),
                   'count': 1,
                   'crs': 'epsg:28992',
                   'BIGTIFF': "IF_SAFER",
                   'transform': Affine(cell_size, 0.0, bounds[0], 0.0, -cell_size, bounds[3]),
                   'tiled': True,
                   'interleave': 'band',
                   'compress': 'deflate',
                   'predictor': 2,
                   'blockxsize': 256,
                   'blockysize': 256}

        cols = int(np.ceil((bounds[2] - bounds[0]) / cell_size / max_size))
        rows = int(np.ceil((bounds[3] - bounds[1]) / cell_size / max_size))

        window_width = int((bounds[2] - bounds[0]) / cols / cell_size)
        window_height = int((bounds[3] - bounds[1]) / rows / cell_size)

        with rasterio.open(tif_file, 'w', **profile) as dst:
            dst.scales = [0.01]
            for row in range(rows):
                for col in range(cols):
                    xmin = bounds[0] + (col * window_width * cell_size)
                    ymax = bounds[3] - (row * window_height * cell_size)
                    xmax = xmin + (window_width * cell_size)
                    ymin = ymax - (window_height * cell_size)

                    bound_poly = Polygon(
                        [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)])
                    if bound_poly.intersects(mask_poly):
                        print(
                            'NGR download: (row: {}/{}, col: {}/{})'.format(row + 1, rows, col + 1, cols))

                        attempt = 1
                        succeed = False

                        while attempt <= attempts and not succeed:
                            try:
                                requestbbox = (xmin, ymin, xmax, ymax)
                                requestwidth = window_width
                                requestheight = window_height

                                gc = wcs.getCoverage(identifier=layer,
                                                     bbox=requestbbox,
                                                     format='GEOTIFF_FLOAT32', width=requestwidth,
                                                     height=requestheight,
                                                     crs='EPSG:28992')
                                with MemoryFile(gc) as memfile:
                                    with memfile.open() as src:
                                        data = src.read(1)
                                        if src_nodata == None:
                                            src_nodata = src.profile['nodata']
                                        data = np.where(
                                            data == src_nodata, nodata, (data * 100).astype(rasterio.int16))
                                        if not bound_poly.within(mask_poly):
                                            geometry = bound_poly.intersection(
                                                mask_poly)
                                            mask = rasterio.features.rasterize(
                                                [(geometry, 1)],
                                                out_shape=data.shape,
                                                transform=src.profile['transform'],
                                                fill=0,
                                                all_touched=True,
                                                dtype=dtype)
                                            data = np.where(
                                                mask == 1, data, nodata)
                                succeed = True
                            except Exception as e:
                                print(
                                    'FAILED ATTEMPT ({}/{}): {} RETRYING 5 SECS'.format(attempt, attempts, e))
                                attempt += 1
                                time.sleep(5)
                                pass

                        dst.write(data.astype(rasterio.int16), window=Window(col * window_width, row * window_height,
                                                                             window_width,
                                                                             window_height), indexes=1)

            if not overviews == None:
                print('creating overviews')
                factors = [int(size / cell_size)
                           for size in overviews if size > cell_size]
                dst.build_overviews(factors, Resampling.average)
                dst.update_tags(ns='rio_overview', resampling='average')


def io(file_in, tif_file, layer='ahn3_05m_dtm', cell_size=0.5, src_nodata=None, overviews=[5, 25]):
    '''downloads and saves BGT data for the extent of a specified polygon-file

    Parameters
    ----------    
    file_in: file with polygon geometries (e.g. ESRI-shapefile) on local drive
    tif_out: location for the tif-file to be written
    layer: ahn-layer to be downloades (default ahn3_05m_dtm)
    cell_size: cell_size of resulting tif (default 0.5m)
    '''

    # mask from file_in
    gdf = gpd.read_file(file_in)
    gdf['dissolve'] = 0
    mask_poly = gdf.dissolve(by='dissolve').geometry[0]

    to_tif(mask_poly, tif_file, layer=layer, cell_size=cell_size,
           src_nodata=src_nodata, overviews=overviews)

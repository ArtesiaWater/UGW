import datetime
import json
import logging
import os
import time

import fiona
import geopandas as gpd
import requests
from osgeo import gdal


def to_gdf(mask_poly, layer="waterdeel", bronhouders=None, end_registration='now', log_level="INFO"):

    logging.basicConfig(level=os.environ.get("LOGLEVEL", log_level))

    if end_registration == 'now':
        end_registration = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    api_url = 'https://download.pdok.io'
    # start download proces
    url = '{}/lv/bgt/api/v1/full/custom'.format(api_url)
    body = {"format": "gmllight",
            "featuretypes": [layer]}

    # else:
    body["geofilter"] = "POLYGON(({}))".format(
        ",".join(['{} {}'.format(coords[0], coords[1])
                  for coords in mask_poly.exterior.coords])
    )

    headers = {'content-type': 'application/json'}

    response = requests.post(url, headers=headers, data=json.dumps(body))

    # check api-status, if completed, download
    if response.status_code in range(200, 300):
        running = True
        url = '{}{}'.format(api_url, response.json()[
                            "_links"]["status"]["href"])

        while running:
            response = requests.get(url)
            if response.status_code in range(200, 300):
                logging.info('{}% ({})'.format(response.json()[
                             "progress"], response.json()['status']))
                status = response.json()['status']
                if status == "COMPLETED":
                    running = False
                else:
                    time.sleep(2)
            else:
                logging.error(response.text)
                running = False
    else:
        logging.error(response.text)

    logging.info('converting to gdf')
    response = requests.get('{}{}'.format(
        api_url, response.json()["_links"]["download"]["href"]))

    vsif = fiona.ogrext.buffer_to_virtual_file(bytes(response.content))
    vsiz = vsif + '.zip'
    gdal.Rename(vsif, vsiz)

    fc = fiona.Collection(vsiz, vsi='zip')
    gdf = gpd.GeoDataFrame.from_features(
        [feature for feature in fc], crs='epsg:28992')

    # select polygons after end_registration
    gdf = gdf[(gdf['eindRegistratie'] != gdf['eindRegistratie'])
              | (gdf['eindRegistratie'] > end_registration)]

    # select polygons of specific bronhouders
    if not bronhouders == None:
        gdf = gdf[gdf['bronhouder'].isin(bronhouders)]

    # select polygons within polygon mask
    gdf = gdf[gdf.intersects(mask_poly)]

    # re-order columns
    columns = [col for col in gdf.columns if not col ==
               'geometry'] + ['geometry']
    gdf = gdf[columns]

    return gdf

# %%


def io(file_in, file_out, layer="waterdeel", bronhouders=None, end_registration='now', log_level="INFO"):
    '''downloads and saves BGT data for the extent of a specified polygon-file

    Parameters
    ----------    
    file_in: file with polygon geometries (e.g. ESRI-shapefile) on local drive
    file_out: file path on local drive for writing the result
    layer: bgt-layer to be downloaded (default: waterdeel). See https://download.pdok.io/lv/bgt/viewer/
    bronhouders: list with codes of bronhouders. If None (default), all will be returned
    end_registration:   string with format %Y-%m-%dT%H:%M:%S or 'now' (default). 
                        Only polygons with eindRegistratie > end_registration will be returned
    log_level: logging level (see https://docs.python.org/3/library/logging.html)
    '''

    # mask from file_in
    gdf = gpd.read_file(file_in)
    gdf['dissolve'] = 0
    mask_poly = gdf.dissolve(by='dissolve').geometry[0]

    # load bgt in data-frame
    gdf = to_gdf(mask_poly,
                 layer=layer,
                 bronhouders=bronhouders,
                 end_registration=end_registration,
                 log_level=log_level)

    # save dataframe to file
    gdf.to_file(file_out)

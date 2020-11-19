import requests
import numpy as np
import pandas as pd
from pyproj import Proj, transform
import os
import zipfile
from owslib.wfs import WebFeatureService
from owslib.etree import etree
from owslib.fes import PropertyIsEqualTo
import json
import geopandas as gpd
from tqdm import tqdm
from shapely.geometry import Point

# https://www.rijkswaterstaat.nl/rws/opendata/
# https://www.rijkswaterstaat.nl/rws/opendata/DistributielaagWebservices-SUM-2v7.pdf
# https://rahuls01.github.io/?python#introduction


def get_catalogus(raw_response=False):
    """Get the catalogus from rws

    Parameters
    ----------
    raw_response : bool, optional
        return the raw json-response without any manipulations.

    Returns
    -------
    catalogus : dict
        A dictionary with the contents of the catalogus.

    """
    url = 'https://waterwebservices.rijkswaterstaat.nl/METADATASERVICES_DBO/OphalenCatalogus/'
    json_input = {"CatalogusFilter": {
        "Grootheden": True,
    }
    }
    catalogus = post_data(url, json_input)
    if raw_response:
        return catalogus
    for aquo_metadata in catalogus['AquoMetadataLijst']:
        flatten(aquo_metadata, 'Grootheid')
    catalogus['LocatieLijst'] = [convert_to_rd(
        x) for x in catalogus['LocatieLijst']]
    return catalogus


def get_locations(catalogus=None, extent=None, grootheid=None, omschrijving=None):
    """Get the measurement locations from the catalogus

    Parameters
    ----------
    catalogus : dict
        A dictionary with the contents of the catalogus.
    extent : list or tuple, optional
        When supplies, only return location within the extent [xmin, xmax, ymin, ymax)].
    grootheid : string, optional
        The code of the parameter that is required to be measured at the locations.
    omschrijving : string, optional
        The omschrijving of the parameter that is required to be measured at the locations.

    Returns
    -------
    catalogus : dict
        A dictionary with the contents of the catalogus

    """
    if catalogus is None:
        catalogus = get_catalogus()
    loc = pd.DataFrame(catalogus['LocatieLijst']
                       ).set_index('Locatie_MessageID')

    if grootheid is not None and omschrijving is not None:
        raise(Exception('Do not supply both grootheid and omschrijving of a parameter'))
    if grootheid is not None or omschrijving is not None:
        par = get_parameters(catalogus)
        if omschrijving is None:
            mask = par['Code'] == grootheid
        else:
            mask = par['Parameter_Wat_Omschrijving'] == omschrijving
        if not mask.any():
            if omschrijving is None:
                raise(Exception('Parameter {} not found'.format(grootheid)))
            else:
                raise(Exception('Parameter {} not found'.format(omschrijving)))
        par_index = par.index[mask][0]

        par_loc = pd.DataFrame(catalogus['AquoMetadataLocatieLijst']).set_index(
            'AquoMetaData_MessageID')
        loc_index = par_loc.loc[par_index]['Locatie_MessageID']
        loc = loc.loc[loc_index]
    if extent is not None:
        mask = (loc['x'] >= extent[0]) & (loc['x'] <= extent[1]) & \
            (loc['y'] >= extent[2]) & (loc['y'] <= extent[3])
        loc = loc[mask]
    return loc


def convert_to_rd(loc):
    inProj = Proj(init='epsg:' + loc['Coordinatenstelsel'])
    outProj = Proj(init='epsg:28992')
    loc['x'], loc['y'] = transform(inProj, outProj, loc['X'], loc['Y'])
    loc['Coordinatenstelsel'] = '28992'
    return loc


def get_parameters(catalogus=None):
    if catalogus is None:
        catalogus = get_catalogus()
    par = pd.DataFrame(catalogus['AquoMetadataLijst']).set_index(
        'AquoMetadata_MessageID')
    return par


def get_grootheid_code(omschrijving, parameters=None, catalogus=None):
    if parameters is None:
        parameters = get_parameters(catalogus)
    mask = parameters['Parameter_Wat_Omschrijving'] == omschrijving
    code = parameters.loc[mask, 'Code'].iloc[0]
    return code


def get_locatie(locations, name=None, code=None):
    if isinstance(locations, pd.Series) and name is None and code is None:
        return {'X': locations.X, 'Y': locations.Y, 'Code': locations.Code}
    if (name is None and code is None) or (name is not None and code is not None):
        raise(Exception('Either supply name or code of  a single location'))
    if name is not None:
        mask = locations['Naam'] == name
    else:
        mask = locations['Code'] == code
    if mask.sum() == 0:
        if name is not None:
            raise(Exception('Location with name {} not found'.format(name)))
        else:
            raise(Exception('Location with code {} not found'.format(name)))
    if mask.sum() > 1:
        raise(Exception('Multiple locations found'))
    locatie = locations.loc[mask, ['X', 'Y', 'Code']].squeeze().to_dict()
    return locatie


def get_locatie_lijst(locations):
    if isinstance(locations, pd.Series):
        return [get_locatie(locations)]
    return [get_locatie(row) for ind, row in locations.iterrows()]


def get_periode(tmin=None, tmax=None):
    if tmax is None:
        tmax = pd.to_datetime('now')
    else:
        tmax = pd.to_datetime(tmax)
    if tmin is None:
        tmin = tmax - pd.DateOffset(months=1)
    else:
        tmin = pd.to_datetime(tmin)
    periode = {"Begindatumtijd": tmin.strftime('%Y-%m-%dT%H:%M:%S.000+01:00'),
               "Einddatumtijd": tmax.strftime('%Y-%m-%dT%H:%M:%S.000+01:00')}
    return periode


def measurementlist2df(measurementlist):
    df = pd.DataFrame(measurementlist)
    df['Tijdstip'] = pd.to_datetime(df['Tijdstip'])
    df = df.set_index('Tijdstip')
    df['Meetwaarde'] = [x['Waarde_Numeriek'] for x in df['Meetwaarde']]
    df.loc[df['Meetwaarde'] == 999999999.0, 'Meetwaarde'] = np.NaN
    if True:
        for column in df.iloc[0]['WaarnemingMetadata']:
            df[column] = [x[column][0] for x in df['WaarnemingMetadata']]
    return df


def flatten(json, veld):
    for key in json[veld]:
        if key in json:
            raise(Exception('Key {} allready existst'.format(key)))
        json[key] = json[veld][key]
    json.pop(veld)
    return json


def ophalen_waarnemingen(locatie, grootheid='WATHTE', tmin=None, tmax=None,
                         raw_response=False):
    """Get some measurements

    Parameters
    ----------
    locatie : dict
        A dictionary with location information, with the keys X, Y and Code
    grootheid : string
        the code for the variable you want measurements from
    tmin : string or pandas.Timestamp
        The start-date for which you which to get measurements
    tmax : string or pandas.Timestamp
        The end-date for which you which to get measurements

    """
    if not isinstance(locatie, dict):
        locatie = get_locatie(locatie)
    url = 'https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/OphalenWaarnemingen'
    periode = get_periode(tmin, tmax)
    AquoMetadata = {"Grootheid": {"Code": grootheid}}
    json_input = {"AquoPlusWaarnemingMetadata": {"AquoMetadata": AquoMetadata},
                  "Locatie": locatie,
                  "Periode": periode}
    json_output = post_data(url, json_input)
    if raw_response:
        return json_output
    if len(json_output['WaarnemingenLijst']) > 1:
        raise(Exception('Meerdere locaties nog niet ondersteund'))
    df = measurementlist2df(
        json_output['WaarnemingenLijst'][0]['MetingenLijst'])
    return df


def ophalen_laatste_waarnemingen(locatie_lijst, grootheid='WATHTE',
                                 raw_response=False):
    """Get the latest measurements

    Parameters
    ----------
    locatie_lijst : list of dicts, pandas DataFrame or pandas Series
        A list of dictionaries with location information, with the keys X, Y and Code
    grootheid : string
        the code for the variable you want measurements from
    raw_response : bool, optional
        return the raw json-response instead of a pandas DataFrame (defaults to
        False)

    """
    if not isinstance(locatie_lijst, list):
        locatie_lijst = get_locatie_lijst(locatie_lijst)
    url = 'https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/OphalenLaatsteWaarnemingen'
    AquoMetadata = {"Grootheid": {"Code": grootheid}}
    json_input = {"AquoPlusWaarnemingMetadataLijst": [{'AquoMetadata': AquoMetadata}],
                  "LocatieLijst": locatie_lijst}
    json_output = post_data(url, json_input)
    if raw_response:
        return json_output

    for waarnemingen in json_output['WaarnemingenLijst']:
        flatten(waarnemingen, 'Locatie')
        flatten(waarnemingen, 'AquoMetadata')

    df = pd.DataFrame(json_output['WaarnemingenLijst'])
    df['MetingenLijst'] = [measurementlist2df(x) for x in df['MetingenLijst']]
    return df


def check_waarnemingen_aanwezig(locatie_lijst, grootheid='WATHTE', tmin=None,
                                tmax=None, raw_response=False):
    if not isinstance(locatie_lijst, list):
        locatie_lijst = get_locatie_lijst(locatie_lijst)
    url = 'https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/CheckWaarnemingenAanwezig'
    periode = get_periode(tmin, tmax)
    json_input = {"AquoMetadataLijst": [{"Grootheid": {"Code": grootheid}}],
                  "LocatieLijst": locatie_lijst,
                  "Periode": periode}
    json_output = post_data(url, json_input)
    if raw_response:
        return json_output
    return json_output['WaarnemingenAanwezig'] == 'true'


def ophalen_aantal_waarnemingen(locatie_lijst, Groeperingsperiode='Week',
                                grootheid='WATHTE', tmin=None, tmax=None):
    if not isinstance(locatie_lijst, list):
        locatie_lijst = get_locatie_lijst(locatie_lijst)
    url = 'https://waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES_DBO/OphalenAantalWaarnemingen'
    periode = get_periode(tmin, tmax)
    json_input = {"AquoMetadataLijst": [{"Grootheid": {"Code": grootheid}}],
                  "Groeperingsperiode": Groeperingsperiode,
                  "LocatieLijst": locatie_lijst,
                  "Periode": periode}
    json_output = post_data(url, json_input)
    return json_output


def aanvragen_bulk_waarnemingen(locatie_lijst, email_to, grootheid='WATHTE',
                                tmin=None, tmax=None, email_from="info@rws.nl"):
    if not isinstance(locatie_lijst, list):
        locatie_lijst = get_locatie_lijst(locatie_lijst)
    url = 'https://waterwebservices.rijkswaterstaat.nl/BULKWAARNEMINGSERVICES_DBO/AanvragenBulkWaarnemingen'
    AquoMetadataLijst = [
        {"Grootheid": {"Code": grootheid}, "Eenheid": {"Code": "cm"}}]
    periode = get_periode(tmin, tmax)

    json_input = {"Zoekvraag": {"AquoMetadataLijst": AquoMetadataLijst,
                                "LocatieLijst": locatie_lijst,
                                "Periode": periode},
                  "Email_succes": {"From": email_from,
                                   "To": email_to,
                                   "Subject": "Aanvraag bestand waarnemingen",
                                   "Body": "Uw bestand met waarnemingen kunt u downloaden via {link_bestand}."},
                  "Email_fout": {"From": email_from,
                                 "To": email_to,
                                 "Subject": "Aanvraag niet gelukt",
                                 "Body": "Uw aanvraag voor het bestand met waarnemingen is mislukt."},
                  "Email_bevestiging": {"From": email_from,
                                        "To": email_to,
                                        "Subject": "Bevestiging van aanvraag",
                                        "Body": "Uw aanvraag is ontvangen. U ontvangt binnen 24 uur een e-mail met daarin een link voor het downloaden van de aanvraag."}}
    json_output = post_data(url, json_input)
    return json_output


def aanvragen_bulk_waarnemingen_in_extent(extent, email_to, grootheid='WATHTE',
                                          tmin=None, tmax=None, separate=True):

    catalogus = get_catalogus()
    locations = get_locations(catalogus, extent=extent, grootheid=grootheid)

    # %%only keep locations with measurements
    locations['keep'] = False
    for ind in tqdm(locations.index, desc='Checking if there are measurements'):
        if check_waarnemingen_aanwezig(locations.loc[ind], grootheid=grootheid,
                                       tmin=tmin, tmax=tmax):
            locations.loc[ind, 'keep'] = True
    locations = locations[locations['keep']]

    # %% generate a download per location for a specific tmin and tmax
    # download all measurements between:
    if separate:
        json_output_list = []
        for i in locations.index:
            jo = aanvragen_bulk_waarnemingen(locations.loc[i],
                                             email_to,
                                             tmin=tmin, tmax=tmax)
            json_output_list.append(jo)
        return json_output_list
    else:
        jo = aanvragen_bulk_waarnemingen(email_to, locatie_lijst=locations,
                                         tmin=tmin, tmax=tmax)
        return jo


def post_data(url, json_input):
    resp = requests.post(url, json=json_input)
    if not resp.ok:
        raise(Exception('{} failed'.format(url.split('/')[-1])))
    json_output = resp.json()
    if not json_output['Succesvol']:
        raise(Exception(json_output['Foutmelding']))
    return json_output


def read_waterinfo_zip(fname, **kwargs):
    """Read a waterinfo zip-file which contains one csv-file"""
    if fname.endswith('.zip'):
        pathname, file = os.path.split(fname)
        name = os.path.splitext(file)[0]
        zf = zipfile.ZipFile(os.path.join(pathname, file))
        df = read_waterinfo_csv(zf.open('{}.csv'.format(name)), **kwargs)
        return df
    else:
        raise(Exception('File has to be a zip-file'))


def read_waterinfo_csv(fname, freq=None, pivot=True, metadata=False,
                       columns='MEETPUNT_IDENTIFICATIE'):
    """Read a waterinfo csv-file"""
    df = pd.read_csv(fname,
                     sep=';',
                     decimal=',',
                     encoding="ISO-8859-1",
                     parse_dates=[['WAARNEMINGDATUM', 'WAARNEMINGTIJD']],
                     dayfirst=True,
                     infer_datetime_format=True,
                     index_col='WAARNEMINGDATUM_WAARNEMINGTIJD',
                     na_values=[999999999, -999999999])

    # from cm to m
    mask = df['EENHEID_CODE'] == 'cm'
    if mask.any():
        df.loc[mask, 'NUMERIEKEWAARDE'] = df.loc[mask, 'NUMERIEKEWAARDE'] / 100.
        df.loc[mask, 'EENHEID_CODE'] = 'm'

    if metadata:
        meta = df.groupby(columns).last()
        assert np.all(meta['EPSG'] == 25831), 'Unexpected reference system'
        meta['x'], meta['y'] = transform(Proj("epsg:25831"),
                                         Proj("epsg:28992"),
                                         meta['X'].values,
                                         meta['Y'].values)
        meta = df2gdf(meta)

    if pivot:
        df = df.pivot_table(index='WAARNEMINGDATUM_WAARNEMINGTIJD',
                            columns=columns,
                            values='NUMERIEKEWAARDE')
        if freq is not None:
            df = df.resample(freq, label='right', closed='right').mean()
    if metadata:
        return df, meta
    return df


def df2gdf(df, xcol='x', ycol='y'):
    """Make a GeoDataFrame from a DataFrame, assuming the geometry are points"""
    geometry = [Point((s[xcol], s[ycol])) for i, s in df.iterrows()]
    gdf = gpd.GeoDataFrame(df.copy(), geometry=geometry)
    return gdf


def get_river_polygons(extent=None, owmnaam=None):
    typename = 'kaderrichtlijn_water:oppervlaktewater_lichamen_v'
    gdf = get_krw_wfs_features(typename, extent=extent, owmnaam=owmnaam)
    return gdf


def get_river_lines(extent=None, owmnaam=None):
    typename = 'kaderrichtlijn_water:oppervlaktewater_lichamen_l'
    gdf = get_krw_wfs_features(typename, extent=extent, owmnaam=owmnaam)
    return gdf


def get_krw_wfs_features(typename, extent=None, owmnaam=None):
    wfs = get_krw_wfs()
    kwargs = {'typename': typename, 'outputFormat': 'json'}
    if extent is not None:
        kwargs['bbox'] = [extent[0], extent[2], extent[1], extent[3]]
    if owmnaam is not None:
        if extent is not None:
            raise(Exception('Either supply extent or ownnaam, not both'))
        filter = PropertyIsEqualTo(propertyname='OWMNAAM', literal=owmnaam)
        kwargs['filter'] = etree.tostring(filter.toXML()).decode("utf-8")
    response = wfs.getfeature(**kwargs).read()
    gdf = gpd.GeoDataFrame.from_features(json.loads(response)['features'])
    gdf = gdf.set_index('OWMNAAM')
    return gdf


def get_krw_wfs():
    url = 'https://geoservices.rijkswaterstaat.nl/apps/geoserver/kaderrichtlijn_water/ows?service=WFS'
    return WebFeatureService(url, version='2.0.0')

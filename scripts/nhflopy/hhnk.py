# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 15:18:06 2020

@author: Artesia
"""
import requests
import pandas as pd
import geopandas as gpd

def get_filters(documentVersion='1.25', documentFormat='PI_JSON'): 
    # format can also be 'PI_XML'
    # documentVersion
    url = 'https://fews.hhnk.nl/FewsPiService/rest/fewspiservice/v1/filters'
    params = {'documentVersion': documentVersion,
              'documentFormat': documentFormat}
    r = requests.get(url, params=params)
    if r.status_code not in [200]:
        raise(Exception(r.text))
    return r.json()['filters'][0]
        
        
def get_locations(filterId='HHNK_Dijkmon_WEB', showAttributes=False,
                  documentVersion='1.25', documentFormat='PI_JSON'):
    url = 'https://fews.hhnk.nl/FewsPiService/rest/fewspiservice/v1/locations'
    params = {'filterId':filterId,
              'showAttributes': showAttributes,
              'documentVersion': documentVersion,
              'documentFormat': documentFormat}
    
    r = requests.get(url, params=params)
    if r.status_code not in [200]:
        raise(Exception(r.text))
    # make a GeoDataFrame
    df = pd.DataFrame(r.json()['locations']).set_index('locationId')
    gdf = result_to_gdf(df)
    return gdf

def get_timeseries(locationIds=None, parameterIds='WATHTE [m][NAP][OW]',
                   filterId='HHNK_Dijkmon_WEB', omitMissing=True,
                   startTime=None, endTime=None,
                   documentVersion='1.25', documentFormat='PI_JSON', **kwargs):
    url = 'https://fews.hhnk.nl/FewsPiService/rest/fewspiservice/v1/timeseries'
    params = {'filterId': filterId,
              'omitMissing': omitMissing,
              'documentVersion': documentVersion,
              'documentFormat': documentFormat}
    if locationIds is not None:
        params['locationIds'] = locationIds
    if parameterIds is not None:
        params['parameterIds'] = parameterIds
    if startTime is not None:
        startTime = pd.to_datetime(startTime)
        params['startTime'] = startTime.strftime('%Y-%m-%dT%H:%M:%SZ')
    if endTime is not None:
        endTime = pd.to_datetime(endTime)
        params['endTime'] = endTime.strftime('%Y-%m-%dT%H:%M:%SZ')
    params.update(kwargs)
    r = requests.get(url, params=params)
    if r.status_code not in [200]:
        raise(Exception(r.text))
    tss = r.json()['timeSeries']
    df = pd.DataFrame([x['header'] for x in tss]).set_index('locationId')
    events = []
    for ts in tss:
        if 'events' in ts:
            event_df = pd.DataFrame(ts['events'])
            for col in ['value', 'flag']:
                if col in event_df:
                    event_df[col] = pd.to_numeric(event_df[col])
            dt = pd.to_datetime(event_df['date'] + ' ' + event_df['time'])
            event_df = event_df.set_index(dt).drop(columns=['date','time'])
            events.append(event_df)
        else:
            events.append(None)
    df['events'] = events
    gdf = result_to_gdf(df)
    return gdf

def result_to_gdf(df):
    for col in ['lat','lon','x','y','z']:
        if col in df:
            df[col] = pd.to_numeric(df[col])
    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['x'], df['y']))

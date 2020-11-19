# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:50:48 2020

@author: oebbe
"""
import numpy as np
import pandas as pd
import xarray as xr


def get_model_ds_time(model_name, start_time, steady_state, steady_start,
                      time_units='DAYS', transient_timesteps=0,
                      perlen=1.0,
                      nstp=1, tsmult=1.0):
    """ Get a model dataset with the time variant data

    Parameters
    ----------
    model_name : str
        name of the model. Cannot have more than 16 characters (modflow)
    start_time : str or datetime
        start time of the model.
    steady_state : bool
        if True the model is steady state with one time step.
    steady_start : bool
        if True the model is transient with a steady state start time step.
    time_units : str, optional
        time unit of the model. The default is 'DAYS'.
    transient_timesteps : int, optional
        number of transient time steps. The default is 0.
    perlen : TYPE, optional
        DESCRIPTION. The default is 1.0.
    nstp : TYPE, optional
        DESCRIPTION. The default is 1.
    tsmult : TYPE, optional
        DESCRIPTION. The default is 1.0.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with time variant model data
    
    """
    if len(model_name)>16:
        raise ValueError('model_name can not have more than 16 characters')
    
    start_time_dt = pd.to_datetime(start_time)

    if steady_state:
        nper = 1
        time_dt = [pd.to_datetime(start_time)]
    elif steady_start:
        nper = 1 + transient_timesteps
        time_dt = pd.date_range(
            start_time_dt - pd.to_timedelta(perlen, unit=time_units),
            start_time_dt + pd.to_timedelta((transient_timesteps - 1) * perlen,
                                            unit=time_units), periods=nper)
        timedelta = pd.to_timedelta(perlen, unit=time_units)
        start_time_dt = start_time_dt - timedelta
    else:
        nper = transient_timesteps
        time_dt = pd.date_range(
            start_time_dt,
            start_time_dt + pd.to_timedelta((transient_timesteps - 1) * perlen,
                                            unit=time_units), periods=nper)

    time_steps = list(range(nper))

    model_ds = xr.Dataset(coords={'time': time_dt})
    
    model_ds['time_steps'] = xr.DataArray(time_steps, dims=('time'),
                                          coords={'time': model_ds.time})
    model_ds.attrs['model_name'] = model_name
    model_ds.attrs['time_units'] = time_units
    model_ds.attrs['start_time'] = str(start_time_dt)
    model_ds.attrs['nper'] = nper
    model_ds.attrs['perlen'] = perlen
    model_ds.attrs['nstp'] = nstp
    model_ds.attrs['tsmult'] = tsmult

    # netcdf files cannot handle booleans
    model_ds.attrs['steady_start'] = int(steady_start)
    model_ds.attrs['steady_state'] = int(steady_state)

    return model_ds

def get_tdis_perioddata(model_ds):
    """ Get tdis_perioddata from model_ds

    Parameters
    ----------
    model_ds : xarray.Dataset
        dataset with time variant model data

    Raises
    ------
    NotImplementedError
        cannot handle timesteps with variable step length.

    Returns
    -------
    tdis_perioddata : [perlen, nstp, tsmult]
        * perlen (double) is the length of a stress period.
        * nstp (integer) is the number of time steps in a stress period.
        * tsmult (double) is the multiplier for the length of successive time
          steps. The length of a time step is calculated by multiplying the
          length of the previous time step by TSMULT. The length of the first
          time step, :math:`\Delta t_1`, is related to PERLEN, NSTP, and
          TSMULT by the relation :math:`\Delta t_1= perlen \frac{tsmult -
          1}{tsmult^{nstp}-1}`.

    """
    
    try:
        float(model_ds.perlen)
        tdis_perioddata = [(model_ds.perlen, model_ds.nstp, model_ds.tsmult)] * int(model_ds.nper)
    except:
        raise NotImplementedError('variable time step length not yet implemented')
    
    # netcdf does not support multi-dimensional array attributes
    #model_ds.attrs['tdis_perioddata'] = tdis_perioddata
    
    return tdis_perioddata
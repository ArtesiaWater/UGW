import pandas as pd
import xarray as xr


def get_model_ts(start_time, steady_state, steady_start,
                 time_units='DAYS', transient_timesteps=0,
                 perlen=1.0, nstp=1, tsmult=1.0):
    """ Get a model dataset with the time variant data

    Parameters
    ----------
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
    model_ts : xarray.Dataset
        dataset with time variant model data
    """

    # hier worden waarden berekend pas dit niet aan
    start_time_dt = pd.to_datetime(start_time)

    if steady_state:
        nper = 1
        time_dt = [pd.to_datetime(start_time)]
    elif steady_start:
        nper = 1 + transient_timesteps
        time_dt = pd.date_range(
            start_time_dt - pd.to_timedelta(1, unit=time_units),
            start_time_dt + pd.to_timedelta((transient_timesteps - 1) * perlen,
                                            unit=time_units), periods=nper)

    else:
        nper = transient_timesteps
        time_dt = pd.date_range(
            start_time_dt,
            start_time_dt + pd.to_timedelta((transient_timesteps - 1) * perlen,
                                            unit=time_units), periods=nper)

    time_steps = list(range(nper))

    model_ts = xr.Dataset(coords={'time': time_dt})
    model_ts['time_steps'] = xr.DataArray(time_steps, dims=('time'),
                                          coords={'time': model_ts.time})
    model_ts.attrs['time_units'] = time_units
    model_ts.attrs['start_time'] = start_time
    model_ts.attrs['nper'] = nper
    model_ts.attrs['perlen'] = perlen
    model_ts.attrs['nstp'] = nstp
    model_ts.attrs['tsmult'] = tsmult

    # netcdf files cannot handle booleans
    model_ts.attrs['steady_start'] = int(steady_start)
    model_ts.attrs['steady_state'] = int(steady_state)

    return model_ts

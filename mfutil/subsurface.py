import numpy as np
import xarray as xr

import mgrid


def add_kh_kv_from_regis_to_dataset(regis_ds, model_ds, anisotropy,
                                    fill_value_kh, fill_value_kv,
                                    verbose=False):
    """ add kh and kv from regis to a model dataset.

    Supports structured and unstructured grids.

    Parameters
    ----------
    regis_ds : xarray.Dataset
        dataset with regis data with kh and kv
    model_ds : xarray.Dataset
        dataset with model data where kh and kv are added to
    anisotropy : int or float
        factor to calculate kv from kh or the other way around
    fill_value_kh : int or float, optional
        use this value for kh if there is no data in regis. The default is 1.0.
    fill_value_kv : int or float, optional
        use this value for kv if there is no data in regis. The default is 1.0.
    verbose : bool, optional
        print additional information. default is False

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data with new kh and kv

    Notes
    -----
    regis also has 'c' and 'kd' values. These are ignored at the moment
    """
    model_ds.attrs['anisotropy'] = anisotropy
    model_ds.attrs['fill_value_kh'] = fill_value_kh
    model_ds.attrs['fill_value_kv'] = fill_value_kv

    kh_arr = regis_ds['kh'].data
    kv_arr = regis_ds['kv'].data

    if verbose:
        print('add kh and kv from regis to modflow model')

    kh, kv = get_kh_kv(kh_arr, kv_arr, anisotropy,
                       fill_value_kh=fill_value_kh,
                       fill_value_kv=fill_value_kv,
                       verbose=verbose)

    model_ds['kh'] = xr.ones_like(model_ds['idomain']) * kh

    model_ds['kv'] = xr.ones_like(model_ds['idomain']) * kv

    return model_ds


def get_kh_kv(kh_in, kv_in, anisotropy,
              fill_value_kh=1.0, fill_value_kv=1.0,
              verbose=False):
    """maak kh en kv rasters voor flopy vanuit een regis raster met nan
    waardes.

    vul kh raster door:
    1. pak kh uit regis, tenzij nan dan:
    2. pak kv uit regis vermenigvuldig met anisotropy, tenzij nan dan:
    3. pak fill_value_kh

    vul kv raster door:
    1. pak kv uit regis, tenzij nan dan:
    2. pak kh uit regis deel door anisotropy, tenzij nan dan:
    3. pak fill_value_kv

    Supports structured and unstructured grids.

    Parameters
    ----------
    kh_in : np.ndarray
        kh from regis with nan values shape(nlay, nrow, ncol) or 
        shape(nlay, len(cid))
    kv_in : np.ndarray
        kv from regis with nan values shape(nlay, nrow, ncol) or 
        shape(nlay, len(cid))
    anisotropy : int or float
        factor to calculate kv from kh or the other way around
    fill_value_kh : int or float, optional
        use this value for kh if there is no data in regis. The default is 1.0.
    fill_value_kv : int or float, optional
        use this value for kv if there is no data in regis. The default is 1.0.
    verbose : bool, optional
        print additional information. default is False

    Returns
    -------
    kh_out : np.ndarray
        kh without nan values (nlay, nrow, ncol) or shape(nlay, len(cid))
    kv_out : np.ndarray
        kv without nan values (nlay, nrow, ncol) or shape(nlay, len(cid))
    """
    kh_out = np.zeros_like(kh_in)
    for i, kh_lay in enumerate(kh_in):
        kh_new = kh_lay.copy()
        kv_new = kv_in[i].copy()
        if ~np.all(np.isnan(kh_new)):
            if verbose:
                print(f'layer {i} has a kh')
            kh_out[i] = np.where(np.isnan(kh_new), kv_new * anisotropy, kh_new)
            kh_out[i] = np.where(np.isnan(kh_out[i]), fill_value_kh, kh_out[i])
        elif ~np.all(np.isnan(kv_new)):
            if verbose:
                print(f'layer {i} has a kv')
            kh_out[i] = np.where(
                np.isnan(kv_new), fill_value_kh, kv_new * anisotropy)
        else:
            if verbose:
                print(f'kv and kh both undefined in layer {i}')
            kh_out[i] = fill_value_kh

    kv_out = np.zeros_like(kv_in)
    for i, kv_lay in enumerate(kv_in):
        kv_new = kv_lay.copy()
        kh_new = kh_in[i].copy()
        if ~np.all(np.isnan(kv_new)):
            if verbose:
                print(f'layer {i} has a kv')
            kv_out[i] = np.where(np.isnan(kv_new), kh_new / anisotropy, kv_new)
            kv_out[i] = np.where(np.isnan(kv_out[i]), fill_value_kv, kv_out[i])
        elif ~np.all(np.isnan(kh_new)):
            if verbose:
                print(f'layer {i} has a kh')
            kv_out[i] = np.where(
                np.isnan(kh_new), fill_value_kv, kh_new / anisotropy)
        else:
            if verbose:
                print(f'kv and kh both undefined in layer {i}')
            kv_out[i] = fill_value_kv

    return kh_out, kv_out


def add_top_bot_to_model_ds(regis_ds, model_ds,
                            nodata=None,
                            max_per_nan_bot=50,
                            gridtype='structured',
                            verbose=False):
    """ add top and bot from a regis dataset to a model dataset.

    Supports structured and unstructured grids.

    Parameters
    ----------
    regis_ds : xarray.Dataset
        dataset with regis data with a top and bottom
    model_ds : xarray.Dataset
        dataset with model data where top and bot are added to
    nodata : int, optional
        if the first_active_layer data array in model_ds has this value,
        it means this cell is inactive in all layers. If nodata is None the
        nodata value in model_ds is used.
        the default is None
    max_per_nan_bot : int or float, optional
        if the percentage of cells that have nan values in all layers is
        higher than this an error is raised. The default is 50.
    gridtype : str, optional
        type of grid, options are 'structured' and 'unstructured'. 
        The default is 'structured'.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data including top and bottom

    """
    if nodata is None:
        nodata = model_ds.attrs['nodata']

    if verbose:
        print('using top and bottom from regis for modflow model')

    if gridtype == 'structured':

        model_ds = add_top_bot_structured(regis_ds, model_ds,
                                          nodata=nodata,
                                          max_per_nan_bot=max_per_nan_bot)

    elif gridtype == 'unstructured':
        model_ds = add_top_bot_unstructured(regis_ds, model_ds,
                                            nodata=nodata,
                                            max_per_nan_bot=max_per_nan_bot)

    return model_ds


def add_top_bot_unstructured(regis_ds, model_ds, nodata=-999,
                             max_per_nan_bot=50):
    """ voeg top en bottom vanuit regis toe aan de model dataset

    Deze functie is bedoeld voor unstructured arrays in modflow 6

    stappen:
    1. Zorg dat de onderste laag altijd een bodemhoogte heeft, als de bodem
    van alle bovenliggende lagen nan is, pak dan 0.
    2. Zorg dat de top van de bovenste laag altijd een waarde heeft, als de
    top van alle onderligende lagen nan is, pak dan 0.
    3. Vul de nan waarden in alle andere lagen door:
        a. pak bodem uit regis, tenzij nan dan:
        b. gebruik bodem van de laag erboven (of de top voor de bovenste laag)

    Supports only unstructured grids.

    Parameters
    ----------
    regis_ds : xarray.Dataset
        dataset with regis data with a top and bottom
    model_ds : xarray.Dataset
        dataset with model data where top and bottom are added to
    nodata : int, optional
        if the first_active_layer data array in model_ds has this value,
        it means this cell is inactive in all layers 
    max_per_nan_bot : int or float, optional
        if the percentage of cells that have nan values in all layers. 
        The default is 50.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data including top and bottom
    """
    # step 1:
    # set nan-value in bottom array
    # set to zero if value is nan in all layers
    # set to minimum value of all layers if there is any value in any layer
    active_domain = model_ds['first_active_layer'].data != nodata

    lowest_bottom = regis_ds['bottom'].data[-1].copy()
    if np.any(active_domain == False):
        percentage = 100 * (active_domain == False).sum() / \
            (active_domain.shape[0])
        if percentage > max_per_nan_bot:
            print(f"{percentage:0.1f}% cells with nan values in every layer"
                  " check your extent")
            raise MemoryError('if you want to'
                              ' ignore this error set max_per_nan_bot higher')

        # set bottom to zero if bottom in a cell is nan in all layers
        lowest_bottom = np.where(active_domain, lowest_bottom, 0)

    if np.any(np.isnan(lowest_bottom)):
        # set bottom in a cell to lowest bottom of all layers
        i_nan = np.where(np.isnan(lowest_bottom))
        for i in i_nan:
            val = np.nanmin(regis_ds['bottom'].data[:, i])
            lowest_bottom[i] = val
            if np.isnan(val):
                raise ValueError(
                    'this should never happen please contact Onno')

    # step 2: get highest top values of all layers without nan values
    highest_top = regis_ds['top'].data[0].copy()
    if np.any(np.isnan(highest_top)):
        highest_top = np.where(active_domain, highest_top, 0)

    if np.any(np.isnan(highest_top)):
        i_nan = np.where(np.isnan(highest_top))
        for i in i_nan:
            val = np.nanmax(regis_ds['top'].data[:, i])
            highest_top[i] = val
            if np.isnan(val):
                raise ValueError(
                    'this should never happen please contact Onno')

    # step 3: fill nans in all layers
    nlay = model_ds.dims['layer']
    top_bot_raw = np.ones((nlay + 1, model_ds.dims['cid']))
    top_bot_raw[0] = highest_top
    top_bot_raw[1:-1] = regis_ds['bottom'].data[:-1].copy()
    top_bot_raw[-1] = lowest_bottom
    top_bot = np.ones_like(top_bot_raw)
    for i_from_bot, blay in enumerate(top_bot_raw[::-1]):
        i_from_top = nlay - i_from_bot
        new_lay = blay.copy()
        if np.any(np.isnan(new_lay)):
            lay_from_bot = i_from_bot
            lay_from_top = nlay - lay_from_bot
            while np.any(np.isnan(new_lay)):
                new_lay = np.where(np.isnan(new_lay),
                                   top_bot_raw[lay_from_top],
                                   new_lay)
                lay_from_bot += 1
                lay_from_top = nlay - lay_from_bot

        top_bot[i_from_top] = new_lay

    model_ds['bot'] = xr.DataArray(top_bot[1:], dims=('layer', 'cid'),
                                   coords={'cid': model_ds.cid.data,
                                           'layer': model_ds.layer.data})
    model_ds['top'] = xr.DataArray(top_bot[0], dims=('cid'),
                                   coords={'cid': model_ds.cid.data})

    return model_ds


def add_top_bot_structured(regis_ds, model_ds, nodata=-999,
                           max_per_nan_bot=50):
    """Voeg top en bottom vanuit regis toe aan de model dataset

    Deze functie is bedoeld voor structured arrays in modflow 6

    stappen:
    1. Zorg dat de onderste laag altijd een bodemhoogte heeft, als de bodem
    van alle bovenliggende lagen nan is, pak dan 0.
    2. Zorg dat de top van de bovenste laag altijd een waarde heeft, als de
    top van alle onderligende lagen nan is, pak dan 0.
    3. Vul de nan waarden in alle andere lagen door:
        a. pak bodem uit regis, tenzij nan dan:
        b. gebruik bodem van de laag erboven (of de top voor de bovenste laag)

    Supports only structured grids.

    Parameters
    ----------
    regis_ds : xarray.Dataset
        dataset with regis data with a top and bottom
    model_ds : xarray.Dataset
        dataset with model data where top and bottom are added to
    nodata : int, optional
        if the first_active_layer data array in model_ds has this value,
        it means this cell is inactive in all layers 
    max_per_nan_bot : int or float, optional
        if the percentage of cells that have nan values in all layers. 
        The default is 50.

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data including top and bottom

    """

    active_domain = model_ds['first_active_layer'].data != nodata

    # step 1:
    # set nan-value in bottom array
    # set to zero if value is nan in all layers
    # set to minimum value of all layers if there is any value in any layer
    lowest_bottom = regis_ds['bottom'].data[-1].copy()
    if np.any(active_domain == False):
        percentage = 100 * (active_domain == False).sum() / \
            (active_domain.shape[0] * active_domain.shape[1])
        if percentage > max_per_nan_bot:
            print(f"{percentage:0.1f}% cells with nan values in every layer"
                  " check your extent")
            raise MemoryError('if you want to'
                              ' ignore this error set max_per_nan_bot higher')

        # set bottom to zero if bottom in a cell is nan in all layers
        lowest_bottom = np.where(active_domain, lowest_bottom, 0)

    if np.any(np.isnan(lowest_bottom)):
        # set bottom in a cell to lowest bottom of all layers
        rc_nan = np.where(np.isnan(lowest_bottom))
        for row, col in zip(rc_nan[0], rc_nan[1]):
            val = np.nanmin(regis_ds['bottom'].data[:, row, col])
            lowest_bottom[row, col] = val
            if np.isnan(val):
                raise ValueError(
                    'this should never happen please contact Onno')

    # step 2: get highest top values of all layers without nan values
    highest_top = regis_ds['top'].data[0].copy()
    if np.any(np.isnan(highest_top)):
        # set top to zero if top in a cell is nan in all layers
        highest_top = np.where(active_domain, highest_top, 0)

    if np.any(np.isnan(highest_top)):
        # set top in a cell to highest top of all layers
        rc_nan = np.where(np.isnan(highest_top))
        for row, col in zip(rc_nan[0], rc_nan[1]):
            val = np.nanmax(regis_ds['top'].data[:, row, col])
            highest_top[row, col] = val
            if np.isnan(val):
                raise ValueError(
                    'this should never happen please contact Onno')

    # step 3: fill nans in all layers
    nlay = model_ds.dims['layer']
    nrow = model_ds.dims['y']
    ncol = model_ds.dims['x']
    top_bot_raw = np.ones((nlay + 1, nrow, ncol))
    top_bot_raw[0] = highest_top
    top_bot_raw[1:-1] = regis_ds['bottom'].data[:-1].copy()
    top_bot_raw[-1] = lowest_bottom
    top_bot = np.ones_like(top_bot_raw)
    for i_from_bot, blay in enumerate(top_bot_raw[::-1]):
        i_from_top = nlay - i_from_bot
        new_lay = blay.copy()
        if np.any(np.isnan(new_lay)):
            lay_from_bot = i_from_bot
            lay_from_top = nlay - lay_from_bot
            while np.any(np.isnan(new_lay)):
                new_lay = np.where(np.isnan(new_lay),
                                   top_bot_raw[lay_from_top],
                                   new_lay)
                lay_from_bot += 1
                lay_from_top = nlay - lay_from_bot

        top_bot[i_from_top] = new_lay

    model_ds['bot'] = xr.DataArray(top_bot[1:], dims=('layer', 'y', 'x'),
                                   coords={'x': model_ds.x.data,
                                           'y': model_ds.y.data,
                                           'layer': model_ds.layer.data})

    model_ds['top'] = xr.DataArray(top_bot[0], dims=('y', 'x'),
                                   coords={'x': model_ds.x.data,
                                           'y': model_ds.y.data})

    return model_ds


def fill_top_bot_kh_kv_at_mask(model_ds, fill_mask,
                               gridtype='structured',
                               gridprops=None):
    """ fill values in top, bot, kh and kv where:
        1. the cell is True in mask
        2. the cell thickness is greater than 0

    fill values:
        top: 0
        bot: minimum of bottom_filled or top
        kh: kh_filled if thickness is greater than 0
        kv: kv_filled if thickness is greater than 0

    Parameters
    ----------
    model_ds : xr.DataSet
        model dataset, should contain 'first_active_layer'
    fill_mask : xr.DataArray
        1 where a cell should be replaced by masked value.
    gridtype : str, optional
        type of grid.        
    gridprops : dictionary, optional
        dictionary with grid properties output from gridgen. Default is None

    Returns
    -------
    model_ds : xr.DataSet
        model dataset with adjusted data variables:
            'top', 'bot', 'kh', 'kv'
    """

    # zee cellen hebben altijd een top gelijk aan 0
    model_ds['top'] = xr.where(fill_mask, 0, model_ds['top'])

    if gridtype == 'structured':
        fill_function = mgrid.fillnan_dataarray_structured_grid
        fill_function_kwargs = {}
    elif gridtype == 'unstructured':
        fill_function = mgrid.fillnan_dataarray_unstructured_grid
        fill_function_kwargs = {'gridprops': gridprops}

    for lay in range(model_ds.dims['layer']):
        bottom_nan = xr.where(fill_mask, np.nan, model_ds['bot'][lay])
        bottom_filled = fill_function(bottom_nan, **fill_function_kwargs)[0]

        kh_nan = xr.where(fill_mask, np.nan, model_ds['kh'][lay])
        kh_filled = fill_function(kh_nan, **fill_function_kwargs)[0]

        kv_nan = xr.where(fill_mask, np.nan, model_ds['kv'][lay])
        kv_filled = fill_function(kv_nan, **fill_function_kwargs)[0]

        if lay == 0:
            # top ligt onder bottom_filled -> laagdikte wordt 0
            # top ligt boven bottom_filled -> laagdikte o.b.v. bottom_filled
            mask_top = model_ds['top'] < bottom_filled
            model_ds['bot'][lay] = xr.where(fill_mask * mask_top,
                                            model_ds['top'],
                                            bottom_filled)
            model_ds['kh'][lay] = xr.where(fill_mask * mask_top,
                                           model_ds['kh'][lay],
                                           kh_filled)
            model_ds['kv'][lay] = xr.where(fill_mask * mask_top,
                                           model_ds['kv'][lay],
                                           kv_filled)

        else:
            # top ligt onder bottom_filled -> laagdikte wordt 0
            # top ligt boven bottom_filled -> laagdikte o.b.v. bottom_filled
            mask_top = model_ds['bot'][lay - 1] < bottom_filled
            model_ds['bot'][lay] = xr.where(fill_mask * mask_top,
                                            model_ds['bot'][lay - 1],
                                            bottom_filled)
            model_ds['kh'][lay] = xr.where(fill_mask * mask_top,
                                           model_ds['kh'][lay],
                                           kh_filled)
            model_ds['kv'][lay] = xr.where(fill_mask * mask_top,
                                           model_ds['kv'][lay],
                                           kv_filled)

    return model_ds


def add_bathymetry_to_top_bot_kh_kv(model_ds, bathymetry,
                                    fill_mask,
                                    kh_sea=10,
                                    kv_sea=10):
    """ add bathymetry to the top and bot of each layer for all cells with
    fill_mask.

    Parameters
    ----------
    model_ds : xarray.Dataset
        dataset with model data, should 
    fill_mask : xr.DataArray
        cell value is 1 if you want to add bathymetry

    Returns
    -------
    model_ds : xarray.Dataset
        dataset with model data where the top, bot, kh and kv are changed

    """
    model_ds['top'] = xr.where(fill_mask,
                               0.0,
                               model_ds['top'])

    lay = 0
    model_ds['bot'][lay] = xr.where(fill_mask,
                                    bathymetry,
                                    model_ds['bot'][lay])

    model_ds['kh'][lay] = xr.where(fill_mask,
                                   kh_sea,
                                   model_ds['kh'][lay])

    model_ds['kv'][lay] = xr.where(fill_mask,
                                   kv_sea,
                                   model_ds['kv'][lay])

    # reset bot for all layers based on bathymetrie
    for lay in range(1, model_ds.dims['layer']):
        model_ds['bot'][lay] = np.where(model_ds['bot'][lay] > model_ds['bot'][lay - 1],
                                        model_ds['bot'][lay - 1],
                                        model_ds['bot'][lay])

    return model_ds

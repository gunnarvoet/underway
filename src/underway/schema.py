"""Common variable names, units, and dataset helpers."""

import numpy as np
import xarray as xr

# core name -> (long_name, units)
CORE = {
    "lon": ("longitude", "degrees_east"),
    "lat": ("latitude", "degrees_north"),
    "heading": ("ship heading", "deg"),
    "cog": ("course over ground", "deg"),
    "sog": ("speed over ground", "m/s"),
    "sst": ("sea surface temperature", "°C"),
    "sss": ("sea surface salinity", "psu"),
    "wind_speed": ("true wind speed", "m/s"),
    "wind_direction": ("true wind direction (coming from)", "deg"),
    "air_temperature": ("air temperature", "°C"),
    "air_pressure": ("air pressure", "hPa"),
    "relative_humidity": ("relative humidity", "%"),
}

ANGULAR = ("heading", "cog", "wind_direction")

KNOTS_TO_MS = 1852 / 3600


def conform(ds, extra=None):
    """Set `long_name` and `units` on core and extra variables.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with variables already renamed and converted to schema units.
    extra : dict, optional
        Ship-specific variables as ``{name: (long_name, units)}``. Names not
        present in `ds` are ignored.

    Returns
    -------
    xr.Dataset
    """
    table = dict(extra or {})
    table.update(CORE)
    for name, (long_name, units) in table.items():
        if name in ds.variables:
            ds[name].attrs = {"long_name": long_name, "units": units}
    return ds


def empty():
    """Return a dataset with a zero-length time coordinate."""
    return xr.Dataset(coords={"time": np.array([], dtype="datetime64[ns]")})


def combine(datasets):
    """Concatenate along time, sort, and keep each time stamp once.

    Parameters
    ----------
    datasets : list of xr.Dataset
        Datasets with a `time` dimension. Empty ones are skipped.

    Returns
    -------
    xr.Dataset
        Sorted by time, without `NaT` and without duplicate time stamps.
    """
    datasets = [ds for ds in datasets if ds.sizes.get("time", 0) > 0]
    if not datasets:
        return empty()
    ds = xr.concat(
        datasets,
        dim="time",
        data_vars="all",
        coords="minimal",
        compat="override",
        join="outer",
    )
    ds = ds.isel(time=~np.isnat(ds.time.values))
    _, index = np.unique(ds.time.values, return_index=True)
    return ds.isel(time=index)

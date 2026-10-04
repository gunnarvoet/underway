"""RRS Discovery TechSAS netCDF files (GPS, SURFMET, TSG)."""

import logging

import xarray as xr

from .. import schema
from ._common import as_paths

log = logging.getLogger("underway")

# TechSAS name -> schema or ship-specific name
NAMES = {
    # position-POSMV_GPS.gps
    "long": "lon",
    "gndcourse": "cog",
    "gndspeed": "sog",
    # MET-SURFMET.SURFMETv3 (apparent wind only)
    "airtemp": "air_temperature",
    "humid": "relative_humidity",
    "direct": "apparent_wind_direction",
    "speed": "apparent_wind_speed",
    # Light-SURFMET.SURFMETv3
    "pres": "air_pressure",
    # SBE45-SBE45.TSG
    "salin": "sss",
    "temp_r": "sst",
    "temp_h": "tsg_temperature",
    "cond": "conductivity",
    "sndspeed": "sound_speed",
}

EXTRA = {
    "apparent_wind_direction": ("apparent wind direction", "deg"),
    "apparent_wind_speed": ("apparent wind speed", "m/s"),
    "tsg_temperature": ("TSG housing temperature", "°C"),
    "conductivity": ("conductivity", "S/m"),
    "sound_speed": ("sound speed", "m/s"),
}


def _read_file(file):
    if file.stat().st_size == 0:
        log.warning("%s: zero bytes, skipped", file.name)
        return schema.empty()
    with xr.open_dataset(file, engine="netcdf4") as ds:
        ds = ds.load()
    ds = ds.reset_coords()
    ds = ds.drop_vars("measureTS", errors="ignore")
    ds = ds.rename({k: v for k, v in NAMES.items() if k in ds.variables})
    if "sog" in ds:
        ds["sog"] = ds["sog"] * schema.KNOTS_TO_MS
    return schema.conform(ds, extra=EXTRA)


def read(files):
    """Read TechSAS netCDF files of one kind.

    Parameters
    ----------
    files : path or iterable of paths
        Files of a single kind, e.g. all ``*position-POSMV_GPS.gps``.

    Returns
    -------
    xr.Dataset
        Core variables renamed and converted to schema units. Variables with
        no schema name keep the TechSAS name and attributes.

    Notes
    -----
    Wind in the SURFMET files is apparent wind. The schema names
    `wind_speed` and `wind_direction` are therefore absent.
    """
    return schema.combine([_read_file(f) for f in as_paths(files)])

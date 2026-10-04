"""R/V Roger Revelle MetAcq files (``<yymmdd>.MET``)."""

import logging

import pandas as pd

from .. import schema
from ._common import as_paths

log = logging.getLogger("underway")

MISSING = -99

# MetAcq tag -> schema or ship-specific name. Tags and units are listed in
# MetAcq.pdf, Appendix A. Tags not listed here are dropped.
NAMES = {
    "LA": "lat",
    "LO": "lon",
    "GY": "heading",
    "CR": "cog",
    "SP": "sog",
    "AT": "air_temperature",
    "BP": "air_pressure",
    "RH": "relative_humidity",
    "TW": "wind_speed",
    "TI": "wind_direction",
    "SA": "sss",
    "WS": "relative_wind_speed",
    "WD": "relative_wind_direction",
    "SW": "shortwave_radiation",
    "LW": "longwave_radiation",
    "PA": "par",
    "PR": "precipitation",
    "FL": "fluorometer",
    "MB": "water_depth_multibeam",
    "LF": "water_depth_3_5khz",
    "HF": "water_depth_12khz",
}

EXTRA = {
    "relative_wind_speed": ("relative wind speed", "m/s"),
    "relative_wind_direction": ("relative wind direction (coming from)", "deg"),
    "shortwave_radiation": ("shortwave radiation", "W/m^2"),
    "longwave_radiation": ("longwave radiation", "W/m^2"),
    "par": ("surface PAR", "uE/s/m^2"),
    "precipitation": ("precipitation", "mm"),
    "fluorometer": ("fluorometer", "ug/l"),
    "water_depth_multibeam": ("multibeam water depth", "m"),
    "water_depth_3_5khz": ("3.5 kHz water depth", "m"),
    "water_depth_12khz": ("12 kHz water depth", "m"),
}


def _read_start_date(file):
    """Return the date on header line 2, e.g. ``# Wed 12-Apr-17  02:12:38``."""
    with open(file) as f:
        f.readline()
        fields = f.readline().split()
    return pd.to_datetime(fields[2], format="%d-%b-%y")


def _read_file(file):
    if file.stat().st_size == 0:
        return schema.empty()
    date = _read_start_date(file)
    df = pd.read_csv(
        file, sep=r"\s+", skiprows=3, dtype={"#Time": str}, on_bad_lines="skip"
    )
    if df.empty:
        return schema.empty()
    hhmmss = df["#Time"]
    seconds = (
        pd.to_numeric(hhmmss.str[:2], errors="coerce") * 3600
        + pd.to_numeric(hhmmss.str[2:4], errors="coerce") * 60
        + pd.to_numeric(hhmmss.str[4:6], errors="coerce")
    )
    # a drop of more than 12 h in time of day marks the next UTC day
    day = (seconds.diff() < -43200).cumsum()
    time = date + pd.to_timedelta(day, unit="D") + pd.to_timedelta(seconds, unit="s")
    # sea surface temperature tag is ST where present, else the TSG tag TT
    sst_tag = "ST" if "ST" in df.columns else "TT"
    names = dict(NAMES)
    names[sst_tag] = "sst"
    df = df[[c for c in df.columns if c in names]]
    df = df.apply(pd.to_numeric, errors="coerce").rename(columns=names)
    df = df.where(df != MISSING)
    df["time"] = time
    dropped = int(df["time"].isna().sum())
    if dropped:
        log.warning("%s: dropped %d malformed lines", file.name, dropped)
    df = df.dropna(subset=["time"])
    if df.empty:
        return schema.empty()
    if "sog" in df:
        df["sog"] = df["sog"] * schema.KNOTS_TO_MS
    ds = df.set_index("time").to_xarray()
    return schema.conform(ds, extra=EXTRA)


def read(files):
    """Read Revelle MetAcq files.

    Parameters
    ----------
    files : path or iterable of paths
        Daily ``<yymmdd>.MET`` files.

    Returns
    -------
    xr.Dataset
        Tags listed in `NAMES` renamed and in schema units. `sst` comes from
        tag ``ST`` where the file has it and from ``TT`` otherwise. ``-99``
        is replaced by `NaN`.
    """
    return schema.combine([_read_file(f) for f in as_paths(files)])

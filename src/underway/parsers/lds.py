"""Sikuliaq LDS text streams (TSG, wind, air)."""

import logging
from dataclasses import dataclass, field

import pandas as pd

from .. import schema
from ._common import as_paths, read_lds_lines

log = logging.getLogger("underway")


@dataclass(frozen=True)
class Stream:
    """Column layout of one LDS stream.

    Parameters
    ----------
    columns : dict
        Index of the comma-separated payload field -> variable name.
    scale : dict
        Variable name -> factor that converts the logged value to schema units.
    extra : dict
        Ship-specific variables as ``{name: (long_name, units)}``.
    """

    columns: dict
    scale: dict = field(default_factory=dict)
    extra: dict = field(default_factory=dict)


STREAMS = {
    # " 26.9030,  5.51282,  35.0429, 1538.914"
    "tsg": Stream(
        columns={0: "sst", 1: "conductivity", 2: "sss", 3: "sound_speed"},
        extra={
            "conductivity": ("conductivity", "S/m"),
            "sound_speed": ("sound speed", "m/s"),
        },
    ),
    # "$WIMWD,72.6,T,,M,8.6,N,4.4,M*49"
    "wind": Stream(columns={1: "wind_direction", 7: "wind_speed"}),
    # "$WIXDR,PRESS,1.012158,bar,s/n146582,TEMP,27.41,C,RH,55.36,%RH,..."
    "air": Stream(
        columns={2: "air_pressure", 6: "air_temperature", 9: "relative_humidity"},
        scale={"air_pressure": 1000.0},  # bar -> hPa
    ),
}


def _read_file(file, stream):
    spec = STREAMS[stream]
    lines = read_lds_lines(file)
    if lines.empty:
        return schema.empty()
    time = pd.to_datetime(lines["time"], format="ISO8601", utc=True, errors="coerce")
    fields = lines["payload"].str.split("*", n=1).str[0].str.split(",", expand=True)
    out = pd.DataFrame({"time": time.dt.tz_convert(None)})
    for column, name in spec.columns.items():
        if column in fields.columns:
            values = pd.to_numeric(fields[column].str.strip(), errors="coerce")
        else:
            values = float("nan")
        out[name] = values * spec.scale.get(name, 1.0)
    good = out.dropna()
    dropped = len(out) - len(good)
    if dropped:
        log.warning("%s: dropped %d malformed lines", file.name, dropped)
    if good.empty:
        return schema.empty()
    ds = good.set_index("time").to_xarray()
    return schema.conform(ds, extra=spec.extra)


def read(files, stream):
    """Read Sikuliaq LDS text files of one stream.

    Parameters
    ----------
    files : path or iterable of paths
        Raw LDS files, e.g. ``tsg_sbe45_fwd.20251127T0000Z``.
    stream : {"tsg", "wind", "air"}
        Selects the column layout in `STREAMS`.

    Returns
    -------
    xr.Dataset
        Variables in schema units on a `time` dimension.

    Examples
    --------
    >>> ds = read(sorted(raw_dir.glob("tsg_sbe45_fwd.*")), stream="tsg")
    """
    if stream not in STREAMS:
        raise ValueError(f"unknown stream {stream!r}, choose from {sorted(STREAMS)}")
    return schema.combine([_read_file(f, stream) for f in as_paths(files)])

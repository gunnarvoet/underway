"""R/V Neil Armstrong dsLog met csv files and raw CNAV GPS files."""

import io
import logging

import pandas as pd

from .. import schema
from . import nmea
from ._common import as_paths, read_complete

log = logging.getLogger("underway")

# dsLog column -> schema or ship-specific name. Core names take the starboard
# sensor (WXTS), port values carry the suffix _port.
MET_NAMES = {
    "Dec_LAT": "lat",
    "Dec_LON": "lon",
    "SPD": "spd",
    "HDT": "heading",
    "COG": "cog",
    "SOG": "sog",
    "WXTS_Ta": "air_temperature",
    "WXTP_Ta": "air_temperature_port",
    "WXTS_Pa": "air_pressure",
    "WXTP_Pa": "air_pressure_port",
    "WXTS_Ri": "rain_intensity",
    "WXTP_Ri": "rain_intensity_port",
    "WXTS_Rc": "rain_accumulation",
    "WXTP_Rc": "rain_accumulation_port",
    "WXTS_Dm": "relative_wind_direction",
    "WXTP_Dm": "relative_wind_direction_port",
    "WXTS_Sm": "relative_wind_speed",
    "WXTP_Sm": "relative_wind_speed_port",
    "WXTS_Ua": "relative_humidity",
    "WXTP_Ua": "relative_humidity_port",
    "WXTS_TS": "wind_speed",
    "WXTP_TS": "wind_speed_port",
    "WXTS_TD": "wind_direction",
    "WXTP_TD": "wind_direction_port",
    "BAROM_S": "barometric_pressure",
    "BAROM_P": "barometric_pressure_port",
    "RAD_SW": "shortwave_radiation",
    "RAD_LW": "longwave_radiation",
    "PAR": "par",
    "SBE45S": "sss",
    "SBE48T": "sst",
    "FLR": "fluorometer",
    "FLOW": "flow",
    "SSVdslog": "sound_speed",
    "Depth12": "water_depth_12khz",
    "Depth35": "water_depth_3_5khz",
    "EM122": "water_depth_em122",
    "EM710": "water_depth_em710",
}

# The csv files carry no units. SOG and SPD in knots and wind speed in m/s
# were checked against positions in AR73 data. The other units follow the
# sensor types and are unverified.
MET_EXTRA = {
    "spd": ("ship speed", "knots"),
    "air_temperature_port": ("air temperature port", "°C"),
    "air_pressure_port": ("air pressure port", "hPa"),
    "rain_intensity": ("rain intensity", "mm/h"),
    "rain_intensity_port": ("rain intensity port", "mm/h"),
    "rain_accumulation": ("rain accumulation", "mm"),
    "rain_accumulation_port": ("rain accumulation port", "mm"),
    "relative_wind_direction": ("relative wind direction", "deg"),
    "relative_wind_direction_port": ("relative wind direction port", "deg"),
    "relative_wind_speed": ("relative wind speed", "m/s"),
    "relative_wind_speed_port": ("relative wind speed port", "m/s"),
    "relative_humidity_port": ("relative humidity port", "%"),
    "wind_speed_port": ("true wind speed port", "m/s"),
    "wind_direction_port": ("true wind direction port (coming from)", "deg"),
    "barometric_pressure": ("barometric pressure", "hPa"),
    "barometric_pressure_port": ("barometric pressure port", "hPa"),
    "shortwave_radiation": ("shortwave radiation", "W/m^2"),
    "longwave_radiation": ("longwave radiation", "W/m^2"),
    "par": ("photosynthetically active radiation", "uE/m^2/s"),
    "fluorometer": ("fluorometer", "mV"),
    "flow": ("flow", "unknown"),
    "sound_speed": ("sea surface sound speed", "m/s"),
    "water_depth_12khz": ("12 kHz water depth", "m"),
    "water_depth_3_5khz": ("3.5 kHz water depth", "m"),
    "water_depth_em122": ("EM122 multibeam water depth", "m"),
    "water_depth_em710": ("EM710 multibeam water depth", "m"),
}


def _read_met_file(file):
    data, partial = read_complete(file)
    try:
        df = pd.read_csv(
            io.BytesIO(data),
            skiprows=1,
            skipinitialspace=True,
            na_values=["NAN", "NODATA"],
            on_bad_lines="skip",
        )
    except pd.errors.EmptyDataError:
        return schema.empty()
    if df.empty:
        return schema.empty()
    time = pd.to_datetime(
        df["DATE_GMT"] + " " + df["TIME_GMT"],
        format="%Y/%m/%d %H:%M:%S.%f",
        errors="coerce",
    )
    df = df.drop(columns=["DATE_GMT", "TIME_GMT"])
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.rename(columns=MET_NAMES)
    df["time"] = time
    dropped = int(df["time"].isna().sum()) + partial
    if dropped:
        log.warning("%s: dropped %d malformed lines", file.name, dropped)
    df = df.dropna(subset=["time"])
    if df.empty:
        return schema.empty()
    df["sog"] = df["sog"] * schema.KNOTS_TO_MS
    ds = df.set_index("time").to_xarray()
    return schema.conform(ds, extra=MET_EXTRA)


def read_met(files):
    """Read Armstrong dsLog met files (``AR<yymmdd>_0000.csv``).

    Parameters
    ----------
    files : path or iterable of paths

    Returns
    -------
    xr.Dataset
        One-minute met and navigation record. Core names hold the starboard
        sensor. Port sensor values carry the suffix ``_port``.
    """
    return schema.combine([_read_met_file(f) for f in as_paths(files)])


def _read_gps_file(file):
    data, partial = read_complete(file)
    lines = pd.Series(data.decode(errors="replace").splitlines(), dtype=str)
    sentences = lines.str.split(" CNAV ", n=1).str[1]
    return nmea.parse(sentences, name=file.name, dropped=partial)


def read_gps(files):
    """Read Armstrong raw CNAV GPS files (``*.CNAV_3050``).

    Parameters
    ----------
    files : path or iterable of paths

    Returns
    -------
    xr.Dataset
        `lat`, `lon`, `cog`, `sog` at 1 Hz.
    """
    return schema.combine([_read_gps_file(f) for f in as_paths(files)])

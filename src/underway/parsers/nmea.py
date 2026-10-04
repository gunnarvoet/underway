"""Vectorized parsing of NMEA 0183 sentences into one row per position fix."""

import logging

import numpy as np
import pandas as pd

from .. import schema

log = logging.getLogger("underway")

EXTRA = {
    "roll": ("roll", "deg"),
    "pitch": ("pitch", "deg"),
    "heave": ("heave", "m"),
}


def _num(series):
    return pd.to_numeric(series, errors="coerce")


def _seconds_of_day(hhmmss):
    """Convert ``hhmmss[.ss]`` strings to seconds since midnight.

    Values outside a valid time of day give `NaN`.
    """
    hours = _num(hhmmss.str[:2])
    minutes = _num(hhmmss.str[2:4])
    seconds = _num(hhmmss.str[4:])
    valid = (hours < 24) & (minutes < 60) & (seconds < 61)
    return (hours * 3600 + minutes * 60 + seconds).where(valid)


def _degrees(value, hemisphere, ndeg):
    """Convert ``d..dmm.mmmm`` strings with `ndeg` degree digits to degrees."""
    degrees = _num(value.str[:ndeg]) + _num(value.str[ndeg:]) / 60
    sign = np.where(hemisphere.isin(["S", "W"]), -1.0, 1.0)
    return degrees * sign


def parse(sentences, name="", dropped=0):
    """Parse NMEA sentences into a dataset with one time step per ZDA sentence.

    Sentences are assigned to the most recent ZDA sentence, so each group of
    sentences must start with ZDA. A sentence missing from one fix leaves
    `NaN` in that fix only. Sentences ahead of the first ZDA and fixes whose
    time cannot be parsed are dropped and counted. The talker id (``GP``,
    ``IN``, ...) is ignored. A GGA position is used only if its time equals
    the ZDA time of its group. If GGA sentences are present and none of them
    match, the sentence order is not ZDA-first and a warning is logged.

    Parameters
    ----------
    sentences : pd.Series of str
        One NMEA sentence per element, in file order, checksum optional.
    name : str, optional
        File name used in the log message.
    dropped : int, optional
        Number of lines the caller already dropped, added to the logged count.

    Returns
    -------
    xr.Dataset
        `lat`, `lon` from GGA, `cog`, `sog` (m/s) from VTG, `heading` from
        HDT, `roll`, `pitch`, `heave` from ``$PSXN,23``. Variables without
        any value are left out.
    """
    s = sentences.dropna().astype(str).str.strip().str.split("*", n=1).str[0]
    s = s.reset_index(drop=True)
    if s.empty:
        return schema.empty()
    f = s.str.split(",", expand=True).reindex(columns=range(10)).astype(object)
    tag = f[0]
    # sentence type without the two-letter talker id, e.g. "$GPZDA" -> "ZDA"
    kind = tag.str[3:].where(tag.str.len() == 6)
    fix = (kind == "ZDA").cumsum()
    dropped += int((fix == 0).sum())

    def pick(mask, columns):
        rows = f.loc[mask & (fix > 0), columns]
        rows.index = fix[rows.index]
        return rows[~rows.index.duplicated()]

    zda = pick(kind == "ZDA", [1, 2, 3, 4])
    date = pd.to_datetime(
        zda[4] + "-" + zda[3] + "-" + zda[2], format="%Y-%m-%d", errors="coerce"
    )
    seconds = _seconds_of_day(zda[1])
    out = pd.DataFrame(index=zda.index)
    out["time"] = date + pd.to_timedelta(seconds, unit="s").dt.round("ms")

    gga = pick(kind == "GGA", [1, 2, 3, 4, 5])
    n_gga = len(gga)
    gga = gga[_num(gga[1]) == _num(zda[1].reindex(gga.index))]
    if n_gga and gga.empty:
        log.warning(
            "%s: no GGA time matches its ZDA time, sentence order not supported",
            name,
        )
    out["lat"] = _degrees(gga[2], gga[3], 2)
    out["lon"] = _degrees(gga[4], gga[5], 3)

    vtg = pick(kind == "VTG", [1, 5])
    out["cog"] = _num(vtg[1])
    out["sog"] = _num(vtg[5]) * schema.KNOTS_TO_MS

    hdt = pick(kind == "HDT", [1])
    out["heading"] = _num(hdt[1])

    psxn = pick((tag == "$PSXN") & (f[1] == "23"), [2, 3, 5])
    out["roll"] = _num(psxn[2])
    out["pitch"] = _num(psxn[3])
    out["heave"] = _num(psxn[5])

    bad_time = int(out["time"].isna().sum())
    dropped += bad_time
    if dropped:
        log.warning("%s: dropped %d malformed lines", name, dropped)
    out = out.dropna(subset=["time"]).dropna(axis=1, how="all")
    if out.empty:
        return schema.empty()
    ds = out.set_index("time").to_xarray()
    return schema.conform(ds, extra=EXTRA)

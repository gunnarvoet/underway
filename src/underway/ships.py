"""Ship servers, drives, and default sources. Data only.

Remote paths are carried over from at-sea use and change between cruises.
Override an entry for one cruise with `Cruise.add_source` under the same name.
In `read_met`, the first source that provides a variable wins, so GPS comes
first.

## Adding a ship

A new ship needs a parser for each raw file format that no existing parser
reads, an entry in `SHIPS`, and tests against real files.

### 1. Check the existing parsers

A format the ship shares with another ship needs no new code.

- `underway.parsers.nmea.parse` reads NMEA sentences (ZDA, GGA, VTG, HDT,
  ``$PSXN,23``) once each line is reduced to the sentence. Each group of
  sentences must start with ZDA. `underway.parsers.seapath` and
  `underway.parsers.armstrong.read_gps` are two wrappers around it.
- `underway.parsers.techsas` reads TechSAS netCDF files.
- `underway.parsers.lds` reads tab-separated text streams. A new stream is
  one more entry in `underway.parsers.lds.STREAMS`.

### 2. Write a parser

A parser is a function ``read(files) -> xr.Dataset`` in a new module under
``underway/parsers/``. It knows nothing about `underway.cruise.Cruise`, the
cache, or the directory layout, so it also works on shore on any list of
files.

- Accept one path or a list of paths with ``_common.as_paths``.
- Read text files with ``_common.read_complete``. It cuts off a last line
  that does not end in a newline, since the file may still be written.
- Rename columns to the core names in `underway.schema.CORE` and convert to
  the schema units (speeds in m/s, pressure in hPa, true wind). Take the
  units from the raw file header or the ship's documentation.
- Variables without a core name keep a ship-specific name. Pass their
  ``(long_name, units)`` to `underway.schema.conform`.
- Drop malformed lines, count them, and log the count on logger
  ``underway``. Return `underway.schema.empty` for a file without data.
- Join the per-file datasets with `underway.schema.combine`.

A parser for a csv file with the columns ``TIME, LAT, LON, SOG_KN, SST, PAR``:

```python
# src/underway/parsers/example.py
import io
import logging

import pandas as pd

from .. import schema
from ._common import as_paths, read_complete

log = logging.getLogger("underway")

NAMES = {"LAT": "lat", "LON": "lon", "SOG_KN": "sog", "SST": "sst", "PAR": "par"}
EXTRA = {"par": ("photosynthetically active radiation", "uE/m^2/s")}


def _read_file(file):
    data, partial = read_complete(file)
    try:
        df = pd.read_csv(io.BytesIO(data))
    except pd.errors.EmptyDataError:
        return schema.empty()
    time = pd.to_datetime(df["TIME"], errors="coerce")
    df = df[list(NAMES)].apply(pd.to_numeric, errors="coerce").rename(columns=NAMES)
    df["sog"] = df["sog"] * schema.KNOTS_TO_MS
    df["time"] = time
    dropped = int(df["time"].isna().sum()) + partial
    if dropped:
        log.warning("%s: dropped %d malformed lines", file.name, dropped)
    df = df.dropna(subset=["time"])
    if df.empty:
        return schema.empty()
    return schema.conform(df.set_index("time").to_xarray(), extra=EXTRA)


def read(files):
    return schema.combine([_read_file(f) for f in as_paths(files)])
```

Add the module to the imports and ``__all__`` in
``underway/parsers/__init__.py``.

### 3. Add the ship to `SHIPS`

```python
"example": Ship(
    name="R/V Example",
    servers={"data": "files.example.edu"},
    sources=(
        Source(
            "met",
            drive="data",
            remote="{cruise_id}/met",
            pattern="*.csv",
            parser=example.read,
            met=True,
        ),
        Source(
            "sadcp",
            drive="data",
            remote="{cruise_id}/adcp/proc",
            pattern="*/contour/*.nc",
        ),
        Source("ctd", drive="data", remote="{cruise_id}/ctd", readonly=True),
    ),
),
```

- ``servers`` maps each share (drive) to its server. A drive appears under
  the mount root, ``/Volumes/<drive>`` on macOS.
- ``remote`` is the directory below the drive. ``{cruise_id}`` is filled in
  from the `underway.cruise.Cruise`.
- ``pattern`` is a glob relative to ``remote`` and may contain directories.
- A source without a parser is synced only.
- Use ``transfer="copy"`` where rsync fails on the share, ``cache=False``
  for raw files that are netCDF already, and ``met=True`` for sources that
  feed `underway.cruise.Cruise.read_met`. List the GPS source first.

All fields are described in `underway.source.Source`.

### 4. Test against real files

- Add a truncated real file to ``tests/make_fixtures.py`` and rebuild
  ``tests/data``.
- Test the first row against values read from the file by eye, a unit
  conversion with a nonzero value, a last line cut inside a field, and a
  zero-byte file. ``tests/test_revelle.py`` is a compact model.
- ``tests/test_ships.py`` checks every entry in `SHIPS` for consistency and
  picks up the new ship without changes.
- Add a check over the full files to ``tests/test_real_data.py``.

### 5. Use it

```python
import underway as uw

c = uw.Cruise("example", "EX2601", "~/data/ex2601")
c.sync("met")
met = c.read("met")
```

A single extra stream on a ship that is already supported needs no change to
the package. Declare it for the cruise with
`underway.cruise.Cruise.add_source`.
"""

from dataclasses import dataclass
from functools import partial

from .parsers import armstrong, lds, revelle, seapath, techsas
from .source import Source


@dataclass(frozen=True)
class Ship:
    """A research vessel.

    Parameters
    ----------
    name : str
        Display name.
    servers : dict
        Drive name -> server name.
    sources : tuple of Source
        Default sources.
    """

    name: str
    servers: dict[str, str]
    sources: tuple[Source, ...]


_LDS = "{cruise_id}/lds/raw"
_TECHSAS = "Ship_Systems/Data/TechSAS/NetCDF"

SHIPS = {
    "sikuliaq": Ship(
        name="R/V Sikuliaq",
        servers={
            "CruiseData": "data.sikuliaq.alaska.edu",
            "science": "files.sikuliaq.alaska.edu",
        },
        sources=(
            Source(
                "gps",
                drive="CruiseData",
                remote=f"{_LDS}/ins_seapath_position",
                pattern="ins_seapath_position.*",
                parser=seapath.read,
                met=True,
            ),
            Source(
                "tsg",
                drive="CruiseData",
                remote=f"{_LDS}/tsg_sbe45_fwd",
                pattern="tsg_sbe45_fwd.*",
                parser=partial(lds.read, stream="tsg"),
                met=True,
            ),
            Source(
                "wind",
                drive="CruiseData",
                remote=f"{_LDS}/wind_gill_fwdmast_true",
                pattern="wind_gill_fwdmast_true.*",
                parser=partial(lds.read, stream="wind"),
                met=True,
            ),
            Source(
                "air",
                drive="CruiseData",
                remote=f"{_LDS}/met_met4a_fwdmast",
                pattern="met_met4a_fwdmast.*",
                parser=partial(lds.read, stream="air"),
                met=True,
            ),
            Source(
                "sadcp",
                drive="CruiseData",
                remote="{cruise_id}/adcp/raw/{cruise_id}/proc",
                pattern="*/contour/*.nc",
                transfer="copy",
            ),
            Source(
                "ctd",
                drive="CruiseData",
                remote="{cruise_id}/ctd/raw",
                readonly=True,
            ),
        ),
    ),
    # TechSAS and ADCP files use transfer="copy". rsync was replaced by a
    # size-compare copy for these in June 2021 (commit cd20c84), files arrive
    # locked (hence chflags in transfer.copy), and shutil.copy2 gave
    # permission errors (hence copyfile). CTD and LADCP stayed on rsync.
    "discovery": Ship(
        name="RRS Discovery",
        servers={
            "current_cruise": "dynetapp.discovery.ad.noc.ac.uk",
            "science_public": "dynetapp.discovery.ad.noc.ac.uk",
        },
        sources=(
            Source(
                "gps",
                drive="current_cruise",
                remote=f"{_TECHSAS}/GPS",
                pattern="*position-POSMV_GPS.gps",
                transfer="copy",
                parser=techsas.read,
                cache=False,
                met=True,
            ),
            Source(
                "surfmet",
                drive="current_cruise",
                remote=f"{_TECHSAS}/SURFMETV3",
                pattern="*MET-SURFMET.SURFMETv3",
                transfer="copy",
                parser=techsas.read,
                cache=False,
                met=True,
            ),
            Source(
                "light",
                drive="current_cruise",
                remote=f"{_TECHSAS}/SURFMETV3",
                pattern="*Light-SURFMET.SURFMETv3",
                transfer="copy",
                parser=techsas.read,
                cache=False,
                met=True,
            ),
            Source(
                "surf",
                drive="current_cruise",
                remote=f"{_TECHSAS}/SURFMETV3",
                pattern="*Surf-SURFMET.SURFMETv3",
                transfer="copy",
                parser=techsas.read,
                cache=False,
            ),
            Source(
                "tsg",
                drive="current_cruise",
                remote=f"{_TECHSAS}/TSG",
                pattern="*SBE45-SBE45.TSG",
                transfer="copy",
                parser=techsas.read,
                cache=False,
                met=True,
            ),
            Source(
                "sadcp",
                drive="current_cruise",
                remote="Ship_Systems/Data/Acoustics/ADCP/proc",
                pattern="os*nb/contour/*",
                transfer="copy",
            ),
            Source(
                "ctd",
                drive="current_cruise",
                remote="Sensors_and_Moorings/CTD/Data/Raw",
            ),
            Source(
                "ladcp",
                drive="current_cruise",
                remote="Sensors_and_Moorings/LADCP/Data",
            ),
        ),
    ),
    "armstrong": Ship(
        name="R/V Neil Armstrong",
        servers={
            "data_on_memory": "10.100.100.30",
            "science_share": "10.100.100.30",
        },
        sources=(
            Source(
                "gps",
                drive="data_on_memory",
                remote="underway/raw",
                pattern="*.CNAV_3050",
                transfer="copy",
                parser=armstrong.read_gps,
                met=True,
            ),
            Source(
                "met",
                drive="data_on_memory",
                remote="underway/proc",
                pattern="AR[0-9]*.csv",
                parser=armstrong.read_met,
                met=True,
            ),
            Source(
                "sadcp",
                drive="data_on_memory",
                remote="adcp/proc",
                pattern="*/contour/*.nc",
            ),
            Source("ctd", drive="data_on_memory", remote="ctd", readonly=True),
        ),
    ),
    "revelle": Ship(
        name="R/V Roger Revelle",
        servers={
            "cruise": "rr-sci-filesvr.ucsd.edu",
            "science_party_share": "rr-sci-filesvr.ucsd.edu",
        },
        sources=(
            Source(
                "met",
                drive="cruise",
                remote="{cruise_id}/metacq/data",
                pattern="*.MET",
                parser=revelle.read,
                met=True,
            ),
            Source(
                "sadcp",
                drive="cruise",
                remote="{cruise_id}/adcp_uhdas/{cruise_id}/proc",
                pattern="*/contour/*.nc",
            ),
        ),
    ),
}

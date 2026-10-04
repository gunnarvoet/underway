"""Research vessel underway data handling: sync shipboard data from the ship
server and parse it into `xarray.Dataset`s. Everything here is under
development and may change at any time.

## Features

- Sync shipboard data (GPS, met, TSG, SADCP, CTD, LADCP) from the ship server
  with `rsync` or a size-compare copy. Mounting the server is automatic on
  macOS. On Linux, mount the shares by hand and pass `mount_root`.
- Parse raw files into `xarray.Dataset` with common variable names and units
  across ships.
- Cache one netCDF product per raw file. A raw file is parsed again only when
  its size changed.
- Parsers are plain functions over file lists and work on any directory
  layout.

## Ships

R/V Sikuliaq, RRS Discovery, R/V Neil Armstrong, R/V Roger Revelle. Servers
and default sources per ship are listed in `underway.ships`.

## Example

```python
import underway as uw

c = uw.Cruise("sikuliaq", "SKQ202521S", "/path/to/your/local/cruise/dir/")
# add a cruise-specific source
c.add_source("ww", drive="science", remote="ww/nc", pattern="*.nc")
# sync data to local computer
c.sync()
# parse what changed and read the full record
gps = c.read("gps")
# core variables of all met sources on a one-minute grid
met = c.read_met()
# local directories
c.path("ctd", "raw")
```

On shore, with raw files in any layout:

```python
gps = uw.parsers.seapath.read(sorted(raw_dir.glob("ins_seapath_position.*")))
```

## Variable names

Core variables are present where the ship provides them: `lon`, `lat`,
`heading`, `cog`, `sog`, `sst`, `sss`, `wind_speed`, `wind_direction`,
`air_temperature`, `air_pressure`, `relative_humidity`. Speeds are in m/s,
pressure in hPa, wind is true wind. All other variables keep ship-specific
names. Names and units are defined in `underway.schema`.

## Adding a ship

A new ship needs a parser for each raw file format that no existing parser
reads, an entry in the ship table, and tests against real files. The steps,
with a worked example, are in `underway.ships`. A single extra stream on a
supported ship needs no change to the package. Declare it for the cruise with
`underway.cruise.Cruise.add_source`.

## Modules

- `underway.cruise`: the `Cruise` object with `sync`, `read`, `read_met`.
- `underway.source`, `underway.ships`: declared sources and the ship tables.
- `underway.parsers`: one module per raw file format.
- `underway.schema`: core variable names, units, and dataset helpers.
- `underway.cache`, `underway.transfer`, `underway.mount`: the netCDF cache,
  file transfer, and mounting of the ship's shares.
"""

import importlib.metadata

from . import cache, mount, parsers, schema, ships, transfer
from .cruise import Cruise
from .source import Source

__all__ = [
    "Cruise",
    "Source",
    "cache",
    "cruise",
    "mount",
    "parsers",
    "schema",
    "ships",
    "source",
    "transfer",
]

__author__ = "Gunnar Voet"
__email__ = "gvoet@ucsd.edu"
# version is defined in pyproject.toml
__version__ = importlib.metadata.version("underway")

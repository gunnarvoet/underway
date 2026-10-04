# underway

Research vessel underway data handling. Everything here is under development and may change at any time.

# Features

* Sync shipboard data (GPS, met, TSG, SADCP, CTD, LADCP) from the ship server with `rsync` or a size-compare copy. Mounting the server is automatic on macOS (`osascript`). On Linux, mount the shares by hand and pass `mount_root`.
* Parse raw files into `xarray.Dataset` with common variable names and units across ships.
* Cache one netCDF product per raw file. A raw file is parsed again only when its size changed.
* Parsers are plain functions over file lists and work on any directory layout.

# Ships

* R/V Sikuliaq
* RRS Discovery
* R/V Neil Armstrong
* R/V Roger Revelle

# Example

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

# Variable names

Core variables present where the ship provides them: `lon`, `lat`, `heading`, `cog`, `sog`, `sst`, `sss`, `wind_speed`, `wind_direction`, `air_temperature`, `air_pressure`, `relative_humidity`. Speeds are in m/s, pressure in hPa, wind is true wind. All other variables keep ship-specific names.

# Documentation

The documentation is at <https://gunnarvoet.github.io/underway/underway.html>. It includes a guide to adding a new ship on the page of the `underway.ships` module.

It is built with [pdoc](https://pdoc.dev/) from the docstrings:

```sh
git submodule update --init   # theme, once after cloning
make docs                     # build into docs/ and open
make servedocs                # live preview
```

# Old interface

The interface up to version 2025.10 (`uw.ship.Sikuliaq` and others) is available at git tag `v2025.10`.

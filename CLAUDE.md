# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```
uv sync                                                    # create env, install package + dev group
uv run pytest                                              # all tests
uv run pytest tests/test_lds.py::test_read_tsg_first_sample   # single test
uv run ruff check src tests                                # lint
uv run ruff format src tests                               # format
uv run python tests/make_fixtures.py                       # rebuild tests/data from the real files
```

- pytest runs with `filterwarnings = ["error"]` and one targeted ignore for the numpy warning that `netCDF4` emits on import. Any other warning fails a test. Pass `data_vars`, `coords`, `compat`, `join` explicitly to `xr.concat` and `xr.merge`.
- `tests/test_real_data.py` reads the full sample files (motive cruise 2, DY202, AR73 and FLEAT on the MOD archive) and skips where they are absent. `tests/make_fixtures.py` needs the same files.
- `tests/data` holds truncated real raw files. `.gitattributes` keeps their line endings byte-exact, since the Armstrong and Revelle files use CRLF.

## Architecture

`underway` syncs research vessel underway data from a ship's file server to a local cruise directory and parses it into `xarray.Dataset`s. The layers, bottom up:

- `parsers/`: one module per raw file format, each with `read(files) -> xr.Dataset`. Parsers know nothing about `Cruise`, the cache, or the directory layout, and must not import the layers above them. `nmea.parse` is shared by the Seapath and Armstrong CNAV readers. `lds.STREAMS` holds the column layout of each Sikuliaq text stream as data.
- `schema.py`: core variable names and units (`CORE`), `conform` to set attributes, `combine` to concatenate files (sorted, unique time stamps, no `NaT`), `bin_average` for the merged met product.
- `cache.py`: one netCDF product per raw file. The product stores the raw file's size. A raw file is parsed again in full only when its size differs. Products are written to a temporary file and renamed.
- `transfer.py`: `rsync` and a size-compare `copy`, both taking a glob relative to the source directory. Failures raise `TransferError`.
- `mount.py`: mounts SMB shares with `osascript` on macOS. On other platforms it reports drives that are missing under the mount root.
- `source.py`, `ships.py`: a `Source` describes one stream of files (drive, remote directory, pattern, transfer method, parser). `SHIPS` is a table of servers and default sources per ship and contains no logic.
- `cruise.py`: `Cruise(ship, cruise_id, local_dir, mount_root)` resolves paths and runs `sync`, `read`, and `read_met` over the source table. `add_source` adds a cruise-specific source or replaces a default one.

### Local layout

`<local_dir>/<source>/raw` for synced files and `<local_dir>/<source>/proc` for products. `Source.local` sets the raw directory exactly, so several sources can share one directory when their patterns differ. Nothing is created until `sync`, `read`, or `path(..., create=True)` needs it.

### Reading

`Cruise.read(name)` updates the product of every raw file whose size changed, then combines all products in `proc`. Products are read even when their raw file is absent (dangling git-annex link). After a parser change, call `read(name, reparse=True)`, since a parser change does not invalidate products on its own. Sources with `cache=False` (TechSAS netCDF) are parsed directly on every call.

`Cruise.read_met(freq)` bin-averages the core variables of every source with `met=True`. Angles are vector-averaged. Where two sources provide the same variable, the source listed first in `ships.py` is used.

## Adding a ship or a file format

1. Add a parser module in `parsers/` with `read(files)`. Map columns onto the core names in `schema.CORE`, convert to schema units, and pass ship-specific names with units to `schema.conform`.
2. Add a truncated real file to `tests/make_fixtures.py`, rebuild `tests/data`, and test the parser against known values from that file.
3. Add a `Ship` entry with its servers and sources to `ships.py`.

Take units from the raw file header or the ship's documentation. Where a unit cannot be verified, say so in a comment, as in `parsers/armstrong.py`.

## Platform constraints

- Mounting works on macOS only. On Linux, mount the shares by hand and pass `mount_root`.
- Server names, drive names, and remote paths in `ships.py` are carried over from at-sea use and cannot be checked ashore. They change between cruises.
- Discovery TechSAS and ADCP files use `transfer="copy"`. See the comment in `ships.py`.

## Conventions

- Date versioning (`YYYY.MM`), set only in `pyproject.toml`. Changes go in `HISTORY.md`.
- Numpy-style docstrings, `pathlib` over `os`, ruff for lint and format.
- Do not access a dataset variable named `roll` as an attribute. `ds.roll` is the xarray method. Use `ds["roll"]`.
- The interface up to 2025.10 is at git tag `v2025.10`.

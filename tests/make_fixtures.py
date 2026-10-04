"""Build tests/data from real raw files.

Run once on a machine that has the source data (elk):

    uv run python tests/make_fixtures.py

The outputs are committed, so the test suite does not need the sources.
"""

from pathlib import Path

import xarray as xr

HOME = Path.home()
OUT = Path(__file__).parent / "data"
MOTIVE = HOME / "Projects/motive/motive/data/cruise2"
DY202 = HOME / "Projects/madeira/data/cruises/dy202/met/raw"
AR73 = Path("/mnt/mod-server/TFO_NESMA/Cruises/ar73/04_ship_data/underway")
FLEAT = Path("/mnt/mod-server/FLEAT/Cruises")

# (source file, output subdirectory, number of lines to keep)
TEXT = [
    (MOTIVE / "gps/raw/ins_seapath_position.20251126T0000Z", "sikuliaq", 88),
    (MOTIVE / "tsg/raw/tsg_sbe45_fwd.20251127T0000Z", "sikuliaq", 38),
    (MOTIVE / "wind/raw/wind_gill_fwdmast_true.20251125T0000Z", "sikuliaq", 45),
    (MOTIVE / "air/raw/met_met4a_fwdmast.20251125T0000Z", "sikuliaq", 40),
    (AR73 / "proc/AR230512_0000.csv", "armstrong", 12),
    (AR73 / "raw/ar20230511_1600.CNAV_3050", "armstrong", 40),
    (FLEAT / "RR1708/data/met/data/170412.MET", "revelle", 14),
    (FLEAT / "RR1607/data/met/data/160530.MET", "revelle", 14),
]

# (file name in DY202, number of time steps to keep)
NETCDF = [
    ("20251024-145312-MET-SURFMET.SURFMETv3", 20),
    ("20251024-145312-Light-SURFMET.SURFMETv3", 20),
    ("20251024-145312-Surf-SURFMET.SURFMETv3", 20),
    ("20251024-145312-position-POSMV_GPS.gps", 20),
    ("20251025-000000-SBE45-SBE45.TSG", 20),
]


def main():
    for src, subdir, nlines in TEXT:
        out = OUT / subdir / src.name
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(src, "rb") as f:
            lines = [f.readline() for _ in range(nlines)]
        out.write_bytes(b"".join(lines))
        print(out)
    for name, nt in NETCDF:
        out = OUT / "discovery" / name
        out.parent.mkdir(parents=True, exist_ok=True)
        with xr.open_dataset(DY202 / name, engine="netcdf4", decode_times=False) as ds:
            ds.isel(time=slice(0, nt)).to_netcdf(out)
        print(out)


if __name__ == "__main__":
    main()

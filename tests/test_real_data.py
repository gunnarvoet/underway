"""Checks against full sample files. Skipped where the files are absent."""

import time
from pathlib import Path

import numpy as np
import pytest

from underway.parsers import armstrong, lds, revelle, seapath, techsas

HOME = Path.home()
MOTIVE = HOME / "Projects/motive/motive/data/cruise2"
DY202 = HOME / "Projects/madeira/data/cruises/dy202/met/raw"
AR73 = Path("/mnt/mod-server/TFO_NESMA/Cruises/ar73/04_ship_data/underway")
FLEAT = Path("/mnt/mod-server/FLEAT/Cruises")

GPS_DAY = MOTIVE / "gps/raw/ins_seapath_position.20251126T0000Z"


def _present(files):
    return sorted(f for f in files if f.is_file())


def _require(files):
    files = _present(files)
    if not files:
        pytest.skip("sample files not present")
    return files


def test_seapath_full_day_fix_count_and_speed():
    file = _require([GPS_DAY])[0]
    start = time.perf_counter()
    ds = seapath.read(file)
    elapsed = time.perf_counter() - start
    assert ds.sizes["time"] == 86400
    assert elapsed < 5.0


def test_seapath_full_day_missing_rmc_leaves_positions_complete():
    file = _require([GPS_DAY])[0]
    ds = seapath.read(file)
    assert int(np.isnan(ds.lat).sum()) == 0


@pytest.mark.parametrize(
    "stream, prefix",
    [
        ("tsg", "tsg_sbe45_fwd"),
        ("wind", "wind_gill_fwdmast_true"),
        ("air", "met_met4a_fwdmast"),
    ],
)
def test_lds_all_files_parse_with_monotonic_time(stream, prefix):
    files = _require((MOTIVE / stream / "raw").glob(f"{prefix}.*"))
    ds = lds.read(files, stream=stream)
    assert ds.sizes["time"] > 1000
    assert (np.diff(ds.time.values) > np.timedelta64(0, "ns")).all()


@pytest.mark.parametrize(
    "pattern",
    ["*position-POSMV_GPS.gps", "*MET-SURFMET.SURFMETv3", "*Light-SURFMET.SURFMETv3"],
)
def test_techsas_all_files_parse_with_monotonic_time(pattern):
    files = _require(DY202.glob(pattern))
    ds = techsas.read(files)
    assert ds.sizes["time"] > 1000
    assert (np.diff(ds.time.values) > np.timedelta64(0, "ns")).all()


def test_techsas_all_tsg_files_empty_in_dy202():
    files = _require(DY202.glob("*SBE45-SBE45.TSG"))
    assert techsas.read(files).sizes["time"] == 0


def test_armstrong_met_all_files_parse():
    files = _require((AR73 / "proc").glob("AR[0-9]*.csv"))
    ds = armstrong.read_met(files)
    assert ds.sizes["time"] > 1000
    assert float(ds.sog.max()) < 10.0


def test_armstrong_gps_first_files_parse():
    files = _require((AR73 / "raw").glob("*.CNAV_3050"))[:5]
    ds = armstrong.read_gps(files)
    assert ds.sizes["time"] > 1000
    assert (np.diff(ds.time.values) > np.timedelta64(0, "ns")).all()


@pytest.mark.parametrize("cruise", ["RR1607", "RR1708"])
def test_revelle_all_files_parse(cruise):
    files = _require((FLEAT / cruise / "data/met/data").glob("*.MET"))
    ds = revelle.read(files)
    assert ds.sizes["time"] > 1000
    assert (np.diff(ds.time.values) > np.timedelta64(0, "ns")).all()
    assert float(ds.lat.min()) >= -90 and float(ds.lat.max()) <= 90

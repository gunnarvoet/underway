import logging

import numpy as np
import pytest

from underway import schema
from underway.parsers import techsas

D = "discovery"


def test_read_gps_first_sample(data):
    first = techsas.read(data / D / "20251024-145312-position-POSMV_GPS.gps").isel(
        time=0
    )
    assert first.time.values.astype("datetime64[s]") == np.datetime64(
        "2025-10-24T14:53:12"
    )
    assert float(first.heading) == pytest.approx(235.6, rel=1e-5)
    assert float(first.cog) == pytest.approx(335.3, rel=1e-5)
    assert float(first.sog) == pytest.approx(0.1 * schema.KNOTS_TO_MS, rel=1e-5)


def test_read_gps_position_variables_renamed(data):
    ds = techsas.read(data / D / "20251024-145312-position-POSMV_GPS.gps")
    assert {"lon", "lat"} <= set(ds.data_vars)
    assert "long" not in ds.variables
    assert "measureTS" not in ds.variables


def test_read_met_first_sample(data):
    first = techsas.read(data / D / "20251024-145312-MET-SURFMET.SURFMETv3").isel(
        time=0
    )
    assert float(first.air_temperature) == pytest.approx(23.22, rel=1e-5)
    assert float(first.relative_humidity) == pytest.approx(74.14, rel=1e-5)
    assert float(first.apparent_wind_speed) == pytest.approx(1.73, rel=1e-5)
    assert float(first.apparent_wind_direction) == pytest.approx(204.0)


def test_read_light_pressure(data):
    first = techsas.read(data / D / "20251024-145312-Light-SURFMET.SURFMETv3").isel(
        time=0
    )
    assert float(first.air_pressure) == pytest.approx(1016.8816, rel=1e-6)


def test_read_file_with_zero_time_steps_returns_empty(data):
    assert techsas.read(data / D / "20251025-000000-SBE45-SBE45.TSG").sizes["time"] == 0


def test_read_zero_byte_file_skipped_with_warning(tmp_path, data, caplog):
    empty = tmp_path / "20251031-000000-position-POSMV_GPS.gps"
    empty.touch()
    with caplog.at_level(logging.WARNING, logger="underway"):
        ds = techsas.read([empty, data / D / "20251024-145312-position-POSMV_GPS.gps"])
    assert ds.sizes["time"] == 20
    assert "zero bytes" in caplog.text


def test_read_every_variable_has_units(data):
    ds = techsas.read(data / D / "20251024-145312-position-POSMV_GPS.gps")
    assert all("units" in ds[v].attrs for v in ds.data_vars)

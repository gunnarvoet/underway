import numpy as np
import pytest

from underway import schema
from underway.parsers import armstrong


def test_read_met_first_row(data):
    first = armstrong.read_met(data / "armstrong/AR230512_0000.csv").isel(time=0)
    assert first.time.values == np.datetime64("2023-05-12T00:00:09.349")
    assert (float(first.lat), float(first.lon)) == pytest.approx((41.524, -70.672))
    assert (float(first.heading), float(first.cog)) == pytest.approx((25.220, 8.100))
    assert float(first.sog) == pytest.approx(0.030 * schema.KNOTS_TO_MS)
    assert (float(first.sss), float(first.sst)) == pytest.approx((0.3367, 13.5880))


def test_read_met_core_uses_starboard_sensor(data):
    first = armstrong.read_met(data / "armstrong/AR230512_0000.csv").isel(time=0)
    assert (float(first.wind_speed), float(first.wind_speed_port)) == pytest.approx(
        (5.1, 4.1)
    )
    assert (
        float(first.wind_direction),
        float(first.wind_direction_port),
    ) == pytest.approx((317.1, 302.0))
    assert (
        float(first.air_temperature),
        float(first.air_temperature_port),
    ) == pytest.approx((15.4, 15.7))
    assert (
        float(first.relative_humidity),
        float(first.relative_humidity_port),
    ) == pytest.approx((58.3, 59.6))
    assert float(first.air_pressure) == pytest.approx(1015.5)


def test_read_met_nan_marker_becomes_nan(data):
    ds = armstrong.read_met(data / "armstrong/AR230512_0000.csv")
    assert np.isnan(ds.water_depth_12khz.values[0])


def test_read_met_row_count(data):
    assert armstrong.read_met(data / "armstrong/AR230512_0000.csv").sizes["time"] == 10


def test_read_met_every_variable_has_units(data):
    ds = armstrong.read_met(data / "armstrong/AR230512_0000.csv")
    assert all("units" in ds[v].attrs for v in ds.data_vars)


def test_read_gps_first_fix(data):
    first = armstrong.read_gps(data / "armstrong/ar20230511_1600.CNAV_3050").isel(
        time=0
    )
    assert first.time.values == np.datetime64("2023-05-11T16:00:00")
    assert float(first.lat) == pytest.approx(41 + 31.441165 / 60)
    assert float(first.lon) == pytest.approx(-(70 + 40.336934 / 60))
    assert float(first.cog) == pytest.approx(331.4)
    assert float(first.sog) == pytest.approx(0.02 * schema.KNOTS_TO_MS)


def test_read_gps_fix_count_and_no_heading(data):
    ds = armstrong.read_gps(data / "armstrong/ar20230511_1600.CNAV_3050")
    assert ds.sizes["time"] == 10
    assert "heading" not in ds

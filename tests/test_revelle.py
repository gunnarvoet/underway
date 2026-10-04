import numpy as np
import pytest

from underway import schema
from underway.parsers import revelle


def test_read_first_row(data):
    first = revelle.read(data / "revelle/170412.MET").isel(time=0)
    assert first.time.values == np.datetime64("2017-04-12T02:12:38")
    assert (float(first.lat), float(first.lon)) == pytest.approx((7.330457, 134.456885))
    assert (float(first.heading), float(first.cog)) == pytest.approx((127.74, 233.4))
    assert float(first.sog) == pytest.approx(0.0 * schema.KNOTS_TO_MS)
    assert (float(first.air_temperature), float(first.air_pressure)) == pytest.approx(
        (29.20, 1007.45)
    )
    assert float(first.relative_humidity) == pytest.approx(73.28)
    assert (float(first.wind_speed), float(first.wind_direction)) == pytest.approx(
        (3.9, 56.5)
    )
    assert (float(first.sst), float(first.sss)) == pytest.approx((31.298, 33.678))


def test_read_missing_value_becomes_nan(data):
    ds = revelle.read(data / "revelle/170412.MET")
    assert np.isnan(ds.water_depth_multibeam.values[0])


def test_read_row_count(data):
    assert revelle.read(data / "revelle/170412.MET").sizes["time"] == 10


def test_read_other_cruise_column_set(data):
    first = revelle.read(data / "revelle/160530.MET").isel(time=0)
    assert first.time.values == np.datetime64("2016-05-30T01:40:58")
    assert float(first.air_temperature) == pytest.approx(29.75)


def test_read_unknown_tags_dropped(data):
    ds = revelle.read(data / "revelle/170412.MET")
    assert "ZO" not in ds and "TT-2" not in ds


def test_read_time_rolls_over_midnight(tmp_path, data):
    lines = (data / "revelle/170412.MET").read_text().splitlines(keepends=True)
    late = lines[4].replace("021238", "235959", 1)
    early = lines[5].replace(lines[5][:6], "000004", 1)
    file = tmp_path / "170412.MET"
    file.write_text("".join(lines[:4] + [late, early]))
    ds = revelle.read(file)
    assert ds.time.values[1] == np.datetime64("2017-04-13T00:00:04")


def test_read_every_variable_has_units(data):
    ds = revelle.read(data / "revelle/170412.MET")
    assert all("units" in ds[v].attrs for v in ds.data_vars)

import logging

import numpy as np
import pytest

from underway import schema
from underway.parsers import seapath

NAME = "ins_seapath_position.20251126T0000Z"


def test_read_first_fix(data):
    first = seapath.read(data / "sikuliaq" / NAME).isel(time=0)
    assert first.time.values == np.datetime64("2025-11-26T00:00:00.300")
    assert float(first.lat) == pytest.approx(21 + 18.947359 / 60)
    assert float(first.lon) == pytest.approx(-(157 + 52.625258 / 60))
    assert float(first.cog) == pytest.approx(4.32)
    assert float(first.sog) == pytest.approx(0.0)
    assert float(first.heading) == pytest.approx(26.63)
    assert (
        float(first["roll"]),
        float(first["pitch"]),
        float(first["heave"]),
    ) == pytest.approx((0.08, 0.24, 0.0))


def test_read_fix_count(data):
    assert seapath.read(data / "sikuliaq" / NAME).sizes["time"] == 10


def test_read_sog_in_meters_per_second(tmp_path, data):
    text = (
        (data / "sikuliaq" / NAME)
        .read_text()
        .replace(
            "$GPVTG,4.32,T,354.83,M,0.0,N,0.0,K,D*2A",
            "$GPVTG,4.32,T,354.83,M,10.0,N,18.5,K,D*2A",
        )
    )
    file = tmp_path / NAME
    file.write_text(text)
    assert float(seapath.read(file).sog[0]) == pytest.approx(10 * schema.KNOTS_TO_MS)


def test_read_missing_sentence_affects_only_its_fix(tmp_path, data):
    reference = seapath.read(data / "sikuliaq" / NAME)
    lines = (data / "sikuliaq" / NAME).read_text().splitlines(keepends=True)
    gga = [i for i, line in enumerate(lines) if "$GPGGA" in line]
    del lines[gga[2]]
    file = tmp_path / NAME
    file.write_text("".join(lines))
    ds = seapath.read(file)
    assert ds.sizes["time"] == 10
    assert np.isnan(ds.lat.values[2])
    np.testing.assert_array_equal(ds.lat.values[3:], reference.lat.values[3:])
    np.testing.assert_array_equal(ds.heading.values, reference.heading.values)


def test_read_lines_before_first_zda_dropped_and_counted(tmp_path, data, caplog):
    lines = (data / "sikuliaq" / NAME).read_text().splitlines(keepends=True)
    first = next(i for i, line in enumerate(lines) if "$GPZDA" in line)
    del lines[first]
    file = tmp_path / NAME
    file.write_text("".join(lines))
    with caplog.at_level(logging.WARNING, logger="underway"):
        ds = seapath.read(file)
    assert ds.sizes["time"] == 9
    assert "dropped 7" in caplog.text


def test_read_truncated_last_line(tmp_path, data):
    text = (data / "sikuliaq" / NAME).read_text()
    file = tmp_path / NAME
    file.write_text(text + "ins_seapath_position\t2025-11-26T00:00:10.33Z\t$GPZDA,0000")
    ds = seapath.read(file)
    assert ds.sizes["time"] == 10


def test_read_zero_byte_file_returns_empty(tmp_path):
    file = tmp_path / NAME
    file.touch()
    assert seapath.read(file).sizes["time"] == 0


def test_read_every_variable_has_units(data):
    ds = seapath.read(data / "sikuliaq" / NAME)
    assert all("units" in ds[v].attrs for v in ds.data_vars)

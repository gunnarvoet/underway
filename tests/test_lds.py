import logging

import numpy as np
import pytest

from underway.parsers import lds


def test_read_tsg_first_sample(data):
    ds = lds.read(data / "sikuliaq/tsg_sbe45_fwd.20251127T0000Z", stream="tsg")
    first = ds.isel(time=0)
    assert first.time.values == np.datetime64("2025-11-27T00:00:03.626800")
    values = [float(first[v]) for v in ["sst", "conductivity", "sss", "sound_speed"]]
    assert values == pytest.approx([26.9030, 5.51282, 35.0429, 1538.914])


def test_read_tsg_row_count(data):
    ds = lds.read(data / "sikuliaq/tsg_sbe45_fwd.20251127T0000Z", stream="tsg")
    assert ds.sizes["time"] == 20


def test_read_wind_first_sample_uses_meters_per_second(data):
    ds = lds.read(
        data / "sikuliaq/wind_gill_fwdmast_true.20251125T0000Z", stream="wind"
    )
    first = ds.isel(time=0)
    assert (float(first.wind_direction), float(first.wind_speed)) == pytest.approx(
        (72.6, 4.4)
    )


def test_read_air_pressure_converted_to_hpa(data):
    ds = lds.read(data / "sikuliaq/met_met4a_fwdmast.20251125T0000Z", stream="air")
    first = ds.isel(time=0)
    assert float(first.air_pressure) == pytest.approx(1012.158)
    assert (
        float(first.air_temperature),
        float(first.relative_humidity),
    ) == pytest.approx((27.41, 55.36))


def test_read_every_variable_has_units(data):
    ds = lds.read(data / "sikuliaq/tsg_sbe45_fwd.20251127T0000Z", stream="tsg")
    assert all("units" in ds[v].attrs for v in ds.data_vars)


def test_read_malformed_and_truncated_lines_dropped_and_counted(tmp_path, caplog):
    file = tmp_path / "tsg_sbe45_fwd.20251127T0000Z"
    file.write_text(
        "# header\n"
        "tsg_sbe45_fwd\t2025-11-27T00:00:03.6268Z\t 26.9030,  5.51282,  35.0429, 1538.914\n"
        "tsg_sbe45_fwd\t2025-11-27T00:00:08.6266Z\t 26.9018,  abc,  35.0435, 1538.911\n"
        "tsg_sbe45_fwd\t2025-11-27T00:00:13.6"
    )
    with caplog.at_level(logging.WARNING, logger="underway"):
        ds = lds.read(file, stream="tsg")
    assert ds.sizes["time"] == 1
    assert "dropped 2" in caplog.text


def test_read_header_only_file_returns_empty(tmp_path):
    file = tmp_path / "tsg_sbe45_fwd.20251127T0000Z"
    file.write_text("# header\n# more header\n")
    assert lds.read(file, stream="tsg").sizes["time"] == 0


def test_read_zero_byte_file_returns_empty(tmp_path):
    file = tmp_path / "tsg_sbe45_fwd.20251127T0000Z"
    file.touch()
    assert lds.read(file, stream="tsg").sizes["time"] == 0


def test_read_unknown_stream_raises(data):
    with pytest.raises(ValueError, match="unknown stream"):
        lds.read(data / "sikuliaq/tsg_sbe45_fwd.20251127T0000Z", stream="nope")

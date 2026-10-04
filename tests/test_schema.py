import numpy as np
import pytest
import xarray as xr

from underway import schema


def _ds(times, **variables):
    time = np.array(times, dtype="datetime64[ns]")
    return xr.Dataset(
        {
            name: ("time", np.asarray(values, dtype=float))
            for name, values in variables.items()
        },
        coords={"time": time},
    )


def test_conform_core_variable_gets_schema_attrs():
    ds = schema.conform(_ds(["2025-01-01"], sog=[1.0]))
    assert ds.sog.attrs == {"long_name": "speed over ground", "units": "m/s"}


def test_conform_extra_variable_gets_given_attrs():
    ds = schema.conform(
        _ds(["2025-01-01"], roll=[1.0]), extra={"roll": ("roll", "deg")}
    )
    assert ds["roll"].attrs == {"long_name": "roll", "units": "deg"}


def test_conform_extra_for_absent_variable_is_ignored():
    ds = schema.conform(_ds(["2025-01-01"], sog=[1.0]), extra={"roll": ("roll", "deg")})
    assert list(ds.data_vars) == ["sog"]


def test_empty_has_zero_length_time():
    assert schema.empty().sizes["time"] == 0


def test_combine_overlapping_files_sorted_and_unique():
    a = _ds(["2025-01-01T00:00:02", "2025-01-01T00:00:00"], sst=[2.0, 0.0])
    b = _ds(["2025-01-01T00:00:02", "2025-01-01T00:00:03"], sst=[2.0, 3.0])
    out = schema.combine([a, b])
    assert out.sst.values.tolist() == [0.0, 2.0, 3.0]


def test_combine_drops_nat():
    a = _ds(["2025-01-01T00:00:00", "NaT"], sst=[0.0, 9.0])
    assert schema.combine([a]).sizes["time"] == 1


def test_combine_skips_empty_datasets():
    a = _ds(["2025-01-01"], sst=[1.0])
    assert schema.combine([schema.empty(), a]).sizes["time"] == 1


def test_combine_no_datasets_returns_empty():
    assert schema.combine([]).sizes["time"] == 0


def _angle_close(value, target):
    difference = abs(value - target) % 360
    return min(difference, 360 - difference) < 1e-6


def test_bin_average_scalar_mean():
    ds = _ds(["2025-01-01T00:00:10", "2025-01-01T00:00:50"], sst=[10.0, 12.0])
    assert float(schema.bin_average(ds, "1min").sst[0]) == pytest.approx(11.0)


def test_bin_average_heading_across_north():
    ds = _ds(["2025-01-01T00:00:10", "2025-01-01T00:00:50"], heading=[350.0, 10.0])
    assert _angle_close(float(schema.bin_average(ds, "1min").heading[0]), 0.0)


def test_bin_average_wind_as_vector_pair():
    ds = _ds(
        ["2025-01-01T00:00:10", "2025-01-01T00:00:50"],
        wind_speed=[5.0, 5.0],
        wind_direction=[350.0, 10.0],
    )
    out = schema.bin_average(ds, "1min")
    assert _angle_close(float(out.wind_direction[0]), 0.0)
    assert float(out.wind_speed[0]) == pytest.approx(5.0 * np.cos(np.deg2rad(10.0)))


def test_bin_average_drops_non_core_variables():
    ds = _ds(["2025-01-01T00:00:10"], sst=[1.0], roll=[2.0])
    assert list(schema.bin_average(ds, "1min").data_vars) == ["sst"]


def test_bin_average_sets_units():
    ds = _ds(["2025-01-01T00:00:10"], sst=[1.0])
    assert schema.bin_average(ds, "1min").sst.attrs["units"] == "°C"


def test_combine_duplicate_time_prefers_row_with_data():
    a = _ds(["2025-01-01T00:00:00"], sst=[np.nan])
    b = _ds(["2025-01-01T00:00:00"], sst=[7.0])
    assert schema.combine([a, b]).sst.values.tolist() == [7.0]


def test_bin_average_direction_below_360():
    ds = _ds(
        ["2025-01-01T00:00:10", "2025-01-01T00:00:50"],
        heading=[350.0, 10.0],
        wind_speed=[5.0, 5.0],
        wind_direction=[350.0, 10.0],
    )
    out = schema.bin_average(ds, "1min")
    assert float(out.heading[0]) < 360.0
    assert float(out.wind_direction[0]) < 360.0

import numpy as np
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

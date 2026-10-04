import numpy as np
import pytest
import xarray as xr

from underway import cache, schema


class CountingParser:
    """Parser stand-in: one time step per line of the raw file."""

    def __init__(self):
        self.calls = 0

    def __call__(self, files):
        self.calls += 1
        n = len(files[0].read_text().splitlines())
        if n == 0:
            return schema.empty()
        time = np.datetime64("2025-01-01", "ns") + np.arange(n) * np.timedelta64(1, "s")
        return xr.Dataset(
            {"sst": ("time", np.arange(n, dtype=float))}, coords={"time": time}
        )


@pytest.fixture
def raw(tmp_path):
    file = tmp_path / "raw" / "tsg.20250101T0000Z"
    file.parent.mkdir()
    file.write_text("a\nb\n")
    return file


@pytest.fixture
def product(tmp_path, raw):
    proc = tmp_path / "proc"
    proc.mkdir()
    return cache.product_path(proc, "SKQ202521S", "tsg", raw)


def test_product_path_keeps_full_raw_name(product):
    assert product.name == "skq202521s_tsg_tsg.20250101T0000Z.nc"


def test_update_new_raw_file_writes_product_with_size(raw, product):
    assert cache.update(raw, product, CountingParser()) is True
    with xr.open_dataset(product) as ds:
        assert (ds.sizes["time"], ds.attrs["raw_size"], ds.attrs["raw_name"]) == (
            2,
            raw.stat().st_size,
            raw.name,
        )


def test_update_unchanged_raw_file_not_parsed_again(raw, product):
    parser = CountingParser()
    cache.update(raw, product, parser)
    assert cache.update(raw, product, parser) is False
    assert parser.calls == 1


def test_update_grown_raw_file_parsed_again(raw, product):
    parser = CountingParser()
    cache.update(raw, product, parser)
    raw.write_text("a\nb\nc\n")
    cache.update(raw, product, parser)
    with xr.open_dataset(product) as ds:
        assert ds.sizes["time"] == 3


def test_update_force_parses_unchanged_file(raw, product):
    parser = CountingParser()
    cache.update(raw, product, parser)
    cache.update(raw, product, parser, force=True)
    assert parser.calls == 2


def test_update_empty_raw_file_cached_as_empty_product(raw, product):
    raw.write_text("")
    parser = CountingParser()
    cache.update(raw, product, parser)
    assert cache.update(raw, product, parser) is False


def test_update_failed_write_keeps_old_product(raw, product, monkeypatch):
    parser = CountingParser()
    cache.update(raw, product, parser)
    raw.write_text("a\nb\nc\n")

    def fail(self, *args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(xr.Dataset, "to_netcdf", fail)
    with pytest.raises(OSError, match="disk full"):
        cache.update(raw, product, parser)
    monkeypatch.undo()
    with xr.open_dataset(product) as ds:
        assert ds.sizes["time"] == 2
    assert list(product.parent.glob("*.tmp")) == []


def test_product_path_differs_for_equal_names_in_two_directories(tmp_path):
    raw_dir = tmp_path / "raw"
    a = cache.product_path(
        tmp_path, "X", "s", raw_dir / "a" / "data.txt", raw_dir=raw_dir
    )
    b = cache.product_path(
        tmp_path, "X", "s", raw_dir / "b" / "data.txt", raw_dir=raw_dir
    )
    assert a != b

from underway import transfer
import logging
import shutil

import pytest

import underway as uw
from underway.cruise import SyncError

needs_rsync = pytest.mark.skipif(
    shutil.which("rsync") is None, reason="rsync not installed"
)

CRUISE_ID = "SKQ202521S"
STREAMS = {
    "ins_seapath_position": "ins_seapath_position.20251126T0000Z",
    "tsg_sbe45_fwd": "tsg_sbe45_fwd.20251127T0000Z",
    "wind_gill_fwdmast_true": "wind_gill_fwdmast_true.20251125T0000Z",
    "met_met4a_fwdmast": "met_met4a_fwdmast.20251125T0000Z",
}


@pytest.fixture
def server(tmp_path, data):
    """Mount root with the Sikuliaq LDS layout, filled from the test data."""
    root = tmp_path / "Volumes"
    lds = root / "CruiseData" / CRUISE_ID / "lds" / "raw"
    for stream, name in STREAMS.items():
        (lds / stream).mkdir(parents=True)
        shutil.copy(data / "sikuliaq" / name, lds / stream / name)
    return root


@pytest.fixture
def cruise(tmp_path, server, monkeypatch):
    # keep sync from calling osascript when the tests run on macOS
    monkeypatch.setattr("underway.mount.sys.platform", "linux")
    return uw.Cruise("sikuliaq", CRUISE_ID, tmp_path / "local", mount_root=server)


def test_constructor_creates_no_directories(tmp_path, cruise):
    assert not (tmp_path / "local").exists()


def test_constructor_unknown_ship_raises(tmp_path):
    with pytest.raises(ValueError, match="unknown ship"):
        uw.Cruise("titanic", "X", tmp_path)


def test_path_default_layout(tmp_path, cruise):
    assert cruise.path("gps", "raw") == tmp_path / "local" / "gps" / "raw"
    assert cruise.path("gps") == tmp_path / "local" / "gps"


def test_path_unknown_source_raises(cruise):
    with pytest.raises(KeyError, match="unknown source"):
        cruise.path("nope")


def test_path_local_is_exact_raw_directory(tmp_path, cruise):
    cruise.add_source("ww", drive="sci", remote="ww", local=tmp_path / "elsewhere")
    assert cruise.path("ww", "raw") == tmp_path / "elsewhere"
    assert cruise.path("ww", "proc") == tmp_path / "local" / "ww" / "proc"


def test_path_create_makes_directory(cruise):
    assert cruise.path("ctd", "proc", create=True).is_dir()


def test_remote_path_fills_cruise_id(server, cruise):
    assert (
        cruise.remote_path("ctd") == server / "CruiseData" / CRUISE_ID / "ctd" / "raw"
    )


def test_add_source_with_server_registers_drive(cruise):
    cruise.add_source("ww", drive="sci", remote="ww", server="sci.example.edu")
    assert cruise.servers["sci"] == "sci.example.edu"


@needs_rsync
def test_sync_selected_source_copies_files(cruise):
    cruise.sync("gps")
    assert (cruise.path("gps", "raw") / STREAMS["ins_seapath_position"]).exists()
    assert not cruise.path("tsg", "raw").exists()


@needs_rsync
def test_sync_missing_remote_reported_after_other_sources_synced(cruise):
    with pytest.raises(SyncError) as err:
        cruise.sync()
    assert "sadcp" in str(err.value) and "ctd" in str(err.value)
    assert "gps" not in str(err.value)
    assert (cruise.path("tsg", "raw") / STREAMS["tsg_sbe45_fwd"]).exists()


@needs_rsync
def test_sync_skip_leaves_source_out(cruise):
    cruise.sync(skip=("sadcp", "ctd"))
    assert (cruise.path("gps", "raw") / STREAMS["ins_seapath_position"]).exists()


def test_sync_unknown_source_raises_before_any_transfer(cruise):
    with pytest.raises(KeyError, match="unknown source"):
        cruise.sync("gps", "nope")
    assert not cruise.path("gps", "raw").exists()


@needs_rsync
def test_read_after_sync_returns_parsed_record(cruise):
    cruise.sync("gps")
    gps = cruise.read("gps")
    assert gps.sizes["time"] == 10
    assert list(cruise.path("gps", "proc").glob("skq202521s_gps_*.nc"))


@needs_rsync
def test_read_result_has_no_cache_attributes(cruise):
    cruise.sync("gps")
    assert "raw_size" not in cruise.read("gps").attrs


def test_read_before_sync_returns_empty(cruise):
    assert cruise.read("gps").sizes["time"] == 0


@needs_rsync
def test_read_grown_raw_file_returns_new_fix(cruise):
    cruise.sync("gps")
    cruise.read("gps")
    raw = cruise.path("gps", "raw") / STREAMS["ins_seapath_position"]
    with open(raw, "a") as f:
        f.write(
            "ins_seapath_position\t2025-11-26T00:00:10.33Z\t$GPZDA,000010.30,26,11,2025,,*64\n"
        )
        f.write("ins_seapath_position\t2025-11-26T00:00:10.62Z\t$GPHDT,27.00,T*34\n")
    assert cruise.read("gps").sizes["time"] == 11


@needs_rsync
def test_read_uses_products_when_raw_files_are_dangling_links(cruise):
    cruise.sync("gps")
    cruise.read("gps")
    raw = cruise.path("gps", "raw") / STREAMS["ins_seapath_position"]
    raw.unlink()
    raw.symlink_to(raw.parent / "annex-object-not-present")
    assert cruise.read("gps").sizes["time"] == 10


@needs_rsync
def test_read_reparse_rewrites_product(cruise):
    cruise.sync("gps")
    cruise.read("gps")
    product = next(cruise.path("gps", "proc").glob("*.nc"))
    before = product.stat().st_mtime_ns
    cruise.read("gps", reparse=True)
    assert product.stat().st_mtime_ns > before


def test_read_source_without_parser_raises(cruise):
    with pytest.raises(ValueError, match="no parser"):
        cruise.read("ctd")


def test_local_dir_tilde_expanded(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    c = uw.Cruise("sikuliaq", CRUISE_ID, "~/cruise")
    assert c.path("gps") == tmp_path / "cruise" / "gps"


@needs_rsync
def test_read_met_merges_core_variables_on_one_minute_grid(cruise):
    cruise.sync("gps", "tsg", "wind", "air")
    met = cruise.read_met()
    assert {"lon", "lat", "heading", "sst", "sss", "wind_speed", "air_pressure"} <= set(
        met.data_vars
    )
    assert "roll" not in met
    assert (met.time.dt.second == 0).all()


def test_read_met_before_sync_returns_empty(cruise):
    assert cruise.read_met().sizes["time"] == 0


def _flaky_parser(files):
    if any(f.name.endswith("bad") for f in files):
        raise ValueError("boom")
    return uw.parsers.seapath.read(files)


@needs_rsync
@pytest.mark.parametrize("cached", [True, False])
def test_read_unparseable_file_logged_and_other_files_returned(cruise, caplog, cached):
    cruise.add_source(
        "gps",
        drive="CruiseData",
        remote=CRUISE_ID + "/lds/raw/ins_seapath_position",
        pattern="ins_seapath_position.*",
        parser=_flaky_parser,
        cache=cached,
    )
    cruise.sync("gps")
    (cruise.path("gps", "raw") / "ins_seapath_position.bad").write_text("x")
    with caplog.at_level(logging.WARNING, logger="underway"):
        gps = cruise.read("gps")
    assert gps.sizes["time"] == 10
    assert "ins_seapath_position.bad" in caplog.text


@needs_rsync
def test_sync_os_error_in_one_source_does_not_stop_the_others(cruise, monkeypatch):
    real = transfer.TRANSFER["rsync"]

    def flaky(src_dir, dst_dir, **kwargs):
        if src_dir.name == "ins_seapath_position":
            raise FileNotFoundError("share dropped")
        return real(src_dir, dst_dir, **kwargs)

    monkeypatch.setitem(transfer.TRANSFER, "rsync", flaky)
    with pytest.raises(SyncError, match="share dropped"):
        cruise.sync("gps", "tsg")
    assert (cruise.path("tsg", "raw") / STREAMS["tsg_sbe45_fwd"]).exists()

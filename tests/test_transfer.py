import shutil

import pytest

from underway import transfer

needs_rsync = pytest.mark.skipif(
    shutil.which("rsync") is None, reason="rsync not installed"
)


@pytest.fixture
def remote(tmp_path):
    src = tmp_path / "remote"
    (src / "wh300" / "contour").mkdir(parents=True)
    (src / "a.20250101").write_text("a")
    (src / "a.20250102").write_text("aa")
    (src / "b.20250101").write_text("b")
    (src / "wh300" / "contour" / "wh300.nc").write_text("nc")
    (src / "wh300" / "contour" / "notes.txt").write_text("txt")
    return src


@pytest.fixture
def local(tmp_path):
    dst = tmp_path / "local"
    dst.mkdir()
    return dst


def _names(directory):
    return sorted(
        p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()
    )


@needs_rsync
def test_rsync_pattern_selects_files(remote, local):
    transfer.rsync(remote, local, pattern="a.*")
    assert _names(local) == ["a.20250101", "a.20250102"]


@needs_rsync
def test_rsync_exclude_skips_file(remote, local):
    transfer.rsync(remote, local, pattern="a.*", exclude=("a.20250102",))
    assert _names(local) == ["a.20250101"]


@needs_rsync
def test_rsync_nested_pattern_keeps_structure(remote, local):
    transfer.rsync(remote, local, pattern="*/contour/*.nc")
    assert _names(local) == ["wh300/contour/wh300.nc"]


@needs_rsync
def test_rsync_flat_pattern_does_not_descend(remote, local):
    (remote / "wh300" / "a.20250103").write_text("nested")
    transfer.rsync(remote, local, pattern="a.*")
    assert _names(local) == ["a.20250101", "a.20250102"]


@needs_rsync
def test_rsync_default_pattern_copies_everything(remote, local):
    transfer.rsync(remote, local)
    assert len(_names(local)) == 5


@needs_rsync
def test_rsync_missing_source_raises(tmp_path, local):
    with pytest.raises(transfer.TransferError, match="rsync"):
        transfer.rsync(tmp_path / "missing", local)


def test_copy_pattern_selects_files(remote, local):
    transfer.copy(remote, local, pattern="a.*")
    assert _names(local) == ["a.20250101", "a.20250102"]


def test_copy_nested_pattern_keeps_structure(remote, local):
    transfer.copy(remote, local, pattern="*/contour/*.nc")
    assert _names(local) == ["wh300/contour/wh300.nc"]


def test_copy_same_size_file_not_copied_again(remote, local):
    transfer.copy(remote, local, pattern="a.*")
    (local / "a.20250101").write_text("x")
    transfer.copy(remote, local, pattern="a.*")
    assert (local / "a.20250101").read_text() == "x"


def test_copy_changed_size_file_copied_again(remote, local):
    transfer.copy(remote, local, pattern="a.*")
    (remote / "a.20250101").write_text("grown")
    transfer.copy(remote, local, pattern="a.*")
    assert (local / "a.20250101").read_text() == "grown"


def test_copy_exclude_skips_file(remote, local):
    transfer.copy(remote, local, pattern="a.*", exclude=("a.20250102",))
    assert _names(local) == ["a.20250101"]


def test_copy_unreadable_file_reported_after_others_copied(remote, local, monkeypatch):
    real = shutil.copyfile

    def flaky(src, dst):
        if src.name == "a.20250101":
            raise PermissionError("denied")
        return real(src, dst)

    monkeypatch.setattr(transfer.shutil, "copyfile", flaky)
    with pytest.raises(transfer.TransferError, match="a.20250101"):
        transfer.copy(remote, local, pattern="a.*")
    assert _names(local) == ["a.20250102"]


def test_make_readonly_sets_mode_0400(remote):
    transfer.make_readonly(remote)
    assert (remote / "a.20250101").stat().st_mode & 0o777 == 0o400


def test_copy_exclude_matches_directory_name(remote, local):
    transfer.copy(remote, local, pattern="*/contour/*", exclude=("wh300",))
    assert _names(local) == []

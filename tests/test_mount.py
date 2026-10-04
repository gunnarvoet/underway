import subprocess

import pytest

from underway import mount


def test_connect_mounted_drive_does_nothing(tmp_path, monkeypatch):
    (tmp_path / "CruiseData").mkdir()
    monkeypatch.setattr(mount.sys, "platform", "linux")
    mount.connect({"CruiseData": "data.example.edu"}, tmp_path)


def test_connect_missing_drive_off_macos_names_drive_and_server(tmp_path, monkeypatch):
    monkeypatch.setattr(mount.sys, "platform", "linux")
    with pytest.raises(mount.MountError, match="//data.example.edu/CruiseData"):
        mount.connect({"CruiseData": "data.example.edu"}, tmp_path)


def test_connect_lists_every_missing_drive(tmp_path, monkeypatch):
    monkeypatch.setattr(mount.sys, "platform", "linux")
    with pytest.raises(mount.MountError) as err:
        mount.connect({"a": "s1", "b": "s2"}, tmp_path)
    assert "//s1/a" in str(err.value) and "//s2/b" in str(err.value)


def test_connect_on_macos_runs_osascript(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(mount.sys, "platform", "darwin")
    monkeypatch.setattr(mount.subprocess, "run", fake_run)
    mount.connect({"CruiseData": "data.example.edu"}, tmp_path)
    assert calls == [
        ["osascript", "-e", 'mount volume "smb://data.example.edu/CruiseData"']
    ]


def test_connect_on_macos_failed_mount_raises(tmp_path, monkeypatch):
    def fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="no route")

    monkeypatch.setattr(mount.sys, "platform", "darwin")
    monkeypatch.setattr(mount.subprocess, "run", fake_run)
    with pytest.raises(mount.MountError, match="no route"):
        mount.connect({"CruiseData": "data.example.edu"}, tmp_path)

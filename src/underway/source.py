"""Description of one stream of files on the ship server."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Source:
    """One stream of files to sync and, optionally, to parse.

    Parameters
    ----------
    name : str
        Key used in `Cruise.sync`, `Cruise.read`, and `Cruise.path`.
    drive : str
        Share name on the ship server, mounted under the cruise's mount root.
    remote : str
        Directory below the drive. May contain ``{cruise_id}``.
    pattern : str, optional
        Glob for files below `remote`. May contain directories.
    exclude : tuple of str, optional
        File name patterns to skip during sync.
    transfer : {"rsync", "copy"}, optional
        Transfer method.
    parser : callable, optional
        Function ``files -> xr.Dataset``. None for sync-only sources.
    cache : bool, optional
        Cache one netCDF product per raw file. False for raw formats that
        are netCDF already.
    met : bool, optional
        Include this source in `Cruise.read_met`.
    readonly : bool, optional
        Set synced files to mode 0o400.
    local : Path, optional
        Exact local directory for the raw files. Defaults to
        ``<local_dir>/<name>/raw``. Several sources may share one directory
        if their patterns differ.
    """

    name: str
    drive: str
    remote: str
    pattern: str = "*"
    exclude: tuple[str, ...] = ()
    transfer: str = "rsync"
    parser: Callable | None = None
    cache: bool = True
    met: bool = False
    readonly: bool = False
    local: Path | None = None

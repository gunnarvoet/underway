"""Cruise: local paths, sync from the ship server, cached reads."""

import logging
from pathlib import Path

import xarray as xr

from . import cache, mount, schema, transfer
from .ships import SHIPS
from .source import Source

log = logging.getLogger("underway")


class SyncError(Exception):
    """One or more sources failed to sync."""


class Cruise:
    """One cruise on one ship.

    Parameters
    ----------
    ship : {"sikuliaq", "discovery", "armstrong", "revelle"}
        Selects the default sources and servers from `underway.ships.SHIPS`.
    cruise_id : str
        Cruise identifier as used in paths on the ship server.
    local_dir : path
        Local cruise data directory. Created when first needed.
    mount_root : path, optional
        Directory under which the ship's drives are mounted.

    Examples
    --------
    >>> c = Cruise("sikuliaq", "SKQ202521S", "~/data/cruise2")
    >>> c.sync("gps", "tsg")
    >>> gps = c.read("gps")
    """

    def __init__(self, ship, cruise_id, local_dir, mount_root="/Volumes"):
        if ship not in SHIPS:
            raise ValueError(f"unknown ship {ship!r}, choose from {sorted(SHIPS)}")
        self.ship = SHIPS[ship]
        self.cruise_id = cruise_id
        self.local_dir = Path(local_dir).expanduser()
        self.mount_root = Path(mount_root)
        self.sources = {source.name: source for source in self.ship.sources}
        self.servers = dict(self.ship.servers)

    def add_source(self, name, drive, remote, server=None, **kwargs):
        """Add a source, or replace a default source of the same name.

        Parameters
        ----------
        name, drive, remote
            See `Source`.
        server : str, optional
            Server that holds `drive`, if the drive is not one of the ship's
            default drives.
        **kwargs
            Remaining `Source` fields.

        Returns
        -------
        Source
        """
        source = Source(name=name, drive=drive, remote=remote, **kwargs)
        self.sources[name] = source
        if server is not None:
            self.servers[drive] = server
        return source

    def _source(self, name):
        if name not in self.sources:
            raise KeyError(
                f"unknown source {name!r}, choose from {sorted(self.sources)}"
            )
        return self.sources[name]

    def path(self, name, kind=None, create=False):
        """Return the local directory of a source.

        Parameters
        ----------
        name : str
            Source name.
        kind : {None, "raw", "proc"}, optional
            ``"raw"`` is `Source.local` if set, else ``<local_dir>/<name>/raw``.
            ``"proc"`` is ``<local_dir>/<name>/proc``. None returns
            ``<local_dir>/<name>``.
        create : bool, optional
            Create the directory if it does not exist.
        """
        source = self._source(name)
        base = self.local_dir / name
        if kind is None:
            directory = base
        elif kind == "raw":
            directory = (
                Path(source.local).expanduser() if source.local else base / "raw"
            )
        elif kind == "proc":
            directory = base / "proc"
        else:
            raise ValueError(f"kind must be 'raw' or 'proc', got {kind!r}")
        if create:
            directory.mkdir(parents=True, exist_ok=True)
        return directory

    def remote_path(self, name):
        """Return the directory of a source on the mounted ship share."""
        source = self._source(name)
        remote = source.remote.format(cruise_id=self.cruise_id)
        return self.mount_root / source.drive / remote

    def connect(self):
        """Mount the ship's drives. See `underway.mount.connect`."""
        mount.connect(self.servers, self.mount_root)

    def sync(self, *names, skip=(), verbose=False):
        """Copy raw files from the ship server.

        Parameters
        ----------
        *names : str
            Sources to sync. All sources if none are given.
        skip : tuple of str, optional
            Sources to leave out.
        verbose : bool, optional
            Print transferred files.

        Raises
        ------
        SyncError
            After every source was attempted, naming each one that failed.
        """
        names = names or tuple(self.sources)
        sources = [self._source(name) for name in names if name not in skip]
        try:
            self.connect()
        except mount.MountError as err:
            log.warning("%s", err)
        failures = {}
        for source in sources:
            remote = self.remote_path(source.name)
            local = self.path(source.name, "raw")
            try:
                if not remote.is_dir():
                    raise transfer.TransferError(f"remote directory {remote} not found")
                local.mkdir(parents=True, exist_ok=True)
                transfer.TRANSFER[source.transfer](
                    remote,
                    local,
                    pattern=source.pattern,
                    exclude=source.exclude,
                    verbose=verbose,
                )
                if source.readonly:
                    transfer.make_readonly(local)
            except transfer.TransferError as err:
                failures[source.name] = str(err)
        if failures:
            lines = [f"{name}: {message}" for name, message in failures.items()]
            raise SyncError("sync failed for\n" + "\n".join(lines))

    def read(self, name, reparse=False):
        """Parse a source and return its full record.

        Raw files whose size changed since the last call are parsed again.
        All other days come from the netCDF products in ``<source>/proc``.

        Parameters
        ----------
        name : str
            Source name. The source needs a parser.
        reparse : bool, optional
            Parse every raw file again, for use after a parser change.

        Returns
        -------
        xr.Dataset
            Sorted by time, each time stamp once. Empty if there is no data.
        """
        source = self._source(name)
        if source.parser is None:
            raise ValueError(f"source {name!r} has no parser")
        raw_dir = self.path(name, "raw")
        raw_files = sorted(f for f in raw_dir.glob(source.pattern) if f.is_file())
        if not source.cache:
            return source.parser(raw_files)
        proc_dir = self.path(name, "proc")
        proc_dir.mkdir(parents=True, exist_ok=True)
        for raw in raw_files:
            product = cache.product_path(proc_dir, self.cruise_id, name, raw)
            cache.update(raw, product, source.parser, force=reparse)
        datasets = []
        for product in sorted(proc_dir.glob(f"{self.cruise_id.lower()}_{name}_*.nc")):
            with xr.open_dataset(product) as ds:
                datasets.append(ds.load())
        out = schema.combine(datasets)
        out.attrs.pop("raw_name", None)
        out.attrs.pop("raw_size", None)
        return out

    def read_met(self, freq="1min"):
        """Return the core variables of all met sources on one time grid.

        Every source with ``met=True`` is read and bin-averaged with
        `underway.schema.bin_average`. Where two sources provide the same
        variable, the source listed first is used.

        Parameters
        ----------
        freq : str, optional
            Bin width as a pandas frequency string.

        Returns
        -------
        xr.Dataset
            Core variables on a regular time grid. Empty if there is no data.
        """
        parts = []
        for name, source in self.sources.items():
            if not source.met or source.parser is None:
                continue
            ds = self.read(name)
            if ds.sizes["time"] > 0:
                parts.append(schema.bin_average(ds, freq))
        if not parts:
            return schema.empty()
        return xr.merge(parts, compat="override", join="outer")

"""Sikuliaq Seapath position files (NMEA sentences in an LDS wrapper)."""

from .. import schema
from . import nmea
from ._common import as_paths, read_lds_lines


def read(files):
    """Read Sikuliaq ``ins_seapath_position.*`` files.

    Parameters
    ----------
    files : path or iterable of paths

    Returns
    -------
    xr.Dataset
        `lat`, `lon`, `cog`, `sog`, `heading`, `roll`, `pitch`, `heave` on a
        `time` dimension, one step per position fix.

    Examples
    --------
    >>> gps = read(sorted(raw_dir.glob("ins_seapath_position.*")))
    """
    datasets = []
    for file in as_paths(files):
        lines, partial = read_lds_lines(file)
        datasets.append(nmea.parse(lines["payload"], name=file.name, dropped=partial))
    return schema.combine(datasets)

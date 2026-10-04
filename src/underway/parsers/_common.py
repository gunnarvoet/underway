"""Helpers shared by the parsers."""

from pathlib import Path

import pandas as pd


def as_paths(files):
    """Return a sorted list of paths from one path or an iterable of paths."""
    if isinstance(files, (str, Path)):
        files = [files]
    return sorted(Path(f) for f in files)


def read_complete(file):
    """Read a text file up to its last newline.

    A file that is still being written can end in a line that is cut inside
    a field. Such a line would parse into wrong values, so it is cut off.

    Parameters
    ----------
    file : path

    Returns
    -------
    data : bytes
        File content up to and including the last newline.
    dropped : int
        1 if a partial last line was cut off, else 0.
    """
    data = Path(file).read_bytes()
    if not data or data.endswith(b"\n"):
        return data, 0
    return data[: data.rfind(b"\n") + 1], 1


def read_lds_lines(file):
    """Read a Sikuliaq LDS text file into string columns `time` and `payload`.

    Lines are ``<log name>\\t<ISO time>\\t<payload>``. Header lines start with
    ``#``. A partial last line is cut off. Lines with fewer than three
    fields give missing values.

    Returns
    -------
    lines : pd.DataFrame
        String columns `time` and `payload`, one row per data line.
    dropped : int
        Number of partial lines cut off (0 or 1).
    """
    data, dropped = read_complete(file)
    text = pd.Series(data.decode(errors="replace").splitlines(), dtype=str)
    text = text[~text.str.startswith("#") & (text.str.len() > 0)]
    if text.empty:
        return pd.DataFrame({"time": [], "payload": []}, dtype=str), dropped
    fields = text.str.split("\t", n=2, expand=True).reindex(columns=range(3))
    # object dtype keeps the .str accessor usable when a column is all missing
    fields = fields.astype(object)
    lines = pd.DataFrame({"time": fields[1], "payload": fields[2]})
    return lines.reset_index(drop=True), dropped

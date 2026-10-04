"""Helpers shared by the parsers."""

from pathlib import Path

import pandas as pd


def as_paths(files):
    """Return a sorted list of paths from one path or an iterable of paths."""
    if isinstance(files, (str, Path)):
        files = [files]
    return sorted(Path(f) for f in files)


def read_lds_lines(file):
    """Read a Sikuliaq LDS text file into string columns `time` and `payload`.

    Lines are ``<log name>\\t<ISO time>\\t<payload>``. Header lines start with
    ``#``. A zero-byte file gives an empty frame.
    """
    try:
        return pd.read_csv(
            file,
            sep="\t",
            comment="#",
            header=None,
            names=["id", "time", "payload"],
            usecols=["time", "payload"],
            dtype=str,
            on_bad_lines="skip",
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame({"time": [], "payload": []}, dtype=str)

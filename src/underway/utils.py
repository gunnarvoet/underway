#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module underway.utils with utility functions"""

import datetime
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker


def concise_date(ax=None, minticks=3, maxticks=10, show_offset=True, **kwargs):
    """
    Better date ticks using matplotlib's ConciseDateFormatter.

    Parameters
    ----------
    ax : axis handle
        Handle to axis (optional).
    minticks : int
        Minimum number of ticks (optional, default 6).
    maxticks : int
        Maximum number of ticks (optional, default 10).
    show_offset : bool, optional
        Show offset string to the right (default True).

    Note
    ----
    Currently only works for x-axis

    See Also
    --------
    matplotlib.mdates.ConciseDateFormatter : For formatting options that
      can be used here.
    """
    if ax is None:
        ax = plt.gca()
    locator = mdates.AutoDateLocator(minticks=minticks, maxticks=maxticks)
    formatter = mdates.ConciseDateFormatter(locator, show_offset=show_offset, **kwargs)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    # remove axis label "time" if present
    if ax.get_xlabel() == "time":
        _ = ax.set_xlabel("")


def quickfig(fs=11, yi=True, w=6, h=4, fgs=None, r=1, c=1, grid=False, **kwargs):
    """
    Quick single pane figure.

    Automatically sets yaxis to be decreasing upwards so
    we can plot against depth.

    Parameters
    ----------
    fs : int, optional
        Fontsize (default 10)
    yi : bool, optional
        Increasing yaxis (default False)
    w : float, optional
        Figure width in inches (default 6)
    h : float, optional
        Figure height in inches (default 4)
    fgs : (float, float)
        Figure size, constructed as (w, h) if not specified here.
    r : int, optional
        Number of rows (default 1)
    c : int, optional
        Number of columns (default 1)
    grid : bool
        Show grid (default False)

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure handle
    ax : matplotlib.axes._subplots.AxesSubplot
        Axis handle
    """
    if fgs is None:
        fgs = (w, h)

    fig, ax = plt.subplots(
        nrows=r,
        ncols=c,
        figsize=fgs,
        constrained_layout=True,
        dpi=75,
        **kwargs,
    )
    if isinstance(ax, np.ndarray):
        [axstyle(axi, fontsize=fs, grid=grid) for axi in ax.flatten()]
    else:
        axstyle(ax, fontsize=fs, grid=grid)
    if yi is False:
        ax.invert_yaxis()
    if r == 1 and c == 1:
        ax.autoscale()

    # some adjustments when using ipympl
    current_backend = mpl.get_backend()
    if current_backend == "module://ipympl.backend_nbagg":
        fig.canvas.header_visible = False
        fig.canvas.toolbar_position = "bottom"
        fig.canvas.resizable = False

    return fig, ax


def determine_header_length(file, comment_char="#"):
    header_length = 1
    with open(file) as f:
        g = f.readline()
        while g[0] == comment_char:
            header_length += 1
            g = f.readline()
    return header_length


def now_utc():
    return np.datetime64(datetime.datetime.now(datetime.UTC).replace(tzinfo=None))


def parse_datetime(dt):
    return np.datetime64(dt.replace(tzinfo=None))


def datetime64_to_str(dt64, unit="D"):
    """Convert numpy datetime64 object or array to str or array of str.

    Parameters
    ----------
    dt64 : np.datetime64 or array-like
        Time in numpy datetime64 format
    unit : str, optional
        Date unit. Defaults to "D".

    Returns
    -------
    str or array of str

    Notes
    -----
    Valid date unit formats are listed at
    https://numpy.org/doc/stable/reference/arrays.datetime.html#arrays-dtypes-dateunits

    """

    return np.datetime_as_string(dt64, unit=unit)


def copy_files(remotedir, local_dir, pattern=None, verbose=False):
    """Copy underway data from ship server."""
    if pattern:
        files = sorted(list(remotedir.glob(pattern)))
    else:
        files = [remotedir]
    for f in files:
        f_local = local_dir.joinpath(f.name)
        if f_local.exists():
            size_remote = f.stat().st_size
            size_local = f_local.stat().st_size
            if size_remote == size_local:
                do_copy = 0
            else:
                do_copy = 1
        else:
            do_copy = 1
        if do_copy:
            # copy2 gave permission errors, switching to copyfile
            # shutil.copy2(f, f_local)
            shutil.copyfile(f, f_local)
            if verbose:
                print(f"copy {f.name}")
    # Unlock files - we are copying file attributes and do not want any locked
    # files, otherwise we are in trouble when syncing again.
    # subprocess.call(["chflags", "-R", "nouchg", f"{local_dir}"])

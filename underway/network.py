#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module underway.network with functions for the local ship network.
"""

from pathlib import Path
import subprocess
import os


def connect_server(server, drive):
    """Connect to a drive on a server. Uses Apple's osascript.

    Parameters
    ----------
    server : str
        Server name or IP address
    drive : str
        Drive name

    Returns
    -------
    output : subprocess returns
    """

    command = 'mount volume "smb://{}/{}"'.format(server, drive)
    output = subprocess.run(["osascript", "-e", command], capture_output=True)
    return output


# == VESSEL SPECIFIC ==

# -> R/V ARMSTRONG
def connect_servers_armstrong():
    """
    Connect to drives data_on_memory and science_share on R/V Armstrong
    """
    vol = Path("/Volumes/")
    server = "10.100.100.30"
    drives = ["data_on_memory", "science_share"]
    drives_local = [vol.joinpath(si) for si in drives]
    # see if drives are connected already
    volg = list(vol.glob("*"))
    conn = {}
    for si in drives:
        sil = vol.joinpath(si)
        conn[si] = False if sil in volg else True

    for k, v in conn.items():
        if v:
            print("connecting to {}...".format(k))
            out = connect_server(server, k)
            if out.returncode == 0:
                print("connected")
            else:
                print(out.stdout.decode("utf-8"))
                print(out.stderr.decode("utf-8"))

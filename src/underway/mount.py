"""Mount SMB shares of the ship server."""

import subprocess
import sys
from pathlib import Path


class MountError(Exception):
    """One or more drives could not be mounted."""


def connect(servers, mount_root):
    """Mount every drive that is not yet present under `mount_root`.

    On macOS the drives are mounted with ``osascript``. On other platforms
    nothing is mounted, and drives that are absent are reported.

    Parameters
    ----------
    servers : dict
        Drive name -> server name.
    mount_root : Path
        Directory under which the drives appear, ``/Volumes`` on macOS.

    Raises
    ------
    MountError
        Lists every drive that is still missing.
    """
    missing = []
    for drive, server in servers.items():
        if (Path(mount_root) / drive).is_dir():
            continue
        if sys.platform != "darwin":
            missing.append(
                f"//{server}/{drive} (mount it at {Path(mount_root) / drive})"
            )
            continue
        command = f'mount volume "smb://{server}/{drive}"'
        result = subprocess.run(
            ["osascript", "-e", command], capture_output=True, text=True, check=False
        )
        if result.returncode != 0:
            missing.append(f"//{server}/{drive} ({result.stderr.strip()})")
    if missing:
        raise MountError("not mounted: " + ", ".join(missing))

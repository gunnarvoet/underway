"""Copy files from a mounted ship share to the local cruise directory."""

import fnmatch
import shutil
import subprocess
import sys


class TransferError(Exception):
    """A transfer failed in part or in full."""


def rsync(src_dir, dst_dir, pattern="*", exclude=(), verbose=False):
    """Sync files matching `pattern` from `src_dir` to `dst_dir` with rsync.

    Parameters
    ----------
    src_dir, dst_dir : Path
        Source and destination directories.
    pattern : str, optional
        Glob relative to `src_dir`, matched at that level only. May contain
        directories, e.g. ``*/contour/*.nc``. The default copies the whole
        tree.
    exclude : tuple of str, optional
        rsync exclude patterns, applied ahead of `pattern`.
    verbose : bool, optional
        Print rsync's file list.

    Raises
    ------
    TransferError
        If rsync exits with a non-zero return code.
    """
    cmd = ["rsync", "-a"]
    if verbose:
        cmd.append("-v")
    for item in exclude:
        cmd += ["--exclude", item]
    if pattern != "*":
        # anchor the pattern at the source root and include only the
        # directories on its path, so rsync does not walk the whole tree
        parts = pattern.split("/")
        for i in range(1, len(parts)):
            cmd += ["--include", "/" + "/".join(parts[:i]) + "/"]
        cmd += ["--include", "/" + pattern, "--exclude", "*"]
    cmd += [f"{src_dir}/", f"{dst_dir}/"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if verbose and result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        raise TransferError(
            f"rsync from {src_dir} failed with code {result.returncode}: "
            f"{result.stderr.strip()}"
        )


def copy(src_dir, dst_dir, pattern="*", exclude=(), verbose=False):
    """Copy files that are missing locally or differ in size.

    For shares where rsync fails. Arguments as for `rsync`, except that the
    default pattern ``"*"`` copies the files directly in `src_dir` and does
    not descend. Files that
    cannot be read are skipped, and one `TransferError` naming them is
    raised after all other files are copied.
    """
    failed = []
    for src in sorted(src_dir.glob(pattern)):
        relative = src.relative_to(src_dir)
        # like rsync, an exclude pattern applies to directories and files
        if any(fnmatch.fnmatch(p, item) for p in relative.parts for item in exclude):
            continue
        dst = dst_dir / relative
        try:
            if not src.is_file():
                continue
            if dst.exists() and dst.stat().st_size == src.stat().st_size:
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            if verbose:
                print(f"copy {src.name}")
        except OSError as err:
            failed.append(f"{src.name} ({err})")
    if sys.platform == "darwin":
        # files copied from the ship share can arrive locked
        subprocess.run(["chflags", "-R", "nouchg", str(dst_dir)], check=False)
    if failed:
        raise TransferError(f"copy from {src_dir} failed for: " + ", ".join(failed))


def make_readonly(directory):
    """Set every file below `directory` to mode 0o400."""
    for file in directory.rglob("*"):
        if file.is_file() and not file.is_symlink():
            file.chmod(0o400)


TRANSFER = {"rsync": rsync, "copy": copy}

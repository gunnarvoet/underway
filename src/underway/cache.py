"""One netCDF product per raw file, re-parsed when the raw file size changes."""

import os

import xarray as xr


def product_path(proc_dir, cruise_id, source, raw, raw_dir=None):
    """Return the product path for a raw file.

    The full raw file name is kept, since some raw names carry the date in
    the suffix (``ins_seapath_position.20251126T0000Z``). With `raw_dir`
    given, directories between `raw_dir` and the file become part of the
    name, so equal file names in two directories give two products.
    """
    name = raw.name
    if raw_dir is not None:
        name = "__".join(raw.relative_to(raw_dir).parts)
    return proc_dir / f"{cruise_id.lower()}_{source}_{name}.nc"


def is_current(product, raw):
    """Return True if `product` was parsed from `raw` at its present size."""
    if not product.exists():
        return False
    with xr.open_dataset(product) as ds:
        return int(ds.attrs.get("raw_size", -1)) == raw.stat().st_size


def write_atomic(ds, path):
    """Write `ds` to a temporary file next to `path`, then rename it into place.

    A failed write leaves an existing file at `path` untouched.
    """
    tmp = path.with_name(path.name + ".tmp")
    try:
        ds.to_netcdf(tmp)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def update(raw, product, parser, force=False):
    """Parse `raw` into `product` unless the product is current.

    Parameters
    ----------
    raw : Path
        Raw file.
    product : Path
        netCDF product for this raw file.
    parser : callable
        Called as ``parser([raw])``, returns an xr.Dataset.
    force : bool, optional
        Parse even when the stored size matches.

    Returns
    -------
    bool
        True if the product was written.
    """
    if not force and is_current(product, raw):
        return False
    # size is read before parsing, so a file that grows during the parse is
    # picked up again on the next call
    size = raw.stat().st_size
    ds = parser([raw])
    ds.attrs["raw_name"] = raw.name
    ds.attrs["raw_size"] = size
    write_atomic(ds, product)
    return True

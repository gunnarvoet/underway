"""Ship servers, drives, and default sources. Data only.

Remote paths are carried over from at-sea use and change between cruises.
Override an entry for one cruise with `Cruise.add_source` under the same name.
In `read_met`, the first source that provides a variable wins, so GPS comes
first.
"""

from dataclasses import dataclass
from functools import partial

from .parsers import armstrong, lds, revelle, seapath, techsas
from .source import Source


@dataclass(frozen=True)
class Ship:
    """A research vessel.

    Parameters
    ----------
    name : str
        Display name.
    servers : dict
        Drive name -> server name.
    sources : tuple of Source
        Default sources.
    """

    name: str
    servers: dict[str, str]
    sources: tuple[Source, ...]


_LDS = "{cruise_id}/lds/raw"
_TECHSAS = "Ship_Systems/Data/TechSAS/NetCDF"

SHIPS = {
    "sikuliaq": Ship(
        name="R/V Sikuliaq",
        servers={
            "CruiseData": "data.sikuliaq.alaska.edu",
            "science": "files.sikuliaq.alaska.edu",
        },
        sources=(
            Source(
                "gps",
                drive="CruiseData",
                remote=f"{_LDS}/ins_seapath_position",
                pattern="ins_seapath_position.*",
                parser=seapath.read,
                met=True,
            ),
            Source(
                "tsg",
                drive="CruiseData",
                remote=f"{_LDS}/tsg_sbe45_fwd",
                pattern="tsg_sbe45_fwd.*",
                parser=partial(lds.read, stream="tsg"),
                met=True,
            ),
            Source(
                "wind",
                drive="CruiseData",
                remote=f"{_LDS}/wind_gill_fwdmast_true",
                pattern="wind_gill_fwdmast_true.*",
                parser=partial(lds.read, stream="wind"),
                met=True,
            ),
            Source(
                "air",
                drive="CruiseData",
                remote=f"{_LDS}/met_met4a_fwdmast",
                pattern="met_met4a_fwdmast.*",
                parser=partial(lds.read, stream="air"),
                met=True,
            ),
            Source(
                "sadcp",
                drive="CruiseData",
                remote="{cruise_id}/adcp/raw/{cruise_id}/proc",
                pattern="*/contour/*.nc",
                transfer="copy",
            ),
            Source(
                "ctd",
                drive="CruiseData",
                remote="{cruise_id}/ctd/raw",
                readonly=True,
            ),
        ),
    ),
    # TechSAS and ADCP files use transfer="copy". rsync was replaced by a
    # size-compare copy for these in June 2021 (commit cd20c84), files arrive
    # locked (hence chflags in transfer.copy), and shutil.copy2 gave
    # permission errors (hence copyfile). CTD and LADCP stayed on rsync.
    "discovery": Ship(
        name="RRS Discovery",
        servers={
            "current_cruise": "dynetapp.discovery.ad.noc.ac.uk",
            "science_public": "dynetapp.discovery.ad.noc.ac.uk",
        },
        sources=(
            Source(
                "gps",
                drive="current_cruise",
                remote=f"{_TECHSAS}/GPS",
                pattern="*position-POSMV_GPS.gps",
                transfer="copy",
                parser=techsas.read,
                cache=False,
                met=True,
            ),
            Source(
                "surfmet",
                drive="current_cruise",
                remote=f"{_TECHSAS}/SURFMETV3",
                pattern="*MET-SURFMET.SURFMETv3",
                transfer="copy",
                parser=techsas.read,
                cache=False,
                met=True,
            ),
            Source(
                "light",
                drive="current_cruise",
                remote=f"{_TECHSAS}/SURFMETV3",
                pattern="*Light-SURFMET.SURFMETv3",
                transfer="copy",
                parser=techsas.read,
                cache=False,
                met=True,
            ),
            Source(
                "surf",
                drive="current_cruise",
                remote=f"{_TECHSAS}/SURFMETV3",
                pattern="*Surf-SURFMET.SURFMETv3",
                transfer="copy",
                parser=techsas.read,
                cache=False,
            ),
            Source(
                "tsg",
                drive="current_cruise",
                remote=f"{_TECHSAS}/TSG",
                pattern="*SBE45-SBE45.TSG",
                transfer="copy",
                parser=techsas.read,
                cache=False,
                met=True,
            ),
            Source(
                "sadcp",
                drive="current_cruise",
                remote="Ship_Systems/Data/Acoustics/ADCP/proc",
                pattern="os*nb/contour/*",
                transfer="copy",
            ),
            Source(
                "ctd",
                drive="current_cruise",
                remote="Sensors_and_Moorings/CTD/Data/Raw",
            ),
            Source(
                "ladcp",
                drive="current_cruise",
                remote="Sensors_and_Moorings/LADCP/Data",
            ),
        ),
    ),
    "armstrong": Ship(
        name="R/V Neil Armstrong",
        servers={
            "data_on_memory": "10.100.100.30",
            "science_share": "10.100.100.30",
        },
        sources=(
            Source(
                "gps",
                drive="data_on_memory",
                remote="underway/raw",
                pattern="*.CNAV_3050",
                transfer="copy",
                parser=armstrong.read_gps,
                met=True,
            ),
            Source(
                "met",
                drive="data_on_memory",
                remote="underway/proc",
                pattern="AR[0-9]*.csv",
                parser=armstrong.read_met,
                met=True,
            ),
            Source(
                "sadcp",
                drive="data_on_memory",
                remote="adcp/proc",
                pattern="*/contour/*.nc",
            ),
            Source("ctd", drive="data_on_memory", remote="ctd", readonly=True),
        ),
    ),
    "revelle": Ship(
        name="R/V Roger Revelle",
        servers={
            "cruise": "rr-sci-filesvr.ucsd.edu",
            "science_party_share": "rr-sci-filesvr.ucsd.edu",
        },
        sources=(
            Source(
                "met",
                drive="cruise",
                remote="{cruise_id}/metacq/data",
                pattern="*.MET",
                parser=revelle.read,
                met=True,
            ),
            Source(
                "sadcp",
                drive="cruise",
                remote="{cruise_id}/adcp_uhdas/{cruise_id}/proc",
                pattern="*/contour/*.nc",
            ),
        ),
    ),
}

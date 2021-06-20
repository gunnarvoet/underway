#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module underway.io with in/out functions
"""

from pathlib import Path
import xarray as xr
import pandas as pd
import subprocess
import os
import shutil

from . import network


class Underway:
    def __init__(self, ship, localdir, atsea=True):
        """Generate underway object.

        Parameters
        ----------
        ship : str ["armstrong", "discovery"]
            Select ship.
        localdir : PosixPath
            Path to local data directory. Will have subdirectories "met", "ctd" and "sadcp".
        """

        self.ship = ship
        self.atsea = atsea
        self.localdir = localdir
        self.local_met = localdir / "met"
        self.local_ctd = localdir / "ctd"
        self.local_sadcp = localdir / "sadcp"
        if ship == "armstrong":
            self.remote_met = Path("/Volumes/data_on_memory/underway/proc")
            self.met_pattern = "AR20*.csv"
            self.read_met_file = read_met_file_armstrong
            self.read_all_met = read_all_met_armstrong
            self.remote_ctd = Path("/Volumes/data_on_memory/ctd/")
            self.connect = network.connect_servers_armstrong
            self.position = network.get_position_armstrong
            # sync shipboard adcp data
            sync_sadcp_armstrong(self.local_sadcp)
            # sync underway met data
            rsync_underway_data(
                self.remote_met, self.local_met, self.met_pattern
            )
        elif ship == "discovery":
            self.connect = network.connect_servers_discovery
            sync_underway_disco(self.local_met)
            sync_sadcp_disco(self.local_sadcp)
            sync_ctd_disco(self.local_ctd)
            self.read_all_met = read_all_met_discovery
        if atsea:
            # connect to ship servers if at sea
            self.connect()
        # read met data into data structure
        self.met = self.read_all_met(self.local_met)

    def sync_ctd_data(self):
        """Sync CTD data from ship server."""
        sync_ctd_data(self.remote_ctd, self.local_ctd)


def rsync_underway_data(remotedir, local_met, pattern):
    """Sync underway data from ship server using rsync."""
    print("syncing data from ship server")
    for f in sorted(list(remotedir.glob(pattern))):
        subprocess.call(["rsync", "-avz", f, local_met])


def copy_underway_data(remotedir, local_met, pattern=None):
    """Copy underway data from ship server."""
    print("syncing data from ship server")
    if pattern:
        files = sorted(list(remotedir.glob(pattern)))
    else:
        files = [remotedir]
    for f in files:
        # subprocess.call(["rsync", "-avz", f, local_met])
        f_local = local_met.joinpath(f.name)
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
            shutil.copy2(f, local_met)
            print(f"copy {f.name}")
    # Unlock files - we are copying file attributes and do not want any locked
    # files, otherwise we are in trouble when syncing again.
    subprocess.call(["chflags", "-R", "nouchg", f"{local_met}"])


def sync_ctd_data(remote_ctd, local_ctd):
    """Sync CTD data from ship server."""
    print("syncing ctd data from ship server")
    subprocess.call(["rsync", "-avzs", remote_ctd.as_posix() + "/", local_ctd])
    ctd_files_chmod(local_ctd)


def ctd_files_chmod(local_ctd):
    ctd_suffix = [
        "dsa",
        "psa",
        "xmlcon",
        "XMLCON",
        "hdr",
        "bl",
        "hex",
        "ros",
        "cnv",
        "asc",
        "btl",
    ]
    for p in ctd_suffix:
        for f in local_ctd.glob("**/*." + p):
            f.chmod(400)


def combine_netcdf(file_list):
    readall = []
    for f in sorted(file_list):
        readall.append(xr.open_dataset(f))
    b = xr.concat(readall, dim="time")
    return b


# === VESSEL SPECIFIC ===

# -> RRS DISCOVERY
def sync_underway_disco(local_met):
    techsas = Path("/Volumes/current_cruise/Ship_Systems/Data/TechSAS/NetCDF/")
    sources = dict(
        met=dict(dir="SURFMETV3", pattern="*.SURFMETv3"),
        gps=dict(dir="GPS", pattern="*position-POSMV_GPS.gps"),
        tsg=dict(dir="TSG", pattern="*SBE45-SBE45.TSG"),
    )
    for src, p in sources.items():
        remote = techsas / p["dir"]
        copy_underway_data(remote, local_met, p["pattern"])


def sync_sadcp_disco(local_sadcp):
    """currently only running bb"""
    adcp = Path(
        "/Volumes/current_cruise/Ship_Systems/Data/Acoustics/ADCP/UHDAS/DY132/proc/"
    )
    sources = ["os75bb", "os150bb"]
    for src in sources:
        file = src + ".nc"
        remote = adcp / src / "contour" / file
        try:
            copy_underway_data(remote, local_sadcp)
        except PermissionError:
            pass


def sync_ctd_disco(local_ctd):
    remote_ctd = Path(
        "/Volumes/current_cruise/Sensors_and_Moorings/DY132/CTD/Data/CTD Raw Data"
    )
    # sync_ctd_data(remote_ctd, local_ctd)
    rsync_underway_data(remote_ctd, local_ctd, "*")


def read_all_met_discovery(local_met):
    sources = dict(
        light="*Light-SURFMET.SURFMETv3",
        met="*MET-SURFMET.SURFMETv3",
        surf="*Surf-SURFMET.SURFMETv3",
        gps="*position-POSMV_GPS.gps",
        tsg="*SBE45-SBE45.TSG",
    )
    out = dict()
    for key, pattern in sources.items():
        files = sorted(local_met.glob(pattern))
        out[key] = combine_netcdf(files)

    # combine met, light, surf (they come with the same time stamps)
    # interpolate gps, tsg to these
    met_all = xr.merge(
        [
            out["light"],
            out["met"],
            out["surf"],
            out["gps"].interp_like(out["met"]),
            out["tsg"].interp_like(out["met"]),
        ]
    )
    return met_all


# -> R/V ARMSTRONG
def read_met_file_armstrong(file):
    met = pd.read_csv(
        file,
        skiprows=1,
        delimiter=", ",
        engine="python",
        na_values=["NODATA", "NAN"],
    )
    met.index = pd.to_datetime(met.DATE_GMT + " " + met.TIME_GMT)
    a = met.to_xarray()
    a = a.rename_dims({"index": "time"})
    a = a.rename_vars({"index": "time"})
    a = rename_variables_armstrong(a)
    return a


def read_all_met_armstrong(local_met):
    metall = local_met.glob("*.csv")
    readall = []
    for f in sorted(metall):
        readall.append(read_met_file_armstrong(f))
    b = xr.concat(readall, dim="time")
    b = add_attributes_armstrong(b)
    return b


def rename_variables_armstrong(a):
    names = dict(
        Dec_LON="lon",
        Dec_LAT="lat",
        SPD="spd",
        HDT="hdg",
        SOG="sog",
        COG="cog",
        WXTS_Ta="air_temperature_stbd",
        WXTP_Ta="air_temperature_port",
        WXTS_Pa="air_pressure_stbd",
        WXTP_Pa="air_pressure_port",
        WXTS_Ri="rain_intensity_stbd",
        WXTP_Ri="rain_intensity_port",
        WXTS_Rc="rain_accumulation_stbd",
        WXTP_Rc="rain_accumulation_port",
        WXTS_Dm="relative_wind_direction_stbd",
        WXTP_Dm="relative_wind_direction_port",
        WXTS_Sm="relative_wind_speed_stbd",
        WXTP_Sm="relative_wind_speed_port",
        WXTS_Ua="relative_humidity_stbd",
        WXTP_Ua="relative_humidity_port",
        WXTS_TS="true_wind_speed_stbd",
        WXTP_TS="true_wind_speed_port",
        WXTS_TD="true_wind_direction_stbd",
        WXTP_TD="true_wind_direction_port",
        BAROM_S="barometric_pressure_stbd",
        BAROM_P="barometric_pressure_port",
        RAD_SW="shortwave_radiation",
        RAD_LW="longwave_radiation",
        PAR="photosynthetically_active_radiation",
        SBE45S="sea_surface_salinity",
        SBE48T="sea_surface_temperature",
        FLR="fluorometer",
        SSVdslog="sea_surface_sound_velocity",
        Depth12="water_depth_12khz",
        Depth35="water_depth_35khz",
        EM122="water_depth_multibeam",
    )
    a = a.rename_vars(names)
    return a


def add_attributes_armstrong(a):
    a.spd.attrs = dict(long_name="ship speed")
    a.hdg.attrs = dict(long_name="ship heading")
    a.cog.attrs = dict(long_name="ship course over ground")
    a.sog.attrs = dict(long_name="ship speed over ground")
    a.air_temperature_stbd.attrs = dict(
        long_name="air temperature starboard", units="°C"
    )
    a.air_temperature_port.attrs = dict(
        long_name="air temperature port", units="°C"
    )
    a.air_pressure_stbd.attrs = dict(long_name="air pressure stbd", units="hPa")
    a.air_pressure_port.attrs = dict(long_name="air pressure port", units="hPa")
    a.rain_intensity_stbd.attrs = dict(
        long_name="rain intensity stbd", units="mm/h"
    )
    a.rain_intensity_port.attrs = dict(
        long_name="rain intensity port", units="mm/h"
    )
    a.rain_accumulation_stbd.attrs = dict(
        long_name="rain accumulation stbd", units="mm/h"
    )
    a.rain_accumulation_port.attrs = dict(
        long_name="rain accumulation port", units="mm/h"
    )
    a.relative_wind_direction_stbd.attrs = dict(
        long_name="relative wind direction stbd", units="deg"
    )
    a.relative_wind_direction_port.attrs = dict(
        long_name="relative wind direction port", units="deg"
    )
    a.relative_wind_speed_stbd.attrs = dict(
        long_name="relative wind speed stbd", units="m/s"
    )
    a.relative_wind_speed_port.attrs = dict(
        long_name="relative wind speed port", units="m/s"
    )
    a.relative_humidity_stbd.attrs = dict(
        long_name="relative humidity stbd", units="%"
    )
    a.relative_humidity_port.attrs = dict(
        long_name="relative humidity port", units="%"
    )
    a.true_wind_speed_stbd.attrs = dict(
        long_name="true wind speed stbd", units="m/s"
    )
    a.true_wind_speed_port.attrs = dict(
        long_name="true wind speed port", units="m/s"
    )
    a.true_wind_direction_stbd.attrs = dict(
        long_name="true wind direction stbd", units="deg"
    )
    a.true_wind_direction_port.attrs = dict(
        long_name="true wind direction port", units="deg"
    )
    a.barometric_pressure_stbd.attrs = dict(
        long_name="barometric pressure stbd", units="hPa"
    )
    a.barometric_pressure_port.attrs = dict(
        long_name="barometric pressure port", units="hPa"
    )
    a.shortwave_radiation.attrs = dict(
        long_name="shortwave radiation", units="W/m$^2$"
    )
    a.longwave_radiation.attrs = dict(
        long_name="longwave radiation", units="W/m$^2$"
    )
    a.photosynthetically_active_radiation.attrs = dict(
        long_name="photosynthetically active radiation", units="uE/m$^2$/s"
    )
    a.sea_surface_salinity.attrs = dict(
        long_name="sea surface salinity (5m)", units="psu"
    )
    a.sea_surface_temperature.attrs = dict(
        long_name="sea surface temperature (5m)", units="°C"
    )
    a.fluorometer.attrs = dict(long_name="fluorometer", units="mV")
    a.sea_surface_sound_velocity.attrs = dict(
        long_name="sea surface sound velocity", units="m/s"
    )
    a.water_depth_12khz.attrs = dict(long_name="12kHz water depth", units="m")
    a.water_depth_35khz.attrs = dict(long_name="35kHz water depth", units="m")
    a.water_depth_multibeam.attrs = dict(
        long_name="multibeam water depth", units="m"
    )
    return a


def sync_sadcp_armstrong(local_sadcp):
    """Sync R/V Armstrong shipboard ADCP data

    Parameters
    ----------
    local_sadcp :


    """

    wh300path = Path("/Volumes/data_on_memory/adcp/proc/wh300/contour/wh300.nc")
    os150nbpath = Path(
        "/Volumes/data_on_memory/adcp/proc/os150nb/contour/os150nb.nc"
    )
    os38nbpath = Path(
        "/Volumes/data_on_memory/adcp/proc/os38nb/contour/os38nb.nc"
    )
    print("syncing sadcp data from ship server")
    subprocess.call(["rsync", "-avz", wh300path.as_posix(), local_sadcp])
    subprocess.call(["rsync", "-avz", os150nbpath.as_posix(), local_sadcp])
    subprocess.call(["rsync", "-avz", os38nbpath.as_posix(), local_sadcp])

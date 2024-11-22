#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module underway.ship defining classes for each research vessel.

This is a rework of underway.io which I am keeping for compatibility with older code. Here, each vessel defines its own class that is based on an abstract base class.
"""

import pathlib
import socket
from itertools import islice
from pathlib import Path
from abc import ABC, abstractmethod
import datetime
import xarray as xr
import pandas as pd
import numpy as np
import subprocess
import shutil
from pynmeagps import NMEAReader
import pynmea2

import underway


class Underway(ABC):
    def __init__(self, local_data, atsea=True, cruise_id=None, bathy=None):
        """Generate underway object.

        Parameters
        ----------
        local_data : PosixPath
            Path to local data directory. Will have subdirectories "met", "ctd"
            and "sadcp".
        atsea : bool, optional
            Set to True to sync data. Defaults to True.
        cruise_id : str, optional
            Cruise ID.
        bathy : PosixPath, optional
            Path to bathymetry file (in netcdf format) for quick plots.
        """
        self.atsea = atsea
        self.cruise_id = cruise_id
        # create directories for met / ctd / sadcp / ladcp / gps
        self.create_base_dirs(local_data)
        # better to do this when needed
        # if self.atsea:
        #     self.sync_data()
        #     self.met = self.read_all_met()

    def create_base_dirs(self, localdir):
        if type(localdir) is not pathlib.PosixPath:
            localdir = Path(localdir)
        self.localdir = localdir
        for dir in ["met", "ctd", "sadcp", "ladcp", "gps"]:
            create_dir = self.localdir / dir
            setattr(self, f"local_{dir}", create_dir)
            create_dir.mkdir(exist_ok=True)
        # make raw and proc directories for ctd, ladcp, gps (met & sadcp
        # already processed)
        for dir in ["ctd", "ladcp", "gps"]:
            for rp in ["raw", "proc"]:
                create_dir = self.localdir / dir / rp
                setattr(self, f"local_{dir}_{rp}", create_dir)
                create_dir.mkdir(exist_ok=True)

    @abstractmethod
    def connect(self):
        pass

    @abstractmethod
    def read_all_met(self):
        pass

    @abstractmethod
    def sync_data(self):
        pass

    def quick_plot_wind_speed(self):
        fig, ax = underway.utils.quickfig(grid=True)
        self.met[self.metvar_wind_speed].plot(ax=ax)
        underway.utils.concise_date(ax)

    def save_met(self):
        savename = "met.nc" if self.cruise_id is None else f"met_{self.cruise_id}.nc"
        metnc = self.local_met.joinpath(savename)
        self.met.to_netcdf(metnc)


class Sikuliaq(Underway):
    def __init__(self, local_data, atsea=True, cruise_id=None, bathy=None):
        """Generate Sikuliaq underway object.

        Parameters
        ----------
        local_data : PosixPath
            Path to local data directory. Will have subdirectories "met", "ctd"
            and "sadcp".
        atsea : bool, optional
            Set to True to sync data. Defaults to True.
        cruise_id : str, optional
            Cruise ID.

        Notes
        -----
        The Sikuliaq doesn't have a nice met system and we have to gather data from various sources.
        """
        super().__init__(local_data, atsea, cruise_id)
        self.ship = "sikuliaq"
        self.set_dataset_ids()
        self.set_paths()

    def set_dataset_ids(self):
        ids = dict(
            wind="wind_gill_fwdmast_true",
            air="met_met4a_fwdmast",
            rain="rain_org815ds",
            tsg="tsg_sbe45_fwd",
            heat="",
            wave="wave_wamos",
        )
        self.ids = ids

    def set_paths(self):
        self.remote_lds = Path("/Volumes/data/SKQ202417S/lds/raw")
        self.remote_gps = self.remote_lds.joinpath("ins_seapath_position")
        # set lds paths
        for name, id in self.ids.items():
            setattr(self, f"remote_{name}", self.remote_lds.joinpath(id))
        self.remote_adcp = Path("/Volumes/data/SKQ202417S/adcp/raw/SKQ202417S/proc")
        self.remote_ww = Path("/Volumes/sci/shipside/SKQ202417S/WireWalkerEmails/nc")
        self.remote_ctd = Path(
            f"/Volumes/data/{self.cruise_id}/ctd/raw/{self.cruise_id}"
        )
        # self.remote_ladcp = Path("/Volumes/science_share/LADCP/")

        # create directory non-standard directories
        dirs = list(self.ids.keys()) + ["ww"]
        for dir in dirs:
            create_dir = self.localdir / dir
            setattr(self, f"local_{dir}", create_dir)
            create_dir.mkdir(exist_ok=True)
            for rp in ["raw", "proc"]:
                create_dir = self.localdir / dir / rp
                setattr(self, f"local_{dir}_{rp}", create_dir)
                create_dir.mkdir(exist_ok=True)

    def connect(self):
        underway.network.connect_servers_sikuliaq()

    def sync_data(self, verbose=False):
        cruise_id = self.cruise_id.upper()
        self.connect()
        if verbose:
            print("-----")
            print("syncing gps data...")
        self.sync_gps_data(verbose=verbose)
        if verbose:
            print("-----")
            print("syncing tsg data...")
        self.sync_tsg_data(verbose=verbose)
        if verbose:
            print("-----")
            print("syncing atmospheric data...")
        self.sync_air_data(verbose=verbose)
        if verbose:
            print("-----")
            print("syncing wind data...")
        self.sync_wind_data(verbose=verbose)
        # if verbose:
        #     print("-----")
        #     print("syncing wave data...")
        # self.sync_wave_data(verbose=verbose)
        if verbose:
            print("-----")
            print("syncing adcp data...")
        self.sync_adcp_data(verbose=verbose)
        if verbose:
            print("-----")
            print("syncing ctd data...")
        self.sync_ctd_data(verbose=verbose)
        # if self.cruise_id == "ar73":
        #     print("syncing ladcp data...")
        #     self.ar73_sync_ladcp()

    def sync_wave_data(self, verbose=False):
        """Sync Wamos data from ship server."""
        rsync_opts = "-av" if verbose else "-a"
        subprocess.call(
            [
                "rsync",
                rsync_opts,
                self.remote_wave.as_posix() + "/",
                self.local_wave_raw.as_posix() + "/",
            ]
        )

    def sync_tsg_data(self, verbose=False):
        """Sync TSG data from ship server."""
        rsync_opts = "-av" if verbose else "-a"
        subprocess.call(
            [
                "rsync",
                rsync_opts,
                self.remote_tsg.as_posix() + "/",
                self.local_tsg_raw.as_posix() + "/",
            ]
        )

    def sync_air_data(self, verbose=False):
        """Sync atmospheric data from ship server."""
        rsync_opts = "-av" if verbose else "-a"
        subprocess.call(
            [
                "rsync",
                rsync_opts,
                "--exclude",
                "met_met4a_fwdmast.20241102T2051Z",
                "--exclude",
                "met_met4a_fwdmast.20241103T0000Z",
                self.remote_air.as_posix() + "/",
                self.local_air_raw.as_posix() + "/",
            ]
        )

    def sync_wind_data(self, verbose=False):
        """Sync wind data from ship server."""
        rsync_opts = "-av" if verbose else "-a"
        subprocess.call(
            [
                "rsync",
                rsync_opts,
                "--exclude",
                "wind_gill_fwdmast_true.20241102T2051Z",
                "--exclude",
                "wind_gill_fwdmast_true.20241103T0000Z",
                self.remote_wind.as_posix() + "/",
                self.local_wind_raw.as_posix() + "/",
            ]
        )

    def sync_gps_data(self, verbose=False):
        """Sync GPS data from ship server."""
        rsync_opts = "-av" if verbose else "-a"
        subprocess.call(
            [
                "rsync",
                rsync_opts,
                self.remote_gps.as_posix() + "/",
                self.local_gps_raw.as_posix() + "/",
            ]
        )

    def sync_adcp_data(self, verbose=False):
        """Sync sADCP data from ship server."""
        adcps = ["wh300", "os38nb", "os75nb"]
        for adcp in adcps:
            copy_underway_data(
                self.remote_adcp.joinpath(f"{adcp}/contour"),
                self.local_sadcp,
                pattern=f"{adcp}.nc",
                verbose=verbose,
            )

    def sync_ctd_data(self, verbose=False):
        """Sync CTD data from ship server."""
        _sync_ctd_data(self.remote_ctd, self.local_ctd_raw, verbose=verbose)

    def sync_ww_loc_data(self, verbose=False):
        """Sync WW location data from ship server."""
        copy_underway_data(
            self.remote_ww,
            self.local_ww_proc,
            pattern=f"*.nc",
            verbose=verbose,
        )

    def read_all_met(self, start_time=None):
        pass

    def read_gps_data(self):
        self._proc_gps_files()
        gps = xr.open_mfdataset(self.local_gps_proc.as_posix() + "/*.nc")
        return gps

    def _parse_gps_data_file(self, file, skiplines=0):
        data = []
        found_zda = False
        with open(file, "r") as f:
            # Using islice allows to quickly skip lines in a file.
            # Set stop to something large to reach the end of the file.
            for line in islice(f, skiplines, skiplines + 24 * 3600 * 8 * 10):
                # skip comments
                if line[0] == "#":
                    continue
                # get rid of the identifier and time tag
                line = line.split("\t")[2]
                if line.startswith("$GPZDA"):
                    # add previous data block to list and start a new packet
                    if found_zda:
                        data.append(packet)
                    found_zda = True
                    packet = dict()
                try:
                    msg = pynmea2.parse(line)
                    if isinstance(msg, pynmea2.types.talker.ZDA):
                        dt = msg.datetime
                        packet["time"] = parse_datetime(msg.datetime)
                    elif isinstance(msg, pynmea2.types.talker.GGA):
                        if dt.time() == msg.timestamp.replace(tzinfo=None):
                            packet["lat"] = msg.latitude
                            packet["lon"] = msg.longitude
                            # packet["antenna_altitude"] = msg.altitude
                            # packet["gps_quality"] = msg.gps_qual
                            # packet["gps_num_sats"] = msg.num_sats
                    elif isinstance(msg, pynmea2.types.talker.VTG):
                        packet["sog"] = float(msg.spd_over_grnd_kts)
                        packet["cog"] = msg.true_track
                    elif isinstance(msg, pynmea2.types.talker.HDT):
                        packet["heading"] = float(msg.heading)
                    elif isinstance(msg, pynmea2.types.proprietary.sxn.SXN23):
                        packet["rol"] = float(msg.roll)
                        packet["pitch"] = float(msg.pitch)
                        packet["heave"] = float(msg.heave)
                except:
                    continue

        df = pd.DataFrame(data)
        if not df.empty:
            df = df.set_index("time")
            ds = df.to_xarray()
            ds = ds.dropna("time")
            return ds
        else:
            return None

    def _proc_gps_files(self):
        utc_today = self._utc_yyyymmdd()
        gps_raw_files = sorted(self.local_gps_raw.glob("ins_seapath*"))
        for file in gps_raw_files:
            yyyymmdd = file.name.split(".")[1].split("T")[0]
            nc_name = f"skq202417s_gps_{yyyymmdd}.nc"
            nc_path = self.local_gps_proc.joinpath(nc_name)
            if not nc_path.exists():
                print(f"parsing {yyyymmdd}")
                ds = self._parse_gps_data_file(file)
                ds.to_netcdf(nc_path)
            elif utc_today == yyyymmdd:
                gps = xr.open_dataset(nc_path)
                header_length = underway.utils.determine_header_length(file)
                existing_time_steps = len(gps.time)
                # only update file if there is at least one minute of new data
                if underway.utils.now_utc() - gps.time[-1].data > np.timedelta64(
                    1, "m"
                ):
                    skiplines = -1 + header_length + existing_time_steps * 8
                    ds = self._parse_gps_data_file(file, skiplines=skiplines)
                    if ds is not None:
                        new = xr.concat([gps, ds], dim="time")
                        nc_path.unlink()
                        new.to_netcdf(nc_path)

    def read_ship_data(self, name):
        self._proc_lds_files(self.ids[name], name)
        proc_dir = getattr(self, f"local_{name}_proc")
        data = xr.open_mfdataset(proc_dir.as_posix() + "/*.nc")
        return data

    def read_air_data(self):
        return self.read_ship_data("air")

    def read_wind_data(self):
        return self.read_ship_data("wind")

    def read_tsg_data(self):
        return self.read_ship_data("tsg")

    def _proc_lds_files(self, id, name):
        utc_today = self._utc_yyyymmdd()
        raw_files = sorted(self.localdir.joinpath(name).joinpath("raw").glob(f"{id}*Z"))
        for file in raw_files:
            yyyymmdd = file.name.split(".")[1].split("T")[0]
            nc_name = f"{self.cruise_id.lower()}_{name}_{yyyymmdd}.nc"
            nc_path = self.localdir.joinpath(name).joinpath("proc").joinpath(nc_name)
            if not nc_path.exists():
                print(f"parsing {name} {yyyymmdd}")
                ds = self._parse_lds_file(file, name)
                ds.to_netcdf(nc_path)
            elif utc_today == yyyymmdd:
                tsg = xr.open_dataset(nc_path)
                existing_time_steps = len(tsg.time)
                try:
                    ds = self._parse_lds_file(
                        file, name, skip_extra_rows=existing_time_steps
                    )
                    if ds is not None:
                        ds = xr.concat([tsg, ds], dim="time")
                except:
                    ds = tsg
                nc_path.unlink()
                ds.to_netcdf(nc_path)

    def _parse_lds_file(self, file, name, skip_extra_rows=0):
        header_length = underway.utils.determine_header_length(file)
        match name:
            case "tsg":
                usecols = [1, 2, 3, 4, 5]
                var_names = ["time", "t", "c", "s", "sspd"]
                long_names = ["temperature", "conductivity", "salinity", "sound speed"]
                units = ["°C", "S/m", "psu", "m/s"]
            case "wind":
                usecols = [1, 3, 7, 9]
                var_names = ["time", "wind_dir", "wind_spd_kn", "wind_spd_ms"]
                long_names = ["wind direction", "wind speed", "wind speed"]
                units = ["deg", "kn", "m/s"]
            case "air":
                usecols = [1, 4, 8, 11]
                var_names = [
                    "time",
                    "barometric_pressure",
                    "air_temperature",
                    "relative_humidity",
                ]
                long_names = [
                    "barometric pressure",
                    "air temperature",
                    "relative humidity",
                ]
                units = ["bar", "°C", "%"]
            case "wave":
                usecols = [1, 3, 4, 5, 6, 7, 8, 9]
                var_names = [
                    "time",
                    "significant_wave_height",
                    "mean_wave_period",
                    "peak_wave_direction",
                    "peak_wave_period",
                    "peak_wave_length",
                    "wamos_surface_current_direction",
                    "wamos_surface_current_speed",
                ]
                long_names = [vname.replace("_", " ") for vname in var_names[1:]]
                units = ["m", "s", "deg (coming from)", "s", "m", "deg", "m/s"]
            case "heat":
                usecols = []
                var_names = []
                long_names = []
                units = []
            case "rain":
                usecols = [1, 3, 4]
                var_names = ["time", "precipitation_rate", "precipitation_accumulation"]
                long_names = ["precipitation rate", "precipitation accumulation"]
                units = ["mm/hr", "mm"]

        dtype = {name: np.float64 for name in var_names[1:]}
        df = pd.read_csv(
            file,
            skiprows=header_length + skip_extra_rows,
            sep=r"[,\t]",
            usecols=usecols,
            names=var_names,
            index_col=0,
            parse_dates=True,
            engine="python",
            dtype=dtype,
        )
        df.index = df.index.tz_convert(None)
        ds = df.to_xarray()
        for name, long_name, unit in zip(var_names[1:], long_names, units):
            ds[name].attrs = dict(long_name=long_name, units=unit)
        return ds

    def _utc_yyyymmdd(self):
        utcnow = underway.utils.now_utc()
        utcstr = underway.utils.datetime64_to_str(utcnow)
        y, m, d = utcstr.split("-")
        yyyymmdd = y + m + d
        return yyyymmdd


class Revelle(Underway):
    # Save the two-digit year as a constant to use later on
    _YY = datetime.datetime.utcnow().strftime("%y")

    def __init__(self, local_data, atsea=True, cruise_id=None, bathy=None):
        """Generate Revelle underway object.

        Parameters
        ----------
        local_data : PosixPath
            Path to local data directory. Will have subdirectories "met", "ctd"
            and "sadcp".
        atsea : bool, optional
            Set to True to sync data. Defaults to True.
        cruise_id : str, optional
            Cruise ID.
        """
        super().__init__(local_data, atsea, cruise_id)
        self.ship = "revelle"
        # self.parse_raw_gps()

    def set_quickplot_variables(self):
        self.metvar_wind_speed = "true_wind_speed_stbd"

    def connect(self):
        underway.network.connect_servers_revelle()

    def position(self):
        # return underway.network.get_position_armstrong()
        pass

    def sync_data(self):
        cruise_id = self.cruise_id.upper()
        # self.remote_ctd = Path("/Volumes/cruise/ctd/")
        # self.remote_ladcp = Path("/Volumes/science_share/LADCP/")
        # set the pattern of the files we are looking for
        # self.gps_pattern = f"AR_GPS10_{self._YY}*.csv"

        self.connect()
        self.sync_sadcp()
        # print("syncing met data...")
        # rsync_underway_data(self.remote_met, self.local_met, self.met_pattern)
        self.sync_met()
        # print("syncing gps data...")
        # rsync_underway_data(self.remote_met, self.local_gps, self.gps_pattern)
        # print("syncing gps raw data...")
        # self.sync_raw_gps()
        # print("syncing ctd data...")
        # self.sync_ctd_data()
        # if self.cruise_id == "ar73":
        #     print("syncing ladcp data...")
        #     self.ar73_sync_ladcp()

    def sync_ctd_data(self, verbose=False):
        """Sync CTD data from ship server."""
        # _sync_ctd_data(self.remote_ctd, self.local_ctd_raw, verbose=verbose)
        pass

    def sync_met(self):
        self.remote_met = Path(f"/Volumes/cruise/{self.cruise_id}/metacq/data")
        # set the pattern of the files we are looking for
        self.met_pattern = "*.MET"
        print("syncing met data...")
        rsync_underway_data(self.remote_met, self.local_met, self.met_pattern)

    def sync_sadcp(self, verbose=False):
        """Sync R/V Revelle shipboard ADCP data

        Parameters
        ----------
        verbose : bool, optional
            Show rsync output.
        """
        # wh300path = Path("/Volumes/data_on_memory/adcp/proc/wh300/contour/wh300.nc")
        cr_id = self.cruise_id.upper()
        os150nbpath = Path(
            f"/Volumes/cruise/{cr_id}/adcp_uhdas/{cr_id}/proc/os150nb/contour/os150nb.nc"
        )
        print("syncing sadcp data...")
        # os38nbpath = Path("/Volumes/data_on_memory/adcp/proc/os38nb/contour/os38nb.nc")
        rsync_opts = "-av" if verbose else "-a"
        # subprocess.call(["rsync", rsync_opts, wh300path.as_posix(), self.local_sadcp])
        subprocess.call(["rsync", rsync_opts, os150nbpath.as_posix(), self.local_sadcp])
        # subprocess.call(["rsync", rsync_opts, os38nbpath.as_posix(), self.local_sadcp])

    def _read_met_date(self, file):
        with open(file) as f:
            # date is on the second line
            f.readline()
            date = f.readline()
        dd = date.split(" ")
        date = dd[2]
        return date

    def _parse_met_time_stamp(self, date, time):
        tstr = f"{time:06d}"
        h = tstr[:2]
        m = tstr[2:4]
        s = tstr[4:6]
        return gv.time.str_to_datetime64(date + f" {h}:{m}:{s}")

    def _read_met_file(self, file):
        date = self._read_met_date(file)
        met = pd.read_csv(
            file,
            skiprows=3,
            sep="\\s+",
            engine="python",
        )
        time = met["#Time"]
        dt = np.array([self._parse_met_time_stamp(date, ti) for ti in time])
        met["time"] = dt
        met = met.set_index("time")
        met = met.drop("#Time", axis=1)
        a = met.to_xarray()
        a = self._rename_variables_revelle(a)
        a = self._add_attributes_revelle(a)
        return a

    def _rename_variables_revelle(self, a):
        names = dict(
            LO="lon",
            LA="lat",
            GY="hdg",
            SP="sog",
            CR="cog",
            AT="air_temperature",
            WD="relative_wind_direction",
            WS="relative_wind_speed",
            RH="relative_humidity",
            TW="true_wind_speed",
            TI="true_wind_direction",
            BP="barometric_pressure",
            SW="shortwave_radiation",
            LW="longwave_radiation",
            SA="sea_surface_salinity",
            ST="sea_surface_temperature",
            FL="fluorometer",
            LF="water_depth_12khz",
            HF="water_depth_35khz",
            MB="water_depth_multibeam",
        )
        for varname, values in a.items():
            if varname not in names.keys():
                a = a.drop(varname)
        a = a.rename_vars(names)
        return a

    def _add_attributes_revelle(self, a):
        # a.spd.attrs = dict(long_name="ship speed", units="knts")
        a.hdg.attrs = dict(long_name="ship heading", units="deg")
        a.cog.attrs = dict(long_name="ship course over ground", units="deg")
        a.sog.attrs = dict(long_name="ship speed over ground", unit="knots")
        a.air_temperature.attrs = dict(long_name="air temperature", units="°C")
        # a.rain_intensity.attrs = dict(
        #     long_name="rain intensity stbd", units="mm/h"
        # )
        # a.rain_accumulation.attrs = dict(
        #     long_name="rain accumulation stbd", units="mm/h"
        # )
        a.relative_wind_direction.attrs = dict(
            long_name="relative wind direction", units="deg"
        )
        a.relative_wind_speed.attrs = dict(long_name="relative wind speed", units="m/s")
        a.relative_humidity.attrs = dict(long_name="relative humidity", units="%")
        a.true_wind_speed.attrs = dict(long_name="true wind speed", units="m/s")
        a.true_wind_direction.attrs = dict(long_name="true wind direction", units="deg")
        a.barometric_pressure.attrs = dict(long_name="barometric pressure", units="mb")
        a.shortwave_radiation.attrs = dict(
            long_name="shortwave radiation", units="W/m$^2$"
        )
        a.longwave_radiation.attrs = dict(
            long_name="longwave radiation", units="W/m$^2$"
        )
        # a.photosynthetically_active_radiation.attrs = dict(
        #     long_name="photosynthetically active radiation", units="uE/m$^2$/s"
        # )
        a.sea_surface_salinity.attrs = dict(
            long_name="sea surface salinity (5m)", units="psu"
        )
        a.sea_surface_temperature.attrs = dict(
            long_name="sea surface temperature (5m)", units="°C"
        )
        a.fluorometer.attrs = dict(long_name="fluorometer", units="ug/l")
        # a.sea_surface_sound_velocity.attrs = dict(
        #     long_name="sea surface sound velocity", units="m/s"
        # )
        a.water_depth_12khz.attrs = dict(long_name="12kHz water depth", units="m")
        a.water_depth_35khz.attrs = dict(long_name="35kHz water depth", units="m")
        a.water_depth_multibeam.attrs = dict(
            long_name="multibeam water depth", units="m"
        )
        return a

    def _utc_yymmdd(self):
        utcnow = gv.time.now_utc()
        utcstr = gv.time.datetime64_to_str(utcnow)
        y, m, d = utcstr.split("-")
        yymmdd = y[2:] + m + d
        return yymmdd

    def _yymmdd_to_datetime64(self, ymd):
        date = f"20{ymd[:2]}-{ymd[2:4]}-{ymd[4:]}"
        return np.datetime64(date)

    def read_all_met(self, start_time=None):
        yymmdd_today = self._utc_yymmdd()
        met_all = self.local_met.glob("*.MET")
        if start_time is not None:
            met_all = [
                file
                for file in met_all
                if self._yymmdd_to_datetime64(file.stem) >= start_time
            ]
        readall = []
        for f in sorted(met_all):
            ncfile = f.with_suffix(".nc")
            if f.stem == yymmdd_today:
                met = self._read_met_file(f)
                # met.to_netcdf(ncfile)
                # met.close()
                readall.append(met)
            elif ncfile.exists():
                # read existing netcdf
                tmp = xr.open_dataset(ncfile)
                if tmp.time.shape[0] > 5758:
                    readall.append(tmp)
                else:
                    print(f"update {f.name}")
                    met = self._read_met_file(f)
                    ncfile.unlink()
                    met.to_netcdf(ncfile)
                    met.close()
                    readall.append(met)
            else:
                print(f.name)
                # read raw and write to netcdf if it does not exist
                met = self._read_met_file(f)
                met.to_netcdf(ncfile)
                met.close()
                readall.append(met)
        b = xr.concat(readall, dim="time")
        b = self._add_attributes_revelle(b)
        b = self._met_remove_outliers(b)
        return b

    def _met_remove_outliers(self, met):
        for v in ["true_wind_speed", "longwave_radiation"]:
            met[v] = met[v].where(met[v] >= 0)
        return met

    def current_location(self):
        return self._read_gps_packet()

    def _parse_gps_packet(self, gps):
        sa = gps.split("\n")
        for si in sa:
            if si.startswith("$GPZDA"):
                s = si.split(",")
                time = s[1]
                dts = f"{s[4]}-{s[3]}-{s[2]} {time[:2]}:{time[2:4]}:{time[4:]}"
                dt64 = np.datetime64(dts)
            if si.startswith("$GPGGA"):
                s = si.split(",")
                lats, latl, lons, lonl = s[2:6]
                lon = np.float64(lons[:3]) + np.float64(lons[3:]) / 60
                lat = np.float64(lats[:2]) + np.float64(lats[2:]) / 60
                lon = -lon if lonl == "W" else lon
        # return dt64, lon, lat
        return lon, lat

    def _read_gps_packet(self):
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        port = 50001
        sock.bind(("", port))
        keep_running = True
        while keep_running:
            s = sock.recv(1024).decode("utf-8")
            if s.endswith("\n"):
                out = ""
                out = [out + sock.recv(1024).decode("utf-8") for i in range(4)]
                out = "".join(out)
                keep_running = False
        sock.close()
        return self._parse_gps_packet(out)


class Armstrong(Underway):
    # Save the two-digit year as a constant to use later on
    _YY = datetime.datetime.utcnow().strftime("%y")

    def __init__(self, local_data, atsea=True, cruise_id=None, bathy=None):
        """Generate Armstrong underway object.

        Parameters
        ----------
        local_data : PosixPath
            Path to local data directory. Will have subdirectories "met", "ctd"
            and "sadcp".
        atsea : bool, optional
            Set to True to sync data. Defaults to True.
        cruise_id : str, optional
            Cruise ID.
        """
        super().__init__(local_data, atsea, cruise_id)
        self.ship = "armstrong"
        self.parse_raw_gps()

    def set_quickplot_variables(self):
        self.metvar_wind_speed = "true_wind_speed_stbd"

    def connect(self):
        underway.network.connect_servers_revelle()

    def position(self):
        return underway.network.get_position_armstrong()

    def sync_data(self):
        self.remote_met = Path("/Volumes/data_on_memory/underway/proc")
        self.remote_ctd = Path("/Volumes/data_on_memory/ctd/")
        self.remote_ladcp = Path("/Volumes/science_share/LADCP/")
        # set the pattern of the files we are looking for
        self.gps_pattern = f"AR_GPS10_{self._YY}*.csv"

        self.connect()
        print("syncing sadcp data...")
        self.sync_sadcp()
        # print("syncing met data...")
        # rsync_underway_data(self.remote_met, self.local_met, self.met_pattern)
        self.sync_met()
        print("syncing gps data...")
        rsync_underway_data(self.remote_met, self.local_gps, self.gps_pattern)
        print("syncing gps raw data...")
        self.sync_raw_gps()
        print("syncing ctd data...")
        self.sync_ctd_data()
        if self.cruise_id == "ar73":
            print("syncing ladcp data...")
            self.ar73_sync_ladcp()

    def sync_ctd_data(self, verbose=False):
        """Sync CTD data from ship server."""
        _sync_ctd_data(self.remote_ctd, self.local_ctd_raw, verbose=verbose)

    def sync_met(self):
        # set the pattern of the files we are looking for
        self.met_pattern = f"AR{self._YY}*.csv"
        print("syncing met data...")
        rsync_underway_data(self.remote_met, self.local_met, self.met_pattern)

    def sync_sadcp(self, verbose=False):
        """Sync R/V Armstrong shipboard ADCP data

        Parameters
        ----------
        verbose : bool, optional
            Show rsync output.
        """
        wh300path = Path("/Volumes/data_on_memory/adcp/proc/wh300/contour/wh300.nc")
        os150nbpath = Path(
            "/Volumes/data_on_memory/adcp/proc/os150nb/contour/os150nb.nc"
        )
        os38nbpath = Path("/Volumes/data_on_memory/adcp/proc/os38nb/contour/os38nb.nc")
        rsync_opts = "-av" if verbose else "-a"
        subprocess.call(["rsync", rsync_opts, wh300path.as_posix(), self.local_sadcp])
        subprocess.call(["rsync", rsync_opts, os150nbpath.as_posix(), self.local_sadcp])
        subprocess.call(["rsync", rsync_opts, os38nbpath.as_posix(), self.local_sadcp])

    def sync_raw_gps(self):
        remote_gps_raw = Path("/Volumes/data_on_memory/underway/raw/")
        sources = dict(
            gps=dict(dir=remote_gps_raw, pattern="*.CNAV_3050"),
        )
        for src, p in sources.items():
            remote = p["dir"]
            copy_underway_data(remote, self.local_gps_raw, p["pattern"])

    def sync_raw_gps_rsync(self):
        local_raw_gps = self.local_gps_raw
        subprocess.call(
            [
                "rsync",
                "-a",
                "--include=*.CNAV_3050",
                "--exclude=*",
                "/Volumes/data_on_memory/underway/raw/",
                local_raw_gps.as_posix() + "/",
            ]
        )

    def _parse_raw_gps_file(self, file):
        with open(file) as file:
            lines = file.readlines()
        time = []
        lon = []
        lat = []
        for line in lines:
            dt, gp = line.split(" CNAV ")
            if gp[:6] == "$GPRMC":
                msg = NMEAReader.parse(gp)
                time.append(pd.to_datetime(dt[4:]).to_numpy())
                lon.append(msg.lon)
                lat.append(msg.lat)
        time = np.array(time)
        lon = np.array(lon)
        lat = np.array(lat)
        gps = xr.Dataset(
            data_vars=dict(lon=(("time"), lon), lat=(("time"), lat)),
            coords=dict(time=(("time"), time)),
        )
        return gps

    def parse_raw_gps(self):
        print("parsing raw gps files")
        raw_gps = self.local_gps.joinpath("raw")
        files = sorted(raw_gps.glob("*.CNAV*"))
        ncfiles = sorted(raw_gps.glob("*.nc"))
        files = sorted(files)
        for file in files:
            if file.with_suffix(".nc") not in ncfiles:
                print(file.stem)
                gps = self._parse_raw_gps_file(file)
                gps.to_netcdf(file.with_suffix(".nc"))
        # how do we deal with the file that is not completed?
        # just redo the last nc file.
        file = ncfiles[-1].with_suffix(".CNAV_3050")
        gps = self._parse_raw_gps_file(file)
        gps.to_netcdf(file.with_suffix(".nc"))

        # read all nc files
        allnav = [xr.open_dataset(file) for file in ncfiles]
        [navi.close() for navi in allnav]
        nav = xr.concat(allnav, dim="time")
        # nav = xr.open_mfdataset(ncfiles).load()
        # nav.close()
        self.nav = nav

    def save_nav(self):
        self.nav.to_netcdf(self.local_gps_proc.joinpath("ar73_gps.nc"))

    def _read_met_file_armstrong(self, file):
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
        a = self._rename_variables_armstrong(a)
        return a

    def _rename_variables_armstrong(self, a):
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

    def read_all_met(self):
        met_all = self.local_met.glob("*.csv")
        readall = []
        for f in sorted(met_all):
            readall.append(self._read_met_file_armstrong(f))
        b = xr.concat(readall, dim="time")
        b = self._add_attributes_armstrong(b)
        return b

    def _add_attributes_armstrong(self, a):
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

    def ar73_sync_ladcp(self, verbose=False):
        _sync_ctd_data(self.remote_ladcp, self.local_ladcp_raw, verbose=verbose)


class Discovery(Underway):
    def __init__(self, local_data, atsea=True, cruise_id=None, bathy=None):
        """Generate Armstrong underway object.

        Parameters
        ----------
        local_data : PosixPath
            Path to local data directory. Will have subdirectories "met", "ctd"
            and "sadcp".
        atsea : bool, optional
            Set to True to sync data. Defaults to True.
        cruise_id : str, optional
            Cruise ID.
        """
        super().__init__(local_data, atsea, cruise_id)
        self.ship = "discovery"

    def set_quickplot_variables(self):
        pass

    def connect(self):
        underway.network.connect_servers_discovery()

    def sync_data(self):
        self.connect()
        self.sync_met_disco()
        sync_sadcp_disco(local_sadcp)
        sync_ctd_disco(local_ctd)
        sync_ladcp_disco(local_ladcp)

    def sync_met_disco(self):
        techsas = Path("/Volumes/current_cruise/Ship_Systems/Data/TechSAS/NetCDF/")
        sources = dict(
            met=dict(dir="SURFMETV3", pattern="*.SURFMETv3"),
            gps=dict(dir="GPS", pattern="*position-POSMV_GPS.gps"),
            tsg=dict(dir="TSG", pattern="*SBE45-SBE45.TSG"),
        )
        for src, p in sources.items():
            remote = techsas / p["dir"]
            copy_underway_data(remote, self.local_met, p["pattern"])

    def sync_gps_disco(self):
        techsas = Path("/Volumes/current_cruise/Ship_Systems/Data/TechSAS/NetCDF/")
        sources = dict(
            gps=dict(dir="GPS", pattern="*position-POSMV_GPS.gps"),
        )
        for src, p in sources.items():
            remote = techsas / p["dir"]
            copy_underway_data(remote, self.local_met, p["pattern"])

    def sync_sadcp_disco(self):
        """currently only running bb"""
        adcp = Path(
            # "/Volumes/current_cruise/Ship_Systems/Data/Acoustics/ADCP/UHDAS/DY132_part2/proc/"
            # "/Volumes/current_cruise/Ship_Systems/Data/Acoustics/ADCP/UHDAS/DY138/proc/"
            "/Volumes/current_cruise/Ship_Systems/Data/Acoustics/ADCP/UHDAS/DY153/proc/"
        )
        sources = ["os75bb", "os150bb", "os75nb", "os150nb"]
        sources = ["os75nb", "os150nb"]
        for src in sources:
            files = [src + ".nc", "contour_uv.mat", "contour_xy.mat"]
            for file in files:
                remote = adcp / src / "contour" / file
                try:
                    copy_underway_data(remote, self.local_sadcp / src)
                except PermissionError:
                    pass

    def sync_ctd_disco(self):
        remote_ctd = Path(
            # "/Volumes/current_cruise/Sensors_and_Moorings/DY132/CTD/Data/CTD Raw Data"
            "/Volumes/current_cruise/Sensors_and_Moorings/DY138/CTD/Data/CTD Raw Data"
        )
        # sync_ctd_data(remote_ctd, local_ctd)
        rsync_underway_data(remote_ctd, self.local_ctd_raw, "*")

    def sync_ladcp_disco(self):
        remote_ladcp = Path(
            # "/Volumes/current_cruise/Sensors_and_Moorings/DY132/LADCP"
            "/Volumes/current_cruise/Sensors_and_Moorings/DY138/LADCP"
        )
        remote_primary = remote_ladcp / "Master data"
        remote_secondary = remote_ladcp / "Slave data"
        # sync ladcp data
        rsync_underway_data(remote_primary, self.local_ladcp / "raw" / "primary", "*")
        rsync_underway_data(
            remote_secondary, self.local_ladcp / "raw" / "secondary", "*"
        )

    def read_all_met(self):
        sources = dict(
            light="*Light-SURFMET.SURFMETv3",
            met="*MET-SURFMET.SURFMETv3",
            surf="*Surf-SURFMET.SURFMETv3",
            gps="*position-POSMV_GPS.gps",
            tsg="*SBE45-SBE45.TSG",
        )
        out = dict()
        for key, pattern in sources.items():
            files = sorted(self.local_met.glob(pattern))
            out[key] = combine_netcdf(files)

        # fix double time stamp bug
        _, ni = np.unique(out["met"].time, return_index=True)
        ni = ~np.isnat(out["met"].time).data
        out["met"] = out["met"].isel(time=ni)
        _, ni = np.unique(out["light"].time, return_index=True)
        out["light"] = out["light"].isel(time=ni)
        _, ni = np.unique(out["surf"].time, return_index=True)
        out["surf"] = out["surf"].isel(time=ni)
        _, ni = np.unique(out["tsg"].time, return_index=True)
        out["tsg"] = out["tsg"].isel(time=ni)
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


def rsync_underway_data(remotedir, local_met, pattern, verbose=False):
    """Sync underway data from ship server using rsync."""
    # print("syncing data from ship server")
    rsync_opts = "-av" if verbose else "-a"
    for f in sorted(list(remotedir.glob(pattern))):
        subprocess.call(["rsync", rsync_opts, f, local_met])


def copy_underway_data(remotedir, local_dir, pattern=None, verbose=False):
    """Copy underway data from ship server."""
    # print("syncing data from ship server")
    if pattern:
        files = sorted(list(remotedir.glob(pattern)))
    else:
        files = [remotedir]
    for f in files:
        # subprocess.call(["rsync", "-avz", f, local_dir])
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
    subprocess.call(["chflags", "-R", "nouchg", f"{local_dir}"])


def _sync_ctd_data(remote_ctd, local_ctd, verbose=False):
    """Sync CTD data from ship server."""
    # print("syncing ctd data from ship server")
    rsync_opts = "-av" if verbose else "-a"
    subprocess.call(["rsync", rsync_opts, remote_ctd.as_posix() + "/", local_ctd])
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


def parse_datetime(dt):
    return np.datetime64(dt.replace(tzinfo=None))

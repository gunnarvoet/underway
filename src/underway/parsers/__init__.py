"""Parsers for raw shipboard files. Each maps a list of files to a Dataset."""

from . import armstrong, lds, nmea, seapath, techsas

__all__ = ["armstrong", "lds", "nmea", "seapath", "techsas"]

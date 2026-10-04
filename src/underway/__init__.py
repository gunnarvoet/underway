import importlib.metadata

from . import parsers
from .cruise import Cruise
from .source import Source

__all__ = ["Cruise", "Source", "parsers"]

__author__ = "Gunnar Voet"
__email__ = "gvoet@ucsd.edu"
# version is defined in pyproject.toml
__version__ = importlib.metadata.version("underway")

import importlib.metadata
from . import io, network, ship, utils

__all__ = ["io", "network", "ship", "utils"]

__author__ = "Gunnar Voet"
__email__ = "gvoet@ucsd.edu"
# version is defined in pyproject.toml
__version__ = importlib.metadata.version("underway")


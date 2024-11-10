__all__ = ["io", "network", "ship", "utils"]

__author__ = "Gunnar Voet"
__email__ = "gvoet@ucsd.edu"
__version__ = "2024.11"

# workaround for when whatever is defined as the default backend is not around:
try:
    import matplotlib.pyplot as plt
except ImportError:
    import matplotlib as mpl

    mpl.use("Agg")

from . import io, network, ship, utils

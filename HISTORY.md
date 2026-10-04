# History

## 2026.10

-   Rewrote the package around a `Cruise` object with declared sources. Ships are tables of sources in `ships.py`. Parsers are plain functions per file format in `parsers`.
-   Common variable names and units across ships.
-   netCDF cache keyed on raw file size replaces the line-count resume and the size thresholds.
-   Parsers for all four ships tested against raw files.
-   Removed `io`, `ship`, `network`, `utils`, the plot helpers, and the live position functions. The old interface is at tag `v2025.10`.
-   Requires Python 3.11.

## 2025.10

-   Updates during cruise DY202 on RRS Discovery.
-   Switch to `uv` for development and backend.

## 2024.11

-   Changed to date versioning
-   Rewrote the code to have separate classes for each research vessel. They are based on an abstract base class that defines general features and specifications for the individual classes, however, I can also do more vessel-specific stuff this way. The new code is in submodule `ship` whereas the old code still lives in `io` in case I need it again.

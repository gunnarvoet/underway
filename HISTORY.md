# History

## 2025.10

-   Updates during cruise DY202 on RRS Discovery.
-   Switch to `uv` for development and backend.

## 2024.11

-   Changed to date versioning
-   Rewrote the code to have separate classes for each research vessel. They are based on an abstract base class that defines general features and specifications for the individual classes, however, I can also do more vessel-specific stuff this way. The new code is in submodule `ship` whereas the old code still lives in `io` in case I need it again.

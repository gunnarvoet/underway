# underway

Research vessel underway data handling. Everything is under heavy development and may change at any time.

# Features
* Connect, reconnect to ship server (currently macos only, uses Apple's `osascript`) 

* Sync shipboard data (met, sadcp, ctd) using `rsync`

* Process met data into xarray.Dataset

* Process sadcp data into xarray.Dataset

# Ships
Currently supporting

* R/V Neil Armstrong
* RRS Discovery
* R/V Roger Revelle
* R/V Sikuliaq

# Example
The following example sets up an object of the `Sikuliaq` class.
It will create subfolders in the `local_data` directory.

```python
import underway as uw
S = uw.ship.Sikuliaq(
    local_data="/path/to/your/local/cruise/dir/",
    cruise_id="SKQ202417S",
)
# connect to ship server
S.connect()
# sync data to local computer
S.sync_data()
# process gps data and read as xarray.Dataset
gps = S.read_gps_data()
```

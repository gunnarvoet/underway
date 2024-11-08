# underway

Research vessel underway data handling. Everything is under heavy development and may change at any time.

* Free software: MIT license

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

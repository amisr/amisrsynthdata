# amisrsynthdata
Module for generating synethetic AMISR data files

## Installation
1. Clone repo
2. In the top directory, install with `pip install .`

## Basic Usage
This package installs the command line tool `amisrsynthdata`, which is used along with a config file to generate an output hdf5 AMISR data file.  The config file specified the ionosphere state and radar configuration that should be used.
```
$ amisrsynthdata config.ini
```

# amisrsynthdata
Module for generating synethetic AMISR data files

## Installation
1. Clone repo
```
$ git clone https://github.com/amisr/amisrsynthdata.git
```
2. Enter repo and switch to the `develop` branch
```
$ cd amisrsynthdata
$ git checkout develop
```
3. Install with pip
```
pip install .
```

## Basic Usage
This package installs the command line tool `amisrsynthdata`, which is used along with a config file to generate an output hdf5 AMISR data file.  The config file specified the ionosphere state and radar configuration that should be used.
```
$ amisrsynthdata config.ini
```
Example configuration files for both PFISR (`pfisr_synth_config.ini`) and RISR-N (`risrn_synth_config.ini`) are provided.

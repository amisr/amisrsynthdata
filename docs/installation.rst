.. installation.rst

Installation
============

The amisrsynthdata package is pure python and can be installed easily with pip::

  pip install git+https://github.com/amisr/amisrsynthdata.git

Alternatively, if you would like to develop on amisrsythdata, first clone the repo, then install::

  git clone https://github.com/amisr/amisrsynthdata.git
  cd amisrsythdata
  pip install -e .


Requirements
------------
The following additional packages are required and will be installed:

- numpy
- h5py
- pymap3d
- apexpy

Additionally, cartopy is required to create the optional summary plots

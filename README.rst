amisrsynthdata
==============

Module for generating synthetic AMISR data files

This module provides tools to let you easily create synthetic data files for the AMISR (Advanced Module Incoherent Scatter Radar) systems.  The files are based on both a specified ionospheric state and a radar configuration.  This can be used to generate synthetic data in the "SRI data format" both for the three existing AMISRs and for hypothetical future "AMISR-like" systems.  Primarily, it was designed to help test the functionality of various inversion algorithms that attempt to create a complete picture of ionospheric state parameters from discrete measurements by creating a way to check the output of these algorithms against known "truth" data.  Please note that this module does NOT attempt to simulate any fundamental ISR theory.

Quick Start
-----------

Installation
************

The amisrsynthdata package is pure python and can be installed easily with pip::

  $ pip install git+https://github.com/amisr/amisrsynthdata.git


Basic Usage
***********

This package installs the command line tool `amisrsynthdata`, which is used along with a config file to generate an output hdf5 AMISR data file.  The config file specified the ionosphere state and radar configuration that should be used::

  $ amisrsynthdata config.ini

Example configuration files for both PFISR and RISR-N are provided.

Limitations
-----------

The following are NOT currently included in the amisrsynthdata module:

1. Any kind of proper treatment or simulation of ISR theory - The module effectively assumes the radar measures plasma parameters perfectly, although empirical errors can be added.
2. Integration over a time period or smearing along the length of pulses, as well as pulse coding.
3. Madrigal data format - Currently files are only generated in the SRI data format.

Contributing
------------

Contributions to this package are welcome and encouraged, particularly to expand the currently set of specified ionospheres.  Create a pull request to submit contributions, or open an issue if you would like to request new features or report a bug.

To add a new ionospheric state function, write the function as a method of the Ionosphere class.  These functions should take three inputs: geodetic latitude (degrees), geodetic longitude (degrees), and geodetic altitude (meters).  All other inputs should be specified in the config file and accessed through the ``params`` dictionaries.

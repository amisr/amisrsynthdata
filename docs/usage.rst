.. usage.rst

Usage
=====

Command Line
------------

To run the AMISR synthetic data generator, simply call `amisrsynthdata` from the command line with a configuration file::

  amisrsynthdata config.yaml


See :ref:`Configuration File` section of the documentation for details on the field that should be in this file.

Python
------

The modules of amisrsynthdata can be imported directly into python programs.  This is most useful to get truth values at arbitrary points and compare with inversion and interpolation codes.  Create YAML config object to initialize the amisrsynthdta modules by loading an existing `config.yaml` file::

  import yaml

  config_file = 'config.yaml'
  with open(config_file, 'r') as cf:
      config = yaml.load(cf, Loader=yaml.FullLoader)

Alternatively, hard-code a dictionary with the appropriate fields.  It is important to note that the module classes are initialized with an actual YAML object or python dictionary, NOT just the filename string as is used with the command line option.

The amisrsynthdata package contains three classes: `Radar`, `Ionosphere`, and `SyntheticData`.  All three are initialized with the YAML config object discussed above, but all fields may not be necessary for the `Radar` and `Ionosphere` classes.  In general, `Ionosphere` contains functions for calculating ionospheric parameters at specified locations, `Radar` contains information about the radar's position and the locations of all beams and range gates, and `SyntheticData` actually calculates ionospheric parameters at each radar gate location and creates output synthetic data files.

`Radar` example::

  from amisrsynthdata import Radar

  rad = Radar(config)

  # Print BeamCode array
  print(rad.beam_codes)

  # Print geodetic latitude of each range gate
  print(rad.lat)

`Ionosphere` example::

  import datetime as dt
  from amisrsynthdata import Ionosphere

  iono = Ionosphere(config)

  # Calculate ionospheric parameters at a particular location and time
  glat = 65.0
  glon = 100.0
  alt = 300000.
  utime = (dt.datetime(2016, 9, 13, 0, 5, 0)-dt.datetime.utcfromtimestamp(0)).total_seconds()
  Ne = iono.density(utime, glat, glon, alt)
  Vi = iono.velocity(utime, glat, glon, alt)
  Te = iono.etemp(utime, glat, glon, alt)
  Ti = iono.itemp(utime, glat, glon, alt)

`SyntheticData` example::

  from amisrsynthdata import SyntheticData

  sd = SyntheticData(config)

  # Get get output range gate positions and electron density
  fov_glat = sd.Geomag['Latitude']
  fov_glon = sd.Geomag['Longitude']
  fov_galt = sd.Geomag['Altitude']
  fov_ne = sd.FittedParams['Ne']

  # Can also get range gate positions from Radar
  glat = sd.radar.lat
  glon = sd.radar.lon
  galt = sd.radar.alt

  # And access Ionosphere functions directly
  Ne = sd.iono.density(utime, glat, glon, galt)

Note that all functionality of the `Radar` and `Ionosphere` classes are available through the `SyntheticData` class.  `SyntheticData` contains an instance of the `Radar` class named `radar` and an instance of the `Ionosphere` class named `iono`.

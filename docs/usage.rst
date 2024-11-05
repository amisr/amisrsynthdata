.. usage.rst

Usage
=====

Command Line
------------

To run the AMISR synthetic data generator, simply call ``amisrsynthdata`` from the command line with a configuration file.

.. code-block::

  amisrsynthdata config.yaml


See :ref:`Configuration File` section of the documentation for details on the fields that should be in this file.

Python
------

The modules of amisrsynthdata can be imported directly into python programs.  This is most useful to get truth values at arbitrary points and compare with inversion and interpolation codes.  Create YAML config object to initialize the amisrsynthdata modules by loading an existing ``config.yaml`` file.

.. code-block:: python

  import yaml

  config_file = 'config.yaml'
  with open(config_file, 'r') as cf:
      config = yaml.load(cf, Loader=yaml.FullLoader)

Alternatively, hard-code a dictionary with the appropriate fields.  It is important to note that the module classes are initialized with an actual YAML object or python dictionary, NOT just the filename string as is used with the command line option.

The amisrsynthdata package contains three classes: ``Radar``, ``Ionosphere``, and ``SyntheticData``.  All three are initialized with the YAML config object discussed above, but all fields may not be necessary for the ``Radar`` and ``Ionosphere`` classes.  In general, ``Ionosphere`` contains functions for calculating ionospheric parameters at specified locations, ``Radar`` contains information about the radar's position and the locations of all beams and range gates, and ``SyntheticData`` actually calculates ionospheric parameters at each radar gate location and creates output synthetic data files.

`Radar` Class Example
*********************

.. code-block:: python

  from amisrsynthdata.radar import Radar

  rad = Radar(config)

  # Print BeamCode array
  print(rad.beam_codes)

  # Print geodetic latitude of each range gate
  print(rad.lat)

`Ionosphere` Class Example
**************************

.. code-block:: python

  import datetime as dt
  from amisrsynthdata.ionosphere import Ionosphere

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

`SyntheticData` Class Example
*****************************

.. code-block:: python

  from amisrsynthdata.syntheticdata import SyntheticData

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

Note that all functionality of the ``Radar`` and ``Ionosphere`` classes are available through the ``SyntheticData`` class.  ``SyntheticData`` contains an instance of the ``Radar`` class named ``radar`` and an instance of the ``Ionosphere`` class named ``iono``.

Benchmarking
------------

The amount of time it takes to run ``amisrsynthdata`` depends on the ionosphere model that is chosen as well as the radar mode. Modes with more beams or finner range or time resolution will generally take longer to compute synthetic data files.  Producing summary plots also increases the time it takes to run ``amisrsynthdata``.  The table below shows rough benchmarking of how long it takes to produce a synthetic data file from the command line for a simple case and a complex case.  The simple case uses the `example configuration file <https://github.com/amisr/amisrsynthdata/blob/develop/example_synth_config.yaml>`_ provided with the package which uses altitude-varying or uniform ionospheric state functions and 6 beams with relatively corse range resolution.  The complex case uses an ionosphere specified from the output of the GEMINI numerical model (the slowest ionosphere option currently available) and mimics the 52 beam imaging mode.  This benchmarking was performed on a laptop workstation and should only be considered approximate.

+--------------+--------+----------+
|              | Simple | Complex  |
+==============+========+==========+
| **No Plots** | 0.77 s |  87.48 s |
+--------------+--------+----------+
| **Plots**    | 9.72 s | 103.24 s |
+--------------+--------+----------+


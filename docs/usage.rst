.. usage.rst

Usage
=====

To run the AMISR synthetic data generator, simpily call `amisrsynthdata` from the command line with a config file::

  amisrsynthdata config.ini


See Configuration File section of the documentation for details on the field that should be in this file.

Alternatively, the modules of amisrsynthdata can be imported directly into python programs.  This is most useful to use the `Ionosphere` class to get truth values at arbitrary points and compare with the results of various inversion techniques.  All classes should be initialized with the name of a config file, similar to how the command line program is run.  Then the methods of the class can be used to determine the truth values at various locations.

Simple code example ::

  from amisrsynthata.Ionosphere import Ionosphere

  config_filg = '/path/to/config_file.ini'
  iono = Ionosphere(config_file)

  glat = 65.0
  glon = 100.0
  alt = 300000.
  ne = iono.density(glat, glon, alt)

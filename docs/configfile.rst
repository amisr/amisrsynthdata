.. configfile.rst

.. _Configuration File:

Configuration File
==================

The configuration file specifies both the state of the ionosphere, the radar mode and configuration, and various other options for creating the synthetic data file.  The file should be in `YAML <https://yaml.org/>`_ and an example configuration file is available in the main project repository.  The configuration file requires six main sections: GENERAL (which specifies various options for how the synthetic data generator should be run), RADAR (which specifies the radar mode and configuration options), and DENSITY, VELOCITY, ITEMP, and ETEMP (which specify the values for electron density, plasma drift velocity, ion temperature, and electron temperature, respectively).  The GENERAL and RADAR sections are composed of several standard parameters that need to be specified (described in the tables below).

**IMPORTANT NOTE**: YAML can natively handle scientific notation, but has stringent formatting requirements.  The mantissa must have a decimal point and the sign of the exponent must be included.  This means `4.0e+11` can be used in the configuration files, but `4e11` would generate errors.

GENERAL
-------

+---------------------+--------------------------------------------------------+----------------------------------+
| Parameter           | Description                                            | Example                          |
+=====================+========================================================+==================================+
| starttime           | File start time (in ISO format)                        | 2021-09-13 00:00:00              |
+---------------------+--------------------------------------------------------+----------------------------------+
| endtime             | File end time (in ISO format)                          | 2021-09-13 00:05:00              |
+---------------------+--------------------------------------------------------+----------------------------------+
| output_filename     | Name out output synthetic data file                    | amisr_synthetic_data_output.h5   |
+---------------------+--------------------------------------------------------+----------------------------------+
| ion_mass [2]_       | List of masses of ions to be included in composition   | [16.,32.,30.,28.,14.]            |
+---------------------+--------------------------------------------------------+----------------------------------+
| err_coef            | Coefficients for r^2 empirical errors (Ne, Ve, Te, Ti) | [1., 1.e-9, 5.e-9, 1.e-9]        |
+---------------------+--------------------------------------------------------+----------------------------------+
| noise               | Whether or not to add random noise                     | false                            |
+---------------------+--------------------------------------------------------+----------------------------------+
| summary_plot        | Filename for output summary plot                       | synthetic_data_summary_risrn.png |
+---------------------+--------------------------------------------------------+----------------------------------+
| summary_plot_time   | Time to plot in the output summary plot                |  2021-09-13 00:03:00             |
+---------------------+--------------------------------------------------------+----------------------------------+


RADAR
-----

+-------------------------+-----------------------------------------------------------+----------------------------+
| Parameter               | Description                                               | Example                    |
+=========================+===========================================================+============================+
| full_name               | Site Name                                                 | Poker Flat                 |
+-------------------------+-----------------------------------------------------------+----------------------------+
| abbreviation            | Site Name                                                 | PFISR                      |
+-------------------------+-----------------------------------------------------------+----------------------------+
| site_coordinates        | Site geodetic coordinates - lat (deg), lon (deg), alt (m) | [65.13, -147.47, 213.]     |
+-------------------------+-----------------------------------------------------------+----------------------------+
| beamcodes [3]_          | Beamcodes to use                                          | [64016, 64157, 64964]      |
+-------------------------+-----------------------------------------------------------+----------------------------+
| beam_azimuth  [3]_      | Azimuth angles of beams to define (deg)                   | [120., 180.]               |
+-------------------------+-----------------------------------------------------------+----------------------------+
| beam_elevation [3]_     | Elevation angles of beams to define (deg)                 | [60., 50.]                 |
+-------------------------+-----------------------------------------------------------+----------------------------+
| acf_slant_range         | Ranges or unfitted gates - start, stop, step (m)          | [80000., 800000., 3000.]   |
+-------------------------+-----------------------------------------------------------+----------------------------+
| altitude_bins [4]_      | Altitude bins to use for fitted data (m)                  | [100000., 800000., 50000.] |
+-------------------------+-----------------------------------------------------------+----------------------------+
| integration_period [1]_ | Integration period (s)                                    | 60.                        |
+-------------------------+-----------------------------------------------------------+----------------------------+

.. [1] Although no actual integration occurs, integration_period determines the time steps in the output file.

.. [2] Most calculations are done assuming only O+, but listing all five major species considered by the fitter will result in arrays of a similar shape.

.. [3] Beams can be specified either by a list of beamcode numbers, a list of azimuth and elevation angles, or some combination of both.

.. [4] Bins should be specified as start, stop, step, but multiple sets can be specified to change step size for different altitude ranges. See example configuration file.



IONOSPHERE
----------

The four ionospheric state sections (DENSITY, VELOCITY, ETEMP, ITEMP) should each contain at least one section defining the function that should be used to define that variable.  These sections should be named by the function that will be used and have sub-parameters that specify whatever parameters that function needs.  As an example, the following specifies density should be a standard Chapman profile with NmF2 = 4.0e+11 m-3, hmF2 = 300000 m, a scale height of 100000 m, and solar zenith angle as 0.::

  DENSITY:
    chapman:
      N0: 4.0e+11
      H: 100000.
      z0: 300000.
      sza: 0.


Which parameters are specified in these sections will vary based on the function selected.  Details about which functions are currently available for each ionospheric state parameters and what inputs they need are available in the :ref:`Ionospheric State` section of the documentation.  Multiple functions can be specified for each of the four state parameters, in which case the package will sum all functions evaluated at each point.  This allows the user to create more complicated patterns, such as a patch on top of a background Chapman layer.

Refer to the API references and the example configuration file for assistance generating the ionospheric state sections, however there are a few general tips to keep in mind for constructing a sensible ionosphere.
  - If a parameter is not important for your use case, it is simplest to just set it to some reasonable, uniform value.
  - However, remember that some codes filter data based on other parameters (i.e., electron density), so make sure any "filler" values appropriate for the use case.
  - When specifying multiple functions for one parameter, the result is the SUM of each function individually, so a Chapman layer on top of a uniform background will increase the peak density of the Chapman layer.
  - Both ion and electron temperature can be specified from the same set of functions from the `Temperature` class, however, different functions, or the same function with different parameters, can be used for each.
  - If a state function does not exist for a particular ionospheric structure, you can write a new one!

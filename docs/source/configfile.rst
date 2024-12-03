.. configfile.rst

.. _Configuration File:

Configuration File
==================

The configuration file specifies both the state of the ionosphere, the radar mode and configuration, and various other options for creating the synthetic data file.  The file should be in `YAML <https://yaml.org/>`_ and an `example configuration file <https://github.com/amisr/amisrsynthdata/blob/main/example_synth_config.yaml>`_ is available in the main project repository.  The contents of this file can be copy/pasted into a local text file or the example file can be downloaded from GitHub or with wget on the command line.

.. code-block::

   wget https://raw.githubusercontent.com/amisr/amisrsynthdata/refs/heads/main/example_synth_config.yaml

This file should work as an example for running amisrsynthdata, but it will have to be modified to synthesize data for different situations.  Note that you can rename this configuration file and keep several on your system.

The configuration file requires six main sections: GENERAL (which specifies various options for how the synthetic data generator should be run), RADAR (which specifies the radar mode and configuration options), and DENSITY, VELOCITY, ITEMP, and ETEMP (which specify the values for electron density, plasma drift velocity, ion temperature, and electron temperature, respectively).  There is an additional optional SUMMARY_PLOT section where plotting parameters for the automatically-generated summary plots can be specified.  The GENERAL and RADAR sections are composed of several standard parameters that need to be specified (described in the tables below).

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
| ion_mass [1]_       | List of masses of ions to be included in composition   | [16., 32., 30., 28., 14.]        |
+---------------------+--------------------------------------------------------+----------------------------------+
| rel_err             | Relative error which will be scaled as r :sup:`2`      | 0.1                              |
+---------------------+--------------------------------------------------------+----------------------------------+
| err_ref_rng         | Reference range to extend the r :sup:`2` error from (m)| 300000.                          |
+---------------------+--------------------------------------------------------+----------------------------------+
| noise               | Whether or not to add random noise                     | False                            |
+---------------------+--------------------------------------------------------+----------------------------------+


RADAR
-----

+-------------------------+-----------------------------------------------------------+----------------------------+
| Parameter               | Description                                               | Example                    |
+=========================+===========================================================+============================+
| full_name               | Site Name                                                 | Poker Flat                 |
+-------------------------+-----------------------------------------------------------+----------------------------+
| abbreviation [2]_       | Site Abbreviation                                         | PFISR                      |
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
| integration_period [5]_ | Integration period (s)                                    | 60.                        |
+-------------------------+-----------------------------------------------------------+----------------------------+


.. [1] Most calculations are done assuming only O+, but listing all five major species considered by the fitter will result in arrays of a similar shape.

.. [2] This is used to look up beamcode tables in addition to being saved to the output datafile metadata.  If a non-AMISR site is listed (not PFISR, RISR-N, or RISR-C), you will not be able to specify beams by their beamcodes and must use the azimuth and elevation options instead.

.. [3] Beams can be specified either by a list of beamcode numbers, a list of azimuth and elevation angles, or some combination of both.  Beams specified by azimuth and elevation are assigned a beamcode starting with 90001 and incrimenting sequentually.

.. [4] Bins should be specified as start, stop, step, but multiple sets can be specified to change step size for different altitude ranges. See example configuration file.

.. [5] Although no actual integration occurs, integration_period determines the time steps in the output file.


IONOSPHERE
----------

The four ionospheric state sections (DENSITY, VELOCITY, ETEMP, ITEMP) should each contain at least one section defining the function that should be used to define that variable.  If multiple functions are listed, the sum of the functions will be used.  These sections should be named by the function that will be used and have sub-parameters that specify whatever parameters that function needs.  As an example, the following specifies density should be a standard Chapman profile with NmF2 = 4.0e+11 m :sup:`-3`, hmF2 = 300000 m, a scale height of 100000 m, and solar zenith angle of 0 degrees.

.. code-block::

  DENSITY:
    - chapman:
        N0: 4.0e+11
        H: 100000.
        z0: 300000.
        sza: 0.


Note the ``-`` prepending the function name but the function parameters are listed as a dash-free indented list.  This syntax is very important to follow exactly so the configuration file is read in correctly.  Which parameters are specified in these sections will vary based on the function selected.  Details about which functions are currently available for each ionospheric state parameters and what inputs they need are available in the :ref:`Ionospheric State` section of the documentation.  Multiple functions can be specified for each of the four state parameters, in which case the package will sum all functions evaluated at each point.  This allows the user to create more complicated patterns, such as a patch on top of a background Chapman layer.

Refer to the API references and the example configuration file for assistance generating the ionospheric state sections, however there are a few general tips to keep in mind for constructing a sensible ionosphere.

  * If a parameter is not important for your use case, it is simplest to just set it to some reasonable, uniform value.
  * However, remember that some codes filter data based on other parameters (i.e., electron density), so make sure any "filler" values appropriate for the use case.
  * When specifying multiple functions for one parameter, the result is the SUM of each function individually, so a Chapman layer on top of a uniform background will increase the peak density of the Chapman layer.
  * Both ion and electron temperature can be specified from the same set of functions from the ``Temperature`` class, however, different functions, or the same function with different parameters, can be used for each.
  * If a state function does not exist for a particular ionospheric structure, you can :ref:`write a new one <New State Functions>`!


SUMMARY_PLOT
------------

This section is optional.  If it is not included, summary plots will not be created.

+-------------------------+-----------------------------------------------------------+----------------------------+
| Parameter               | Description                                               | Example                    |
+=========================+===========================================================+============================+
| output_prefix           | Base file name of output summary plots                    | synthdata_summary          |
+-------------------------+-----------------------------------------------------------+----------------------------+
| plot_time               | Target time for altitude slices and 3D plot               | 2016-09-13 00:10:00        |
+-------------------------+-----------------------------------------------------------+----------------------------+
| plot_beam               | Beamcode for RTI plot                                     | 64157                      |
+-------------------------+-----------------------------------------------------------+----------------------------+
| alt_slices              | Altitudes to use for altitude slices (m)                  | [200000., 300000., 400000.]|
+-------------------------+-----------------------------------------------------------+----------------------------+
| slice_xrng              | E-W limits and step side of altitude slice (m)            | [-500000., 500000., 10000.]|
+-------------------------+-----------------------------------------------------------+----------------------------+
| slice_yrng              | N-S limits and step size of altitude slice (m)            | [-450000., 550000., 10000.]| 
+-------------------------+-----------------------------------------------------------+----------------------------+
| dens_colors             | Limits and color map to use for density plots             | vmin: 0                    |
|                         |                                                           |                            |
|                         |                                                           | vmax: 5.0e+11              |
|                         |                                                           |                            |
|                         |                                                           | cmap: viridis              |
+-------------------------+-----------------------------------------------------------+----------------------------+
| itemp_colors            | Limits and color map to use for ion temperature plots     | vmin: 0                    |
|                         |                                                           |                            |
|                         |                                                           | vmax: 3000.                |
|                         |                                                           |                            |
|                         |                                                           | cmap: magma                |
+-------------------------+-----------------------------------------------------------+----------------------------+
| etemp_colors            | Limits and color map to use for electron temperature plots| vmin: 0                    |
|                         |                                                           |                            |
|                         |                                                           | vmax: 5000.                |
|                         |                                                           |                            |
|                         |                                                           | cmap: inferno              |
+-------------------------+-----------------------------------------------------------+----------------------------+
| vlos_colors             | Limits and color map to use for velocity plots            | vmin: -500.                |
|                         |                                                           |                            |
|                         |                                                           | vmax: 500.                 |
|                         |                                                           |                            |
|                         |                                                           | cmap: bwr                  |
+-------------------------+-----------------------------------------------------------+----------------------------+




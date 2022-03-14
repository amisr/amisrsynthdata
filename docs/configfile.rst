.. configfile.rst

Configuration File
==================

The configuration file specifies both the state of the ionosphere, the radar mode and configuration, and various other options for creating the synthetic data file.  Example configuration files for PFISR and RISR-N are provided in the main project repository.  The configuration file requires six main sections: GENERAL (which specifies various options for how the synthetic data generator should be run), RADAR (which specifies the radar mode and configuration options), and DENSITY, VELOCITY, ITEMP, and ETEMP (which specify the values for electron density, plasma drift velocity, ion temperature, and electron temperature, respectively).  The GENERAL and RADAR sections are composed of several standard parameters that need to be specified (described in the tables below), but the four ionospheric state sections should include the parameter TYPE as the name of the function that should be used to define that variable, then any additional parameters that function requires.  Which parameters are specified in these sections will vary based on the function selected (ie, a uniform density field only needs to be defined by a single value, but a Chapman layer will require the NmF2, hmF2, and the scale height).  Details about which functions are currently available for each ionospheric state parameters and what inputs they need are available in the Ionospheric State section of the documentation.

GENERAL
-------

+---------------------+--------------------------------------------------------+----------------------------------+
| Parameter           | Description                                            | Example                          |
+=====================+========================================================+==================================+
| STARTTIME           | File start time (in ISO format)                        | 2021-09-13 00:00:00              |
+---------------------+--------------------------------------------------------+----------------------------------+
| ENDTIME             | File end time (in ISO format)                          | 2021-09-13 00:05:00              |
+---------------------+--------------------------------------------------------+----------------------------------+
| OUTPUT_FILENAME     | Name out output synthetic data file                    | amisr_synthetic_data_output.h5   |
+---------------------+--------------------------------------------------------+----------------------------------+
| ION_MASS  [2]_      | List of masses of ions to be included in composition   | 16.,32.,30.,28.,14.              |
+---------------------+--------------------------------------------------------+----------------------------------+
| ERR_COEF            | Coefficients for r^2 empirical errors (Ne, Ve, Te, Ti) | 1.,1.e-9,5.e-9,1.e-9             |
+---------------------+--------------------------------------------------------+----------------------------------+
| NOISE               | Whether or not to add random noise                     | FALSE                            |
+---------------------+--------------------------------------------------------+----------------------------------+
| SUMMARY_PLOT        | Filename for output summary plot                       | synthetic_data_summary_risrn.png |
+---------------------+--------------------------------------------------------+----------------------------------+
| SUMMARY_PLOT_TIME   | Time to plot in the output summary plot                |  2021-09-13 00:03:00             |
+---------------------+--------------------------------------------------------+----------------------------------+


RADAR
-----

+-------------------------+-----------------------------------------------------------+----------------------------+
| Parameter               | Description                                               | Example                    |
+=========================+===========================================================+============================+
| NAME                    | Site Name                                                 | Poker Flat                 |
+-------------------------+-----------------------------------------------------------+----------------------------+
| SITE_COORDS             | Site geodetic coordinates - lat (deg), lon (deg), alt (m) | 65.13,-147.47,213.         |
+-------------------------+-----------------------------------------------------------+----------------------------+
| BEAMCODE_FILENAME [3]_  | Name of beamcode look-up file                             | bcotable_pfisr.txt         |
+-------------------------+-----------------------------------------------------------+----------------------------+
| BEAMCODES [4]_          | Beamcodes to use                                          | 64016, 64157, 64964, 65066 |
+-------------------------+-----------------------------------------------------------+----------------------------+
| BEAM_AZIMUTH  [4]_      | Azimuth angles of beams to define (deg)                   | 120., 180.                 |
+-------------------------+-----------------------------------------------------------+----------------------------+
| BEAM_ELEVATION [4]_     | Elevation angles of beams to define (deg)                 | 60., 50.                   |
+-------------------------+-----------------------------------------------------------+----------------------------+
| RANGE_START             | Start range for unfitted gates (m)                        | 80000.                     |
+-------------------------+-----------------------------------------------------------+----------------------------+
| RANGE_END               | End range for unfitted gates (m)                          | 800000.                    |
+-------------------------+-----------------------------------------------------------+----------------------------+
| RANGE_STEP              | Range step for unfitted gates (m)                         | 3000.                      |
+-------------------------+-----------------------------------------------------------+----------------------------+
| ALTBINS                 | Altitude bins to use for fitted data                      | 200000.,400000.,600000.    |
+-------------------------+-----------------------------------------------------------+----------------------------+
| INTEGRATION_PERIOD [1]_ | Integration period (s)                                    | 60.                        |
+-------------------------+-----------------------------------------------------------+----------------------------+

.. [1] Although no actual integration occurs, INTEGRATION_PERIOD determines the time steps in the output file.

.. [2] Most calculations are done assuming only O+, but listing all five major species considered by the fitter will result in arrays of a similar shape.

.. [3] Appropriate beamcode tables can be downloaded from:

PFISR: https://amisr.com/amisr/about/about_pfisr/pfisr-specs/

RISR-N: https://amisr.com/amisr/about/resolute-bay-isrs/risrn-specs/

RISR-C: https://www.ucalgary.ca/aurora/projects/risrc

.. [4] Beams can be specified either by a list of beamcode numbers, a list of azimuth and elevation angles, or some combination of both.

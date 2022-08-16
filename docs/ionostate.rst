
.. _Ionospheric State:

Ionospheric State
=================

Density
-------

.. autoclass:: amisrsynthdata.Density.Density
    :members:
    :undoc-members:
    :show-inheritance:


Temperature
-----------
These general temperature functions are used to specify both ion and electron temperature.

.. autoclass:: amisrsynthdata.Temperature.Temperature
    :members:
    :undoc-members:
    :show-inheritance:


Velocity
--------

.. autoclass:: amisrsynthdata.Velocity.Velocity
    :members:
    :undoc-members:
    :show-inheritance:


Create New State Functions
--------------------------

If the existing functions are not sufficient to create the ionospheric state of interest, you can write your own and add them to the package.  Instructions are provided below for adding these functions to your own local version of the code, but the package developers encourage users to create a pull request to contribute their additions to the main version of the package and improve the functionality.

.. _Function Guidelines:

Function Guidelines
*******************

Functions must have a unique name from functions already used by that class, and MUST take in the following variable in this order:
  - Unix Time (seconds since January 1, 1970)
  - Geodetic Latitude (degrees)
  - Geodetic Longitude (degrees)
  - Geodetic Altitude (meters)

Latitude, longitude, and altitude can be multidimensional, but they will always all be the same shape.  Unix Time will be a 1D array of any length.  Functions should return scalar parameters (density and temperature) as arrays where the first dimension matches the length of the Unix Time array, and subsequent dimensions match the shape of the position arrays.  Velocity functions should return a similar array, but with one final dimension that holds the three velocity vector components.

Density values must be returned in units of m-3, temperature values in units of K, and velocity values in units of m/s.  The geodetic components of velocity vectors must be returned in the order east, north, up.  Any transformation to/from other units or coordinate systems must be handled completely within the function.  Parameters input through the configuration file can technically be in any units, but you are encouraged to use standard SI units (meters, seconds, and Kelvin) to be consistent with other functions.

Functions contributed back to the main package MUST include a `numpy-style docstring <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_ that provides a short description of the function and what it should be used for and lists the parameters that must be specified in the configuration file.  ANY parameter your function needs as input aside from the four position/time parameters listed above MUST be in the configuration file.

Instructions for Adding Function
********************************

  1. Clone the source repository and install following the :ref:`developer instructions <developer installation>`.
  2. Create a new branch to develop your new function.
  3. Add your function to the appropriate class  in `Ionosphere.py` following the :ref:`Function Guidelines`.
  4. Modify your configuration file to use your new function.
  5. TEST!!!
  6. Push your new branch to GitHub and create a pull request to merge it into the develop branch of the main package repository.
  7. If you are not interested in contributing your function to the main repository, simply merge the branch you developed the function on with your local main branch.

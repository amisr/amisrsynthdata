
.. _Ionospheric State:

Ionospheric State
=================

All ionospheric state functions take the same input parameters:

  * utime: Unix Time (seconds from January 1, 1970)
  * glat: Geodetic Latitude (degrees)
  * glon: Geodetic Longitude (degrees)
  * galt: Geodetic Altitude (meters)

Other parameters that are needed to define the state should be specified in the appropriate section of the :ref:`configuration file <Configuration File>`.  The parameters that are needed to define different types of ionospheric state functions are listed below.

Density
-------

.. autoclass:: amisrsynthdata.ionosphere.Density
    :members:
    :undoc-members:
    :show-inheritance:


Temperature
-----------
These general temperature functions are used to specify both ion and electron temperature.

.. autoclass:: amisrsynthdata.ionosphere.Temperature
    :members:
    :undoc-members:
    :show-inheritance:


Velocity
--------

.. autoclass:: amisrsynthdata.ionosphere.Velocity
    :members:
    :undoc-members:
    :show-inheritance:



.. _New State Functions:

Create New State Functions
--------------------------

If the existing functions are not sufficient to create the ionospheric state of interest, you can write your own and add them to the package.  Instructions are provided below for adding these functions to your own local version of the code, but the package developers encourage users to create a pull request to contribute their additions to the main version of the package and improve the functionality.

Overview
********

In the source code the directory ``amisrsynthdata/state_functions`` contains the files ``density.py``, ``temperature.py``, and ``velocity.py``.  Each of these files contains a class which, after instantiation, can be called to return the value of that paricular parameter.  The state functions that describe different structures or features that can be selected for the artificial ionosphere are methods of this class.  As an example, ``density.py`` contains the class ``Density``, which has the methods ``chapman``, ``circle_patch``, and ``wave`` (amongst others) which can be used to create a Chapman layer, a circular density enhancement, or a wave-like structure in plasma density.  Note that there is only one ``Temperature`` class which is used to specify both ion and electron temperature as the functional forms of both are often similar.

New functions should be added as methods to the class of the parameter they represent.  For example, if an ionosphere with a localized ion heating even is desired, a function named something like ``plasma_heating`` should be added to the ``Temperature`` class (in ``temperature.py``) that returns enhanced temperature in a small area.  It technically does not matter where in the class this function is added, but at the bottom or grouped with other similar functions is recommended.

.. _Function Guidelines:

Function Guidelines
*******************

Functions must have a unique name from functions already used by that class, and MUST take in the following variable in this order:

  * Unix Time (seconds since January 1, 1970)
  * Geodetic Latitude (degrees)
  * Geodetic Longitude (degrees)
  * Geodetic Altitude (meters)

Latitude, longitude, and altitude can be scalar or multidimensional arrays, but they will always all be the same shape.  Unix Time will either be a scalar or a 1D array of any length.  Functions should be able to handle any combination of these inputs gracefully.  Functions should return scalar parameters (density and temperature) as arrays where the first dimension matches the length of the Unix Time array, and subsequent dimensions match the shape of the position arrays.  Velocity functions should return a similar array, but with one final dimension that holds the three velocity vector components.

Density values must be returned in units of m-3, temperature values in units of K, and velocity values in units of m/s.  The geodetic components of velocity vectors must be returned in the order east, north, up.  Any transformation to/from other units or coordinate systems must be handled completely within the function.  Parameters input through the configuration file can technically be in any units, but you are encouraged to use standard SI units (meters, seconds, and Kelvin) to be consistent with other functions.

Functions contributed back to the main package MUST include a `numpy-style docstring <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_ that provides a short description of the function and what it should be used for and lists the parameters that must be specified in the configuration file.  ANY parameter your function needs as input aside from the four position/time parameters listed above MUST be in the configuration file.  

Prameters specified in the config file will be available as attributes of the class, so if a parameter `N1` is defined in the config file it will be available for your function to use as `self.N1`.  When multiple state functions are summed for one parameter, an independent insantce of the class is created for each of them, so it is fine to repeat variable names between different functions.

Step-by-step Instructions
*************************

  1. Clone the source repository and install following the :ref:`developer instructions <developer installation>`.
  2. Create a new branch to develop your new function.
  3. Add your function to the appropriate file in ``amisrsynthdata/state_functions`` following the :ref:`Function Guidelines`.
  4. Modify your configuration file to use your new function.
  5. TEST!!!
  6. Push your new branch to GitHub and create a pull request to merge it into the develop branch of the main package repository.
  7. If you are not interested in contributing your function to the main repository, simply merge the branch you developed the function on with your local main branch.

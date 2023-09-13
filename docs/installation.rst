.. installation.rst

Installation
============

The amisrsynthdata package is pure python and can be installed easily with pip.

.. code-block::

  pip install git+https://github.com/amisr/amisrsynthdata.git

.. _developer installation:

Alternatively, if you would like to develop on amisrsythdata, first clone the repo, then install.

.. code-block::

  git clone https://github.com/amisr/amisrsynthdata.git
  cd amisrsythdata
  pip install -e .


Requirements
------------
The following additional packages are required and will be installed:

  * `numpy <https://numpy.org/>`_
  * `h5py <https://docs.h5py.org/en/stable/index.html>`_
  * `pymap3d <https://pypi.org/project/pymap3d/>`_
  * `apexpy <https://apexpy.readthedocs.io/en/latest/>`_
  * `pyyaml <https://pyyaml.org/wiki/PyYAMLDocumentation>`_
  * `importlib-resources <https://pypi.org/project/importlib-resources/>`_

If automatic installation of any of these packages fail, try installing it in your environment manually referring to instructions and tips from that package's documentation.  Additionally, `matplotlib <https://matplotlib.org/>`_ and `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ are required to create the optional summary plots.  Utilizing output from the GEMINI non-linear ionospheric dynamics model requires `pygemini <https://github.com/gemini3d/pygemini>`_ be installed.

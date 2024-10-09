.. installation.rst

Installation
============

The amisrsynthdata package is pure python and can be installed easily with pip.

.. code-block::

  pip install amisrsynthdata

If you would like create the optional summary plots, you will need to install
with the ``plots`` option.

.. code-block::

   pip install amisrsynthdata[plots]

If you would like to use the GEMINI mode to specify the ionosphere, you will
need to install with the ``gemini`` option.

.. code-block::

   pip install amisrsynthdata[gemini]

The amisrsynthdata package depends on `apexpy <https://apexpy.readthedocs.io/en/latest/>`_, which should be installed automatically.  Historically, this package has sometimes been challenging to install on some computer systems.  If the installation of amisrsynthdata fails because apexpy fails to install, or amisrsynthdata crashed due to an apexpy error, consider manually installing apexpy first following the `apexpy installation instructions <https://apexpy.readthedocs.io/en/latest/installation.html>`_, then reinstalling amisrsynthdata as shown above.

.. _developer installation:

Development
-----------

To install amisrsynthdata for development (adding features or fixing bugs):

1. `Fork amisrsynthdata on GitHub <https://github.com/amisr/amisrsynthdata/fork>`_.
2. Clone your fork locally

.. code-block::

    git clone git@github.com:your_name_here/amisrsynthdata.git

3. Install with the editable flag

.. code-block::

  cd amisrsynthdata
  pip install -e .

Note that the installation options discussed above can be used here as well.


Requirements
------------
The following additional packages are required and will be installed:

  * `numpy <https://numpy.org/>`_
  * `h5py <https://docs.h5py.org/en/stable/index.html>`_
  * `pymap3d <https://pypi.org/project/pymap3d/>`_
  * `apexpy <https://apexpy.readthedocs.io/en/latest/>`_
  * `pyyaml <https://pyyaml.org/wiki/PyYAMLDocumentation>`_
  * `importlib-resources <https://pypi.org/project/importlib-resources/>`_ (for python < 3.9)

If automatic installation of any of these packages fail, try installing it in your environment manually referring to instructions and tips from that package's documentation.  Additionally, `matplotlib <https://matplotlib.org/>`_ and `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ are required to create the optional summary plots.  Utilizing output from the GEMINI non-linear ionospheric dynamics model requires `pygemini <https://github.com/gemini3d/pygemini>`_ be installed.

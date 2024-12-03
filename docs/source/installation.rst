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

If installing with the ``plots`` option, the following additional packages are also required and will be installed:


  * `matplotlib <https://matplotlib.org/>`_
  * `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_

If automatic installation of any of these packages fail, try installing them in your environment manually referring to instructions and tips from that package's documentation.

Utilizing output from the GEMINI non-linear ionospheric dynamics model to specify the ionospheric state requires `pygemini <https://github.com/gemini3d/pygemini>`_ to be installed. This is a specialized option and interested users should install this package manually using instructions in its README file before attempting to run `amisrsynthdata` with GEMINI configuration option.


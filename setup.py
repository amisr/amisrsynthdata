############################################################
# This setup.py file ONLY exists to specify dependencies.
# Normally, these should be specified in pyproject.toml, 
# HOWEVER:
#   1. We use autodocs to generate the API documentation.
#   2. ReadTheDocs has to install the package with pip in
#      order to generate autodocs.
#   3. The apexpy dependency cannot install on RTD due to
#      requiring a FORTRAN compiler.
#   4. When the apexpy installation fails, the entire pip 
#      install of amisrsynthdada fails
#   5. There is not presently a flag or option for pip that
#      allows us to exclude the installation of apexpy
#   6. The pyproject.toml standard also does not allow 
#      dynamically selecting dependencies based on the
#      environment.
# Mocking is used for apexpy in RTDs, but if none of it 
# works at all when errors occur during installation.  If
# a more elegant solution to any of these problems becomes
# available, remove this file and go back to all 
# dependencies specified in pyproject.toml.
############################################################

import setuptools
import os

install_requires = [
    "numpy",
    "h5py",
    "apexpy",
    "pymap3d",
    "pyyaml",
    "importlib_resources; python_version < '3.9'",
]

if os.getenv('READTHEDOCS'):
    install_requires.remove("apexpy")

setuptools.setup(install_requires=install_requires)

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
    install_requires.pop("apexpy")

setuptools.setup(install_requires=install_requires)

#!/usr/bin/env python
"""
amisrsynthdata creates synthetic processed AMISR data files (SRI format)
  based on a specifiec ionosphere and radar configuration.  It is useful
  for evaluating high-level data products and inversion techniques.
The full license can be found in LICENSE.txt
"""

import os
import re
import sys
import subprocess
from setuptools import find_packages, setup

here = os.path.abspath(os.path.dirname(__file__))

# # Get the package requirements
# REQSFILE = os.path.join(here, 'requirements.txt')
# with open(REQSFILE, 'r') as f:
#     REQUIREMENTS = f.readlines()
# REQUIREMENTS = '\n'.join(REQUIREMENTS)
REQUIREMENTS = [
    'numpy',
    'h5py',
    'pymap3d',
    'apexpy'
]

# Get the readme text
# README = os.path.join(here, 'README.rst')
with open('README.md', 'r', encoding='utf-8') as f:
    READMETXT = f.read()
# READMETXT = '\n'.join(READMETXT)

# READMETXT = 'Basic readme\n'

# # Get version number from __init__.py
# regex = "(?<=__version__..\s)\S+"
# with open(os.path.join(here,'resolvedvelocities/__init__.py'),'r', encoding='utf-8') as f:
#     text = f.read()
# match = re.findall(regex,text)
# version = match[0].strip("'")

# Package description
DESC = "Tool for generating synthetic AMISR data files."

#############################################################################
# First, check to make sure we are executing
# 'python setup.py install' from the same directory
# as setup.py (root directory)
#############################################################################
# PATH = os.getcwd()
# assert('setup.py' in os.listdir(PATH)), \
#        "You must execute 'python setup.py install' from within the \
# repo root directory."


#############################################################################
# Now execute the setup
#############################################################################
setup(name='amisrsynthdata',
      install_requires=REQUIREMENTS,
      # setup_requires=REQUIREMENTS,
      version='0.0.1',
      description=DESC,
      author="AMISR",
      author_email="leslie.lamarche@sri.com",
      url="https://github.com/amisr/amisrsynthdata",
      download_url="https://github.com/amisr/amisrsynthdata",
      packages=find_packages(),
      long_description=READMETXT,
      zip_safe=False,
      py_modules=['resolvedvelocities'],
      classifiers=["Development Status :: 2.0.0 - Release",
                   "Topic :: Scientific/Engineering",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Natural Language :: English",
                   "Programming Language :: Python",
                  ],
      entry_points={
          'console_scripts': [
              'amisrsynthdata=amisrsynthdata.create_synth_data:main',
        ],
}
      )
